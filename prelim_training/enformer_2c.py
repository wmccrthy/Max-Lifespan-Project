import torch, numpy as np
from enformer_pytorch import from_pretrained
from enformer_pytorch.finetune import HeadAdapterWrapper
from torch import nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset
from torch import Tensor
import pytorch_lightning as L
from pytorch_lightning.loggers import WandbLogger
import wandb, os
from torch.nn import Dropout
import re
import simplejson as json

# dataset class for bigger dataset
# read from tokenized csv

def loadTokenized(num_data=829124):
    data_path = "/data/rsg/chemistry/maggiejl/P2/prelim_training/data/tokenized_enformer_829124_training.csv"
    data = []
    labels = []
    ind = 0
    non_line = 0
    pattern = re.compile(r'[^,]+,[^,]+,([^,]+),"(.*)"')
    with open(data_path) as to_read:
        for line in to_read:
            ind += 1
            if len(data) >= num_data: return data, labels
            if len(data) % 1000 == 0: print("yay for ", len(data), "data points")
            match = pattern.match(line)
            if not match:
                non_line += 1
                continue
            lifespan, seq = match.groups()
            labels.append(float(lifespan))
            data.append(seq)
    print("finished loading in!")
    print("total data: ", ind)
    if non_line > 0:
        print("did not match: ", non_line)
    return data, labels

class EnformerDataset(Dataset):
    """Dataset for supervised fine-tuning."""
    def __init__(self, seq, labels):
        super(EnformerDataset, self).__init__()
        self.labels = labels
        self.seq = seq

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        seq = self.seq[idx].strip()
        seq = json.loads(seq)
        label = self.labels[idx]
        return torch.tensor(seq, dtype=torch.int64), torch.Tensor([label])

class Enformer_Model(L.LightningModule):
    def __init__(self, batch_size, max_len) -> None:
        super().__init__()
        en_pretrain = from_pretrained("EleutherAI/enformer-official-rough", 
                                                  #target_length=max_len,
                                                   use_tf_gamma = False,
        )
        self.model = HeadAdapterWrapper(
            enformer = en_pretrain,
            num_tracks = 128,
            post_transformer_embed = False,
        )
        self.fc1 = nn.Linear(128, 1)
        self.dropout1 = Dropout(p=0.1)
        self.fc2 = nn.Linear(896, 1)
        self.droupout2 = Dropout(p=0.1)

    def forward(self, inputs):
        #print("inputs shape:", inputs.shape) #dim = (batch size, max_len) #usually (batch size, 1, max_len)
        outputs = self.model(inputs) #embed first seq in batch
        #print("outputs shape:", outputs.shape) # dim = (batch size, 896, 128) # usually (batch size, max_len, 768)

        output1 = self.fc1(outputs) #dim = (batch size, 896, 1)
        output1 = self.dropout1(output1)
        output1 = output1.squeeze(-1)
        final_output = self.fc2(output1) #dim = (batch size, 1)
        final_output = self.droupout2(final_output)
        #print("final output shape:", final_output.shape)
        return final_output 

    def training_step(self, batch, batch_idx): 
        # breakpoint()
        loss = self._shared_eval_step(batch, batch_idx)
        self.log("train_loss", loss, sync_dist=True) 
        return loss 

    def validation_step(self, batch, batch_idx):
        loss = self._shared_eval_step(batch, batch_idx)
        self.log("val_loss", loss, sync_dist=True)

    def _shared_eval_step(self, batch, batch_idx):
        inputs, target = batch
        output = self.forward(inputs)  
        #print("eval input shape:",inputs.shape, "target shape :",target.shape,"output:",output.shape)
        #print("new target shape: ", reshaped_target.shape)

        criterion = F.mse_loss(output, target)
        # print("criterion:", criterion)
        loss = torch.sqrt(criterion)
        # print("loss:", loss) 
        return loss

    def configure_optimizers(self):
        optimizer = torch.optim.SGD(self.parameters(), lr=1e-4, weight_decay=1e-5)
        return {
            "optimizer": optimizer,
            "lr_scheduler": {
                "scheduler": torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer),
                "monitor": "val_loss",
                "frequency": 1
            }
        }


if __name__ == "__main__":
    batch_size, num_epochs, num_dev = int(2), int(70), int(1)
    max_len = max(22785, 196608)

    training_inps, training_labels = loadTokenized()
    print("prepared data!")
    print("length of data: ", len(training_inps), len(training_labels))
    
    #avoid hugging face error 
    os.environ["TOKENIZERS_PARALLELISM"] = "false"      

    #split data into training and validation set
    indices = np.arange(len(training_labels))
    np.random.shuffle(indices)
    shuffled_training_inps = [training_inps[i] for i in indices]
    shuffled_training_labels = [training_labels[i] for i in indices]
    
    train_data_size = int(len(training_labels) * .9)
    train_dataset = shuffled_training_inps[:train_data_size]
    train_labels = shuffled_training_labels[:train_data_size]
    valid_dataset = shuffled_training_inps[train_data_size:]
    valid_labels = shuffled_training_labels[train_data_size:]
    
    train_dataset = EnformerDataset(train_dataset, train_labels)
    valid_dataset = EnformerDataset(valid_dataset, valid_labels)
    training_data = DataLoader(train_dataset, batch_size=batch_size, num_workers=1, drop_last = True)
    valid_data = DataLoader(valid_dataset, batch_size=batch_size, num_workers=1, drop_last = True)
    
    print("initializing model")
    model = Enformer_Model(batch_size, max_len)

    # add wandb logger 
    os.environ["WANDB_DIR"]=os.path.abspath("/data/rsg/chemistry/maggiejl/wandb_dir")
    wandb_logger = WandbLogger(log_model="all", project="training-1", dir = "/data/rsg/chemistry/maggiejl/P2/prelim_training")

    print("ready to train!")
    trainer = L.Trainer(max_epochs = num_epochs,
                        accelerator='gpu',
                        devices=num_dev,
                        logger=wandb_logger,
                        ) 
    trainer.fit(model, training_data, valid_data)