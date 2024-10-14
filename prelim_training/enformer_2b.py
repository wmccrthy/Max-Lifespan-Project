import torch, numpy as np, pandas as pd
from enformer_pytorch import Enformer, from_pretrained
from enformer_pytorch.finetune import HeadAdapterWrapper
from torch import nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset, TensorDataset, random_split
from torch import Tensor
import pytorch_lightning as L
from transformers import AutoTokenizer, AutoModel
import sys 
from functools import reduce 
from pytorch_lightning.loggers import WandbLogger
from pytorch_lightning.callbacks import ModelCheckpoint
import wandb, os
from torch.nn import Dropout
import simplejson as json
import re
    
#adjustments from 1b for a big dataset -- does not work
#could try again with longtmux     

#learning rate scheduler 
#optional dropout/weight decay

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
            seq = seq.strip()
            seq = json.loads(seq)
            labels.append(float(lifespan))
            data.append(torch.tensor(seq))
    print("finished loading in!")
    print("total data: ", ind)
    print("did not match: ", non_line)
    return data, labels

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
        print("inputs shape:", inputs.shape) #dim = (batch size, max_len)
        outputs = self.model(inputs) #embed first seq in batch
        print("outputs shape:", outputs.shape) # dim = (batch size, 896, 128) # usually (batch size, max_len, 768)

        output1 = self.fc1(outputs) #dim = (batch size, 896, 1)
        output1 = self.dropout1(output1)
        output1 = output1.squeeze(-1)
        final_output = self.fc2(output1) #dim = (batch size, 1)
        final_output = self.droupout2(final_output)
        print("final output shape:", final_output.shape)
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
        reshaped_target = torch.unsqueeze(target, dim=1)
        #print("new target shape: ", reshaped_target.shape)

        criterion = F.mse_loss(output, reshaped_target)
        # print("criterion:", criterion)
        loss = torch.sqrt(criterion)
        # print("loss:", loss) 
        return loss

    def configure_optimizers(self):
        optimizer = torch.optim.Adam(self.parameters(), lr=1e-4, weight_decay=1e-5)
        return {
            "optimizer": optimizer,
            "lr_scheduler": {
                "scheduler": torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer),
                "monitor": "val_loss",
                "frequency": 1
            }
        }

def train():
    batch_size, num_epochs, num_dev = int(2), int(70), int(1)
    max_len = max(22785, 196608)

    training_inps, training_labels = loadTokenized()
    print("prepared data!")
    print("length of data: ", len(training_inps), len(training_labels))
    print(training_inps[0][:20], type(training_inps[0]))
    #training_inps = [torch.Tensor(np.fromstring(entry[1:-1], sep=',', dtype=np.float32)) for entry in training_inps]
    print(training_inps[0][:20], type(training_inps[0][0]))
    training_inps = torch.stack(training_inps)
    training_labels = torch.Tensor(training_labels)
    print("done with stacking")
    
    #avoid hugging face error 
    os.environ["TOKENIZERS_PARALLELISM"] = "false"      

    training_data = TensorDataset(training_inps, training_labels)
    #split data into training and validation set
    train_data_size = int(len(training_data) * .8)
    valid_data_size = len(training_data) - train_data_size

    training_data, valid_data = random_split(training_data, [train_data_size, valid_data_size])
    training_data = DataLoader(training_data, batch_size=batch_size, num_workers=0, drop_last = True)
    valid_data = DataLoader(valid_data, batch_size=batch_size, num_workers=0, drop_last = True)
    print("loaded data")
    model = Enformer_Model(batch_size, max_len)
    print('initialized model')

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


if __name__ == "__main__":
    args = sys.argv
    globals()[args[1]](*args[2:])