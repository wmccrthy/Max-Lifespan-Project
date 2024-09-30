import torch, numpy as np, pandas as pd
from enformer_pytorch import Enformer, from_pretrained
from enformer_pytorch.finetune import HeadAdapterWrapper
from torch import nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset, TensorDataset, random_split
from torch import Tensor
# import lightning as L
import pytorch_lightning as L
from transformers import AutoTokenizer, AutoModel
import sys 
from functools import reduce 
from pytorch_lightning.loggers import WandbLogger
from pytorch_lightning.callbacks import ModelCheckpoint
import wandb, os
from torch.nn import Dropout
           
#learning rate scheduler 
#optional dropout/weight decay

def get_max_len(file_path, passed_in_max_len=None):
    if passed_in_max_len is not None:
        return passed_in_max_len
    max_len = 0
    with open(file_path) as to_read:
        for line in to_read:
            line = line.split(",")
            if len(line) < 2: continue 
            lifespan, seq = line[1], line[2].strip().replace('"',"").replace("\n", "") #make sure indices of seq and lifespan are modified according to what training data is passed 
            if len(lifespan) > 5 or len(seq) == 0: continue #indicates lifespan is 'Unknown' 
            max_len = max(max_len, len(seq))
    print("Pad up to:", max_len)
    return min(max_len, 196608) #max length of enformer model
 
def transform4Enformer(seq, max_len):
    """
    usually in format TCAACATGGAGT
    needs to be in format torch.randint(0, 5, (1, max_len)) for ACGTN
    -1 for padding
    """
    if len(seq) > max_len: seq = seq[:max_len]
    newseq = []
    for letter in seq:
        if letter == "A": newseq.append(0)
        elif letter == "C": newseq.append(1)
        elif letter == "G": newseq.append(2)
        elif letter == "T": newseq.append(3)
        elif letter == "N": newseq.append(4)
    newseq = newseq + [-1]*(max_len - len(newseq))
    newseq = torch.tensor(newseq)
    return newseq

def prepare4Enformer(file_path, max_len):
    """
    GIVEN SEQUENCE FILE, RETURNS MAX TOKENIZED LENGTH OF SEQ FROM THAT FILE S.T WE KNOW WHERE TO PAD TO
    """
    data, labels = [], []
    with open(file_path) as to_read:
        for line in to_read:
            line = line.split(",")
            if len(line) < 2: continue 
            lifespan, seq = line[1], line[-1].strip().replace('"',"").replace("\n", "") #make sure indices of seq and lifespan are modified according to what training data is passed 
            seq = seq.replace('"', "").replace(" ", "")
            if len(lifespan) > 5 or len(seq) == 0: 
                print("indicates lifespan is unknown")
                continue #indicates lifespan is 'Unknown' 
            seq = transform4Enformer(seq, max_len)
            data.append(seq)
            labels.append(float(lifespan))
    return data, torch.Tensor(labels)

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
        #print("inputs type: ", inputs.dtype)
        outputs = self.model(inputs) #embed first seq in batch
        #print("outputs shape:", outputs.shape) # dim = (batch size, 896, 128) # usually (batch size, max_len, 768)

        output1 = self.fc1(outputs) #dim = (batch size, 896, 1)
        #output1 = self.dropout1(output1)
        output1 = output1.squeeze(-1)
        final_output = self.fc2(output1) #dim = (batch size, 1)
        #final_output = self.droupout2(final_output)
        #print("final output shape:", final_output.shape)
        return final_output 

    def training_step(self, batch, batch_idx): 
        # breakpoint()
        loss = self._shared_eval_step(batch, batch_idx)
        self.log("train_loss", loss, sync_dist=True) 
        return loss 

    def validation_step(self, batch, batch_idx):
        loss = self._temp_shared_eval_step(batch, batch_idx)
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
    
    def _temp_shared_eval_step(self, batch, batch_idx):
        inputs, target = batch
        output = self.forward(inputs)
        #print("eval input shape:",inputs.shape, "target shape :",target.shape,"output:",output.shape)
        reshaped_target = torch.unsqueeze(target, dim=1)
        print(reshaped_target[:5])
        #print("new target shape: ", reshaped_target.shape)

        criterion = F.mse_loss(output, reshaped_target)
        # print("criterion:", criterion)
        loss = torch.sqrt(criterion)
        # print("loss:", loss) 
        return loss

    def configure_optimizers(self):
        optimizer = torch.optim.Adam(self.parameters(), lr=1e-4) #, weight_decay=1e-5)
        return {
            "optimizer": optimizer,
            "lr_scheduler": {
                "scheduler": torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer),
                "monitor": "val_loss",
                "frequency": 1
            }
        }


def train():
    training_data = "/data/rsg/chemistry/maggiejl/P2/prelim_training/data/arbitrary_1000_training.csv"
    batch_size, num_epochs, num_dev = int(2), int(30), int(1)

    max_len = get_max_len(training_data, 196608)
    print("got max len: ", max_len)
    training_inps, training_labels = prepare4Enformer(training_data, int(max_len))
    print("prepared data!")
    print("length of data: ", len(training_inps), len(training_labels))
    training_inps = torch.stack(training_inps)
    
    #avoid hugging face error 
    os.environ["TOKENIZERS_PARALLELISM"] = "false"      

    training_data = TensorDataset(training_inps, training_labels)
    #split data into training and validation set
    train_data_size = int(len(training_data) * .8)
    valid_data_size = len(training_data) - train_data_size

    training_data, valid_data = random_split(training_data, [train_data_size, valid_data_size])
    training_data = DataLoader(training_data, batch_size=batch_size, num_workers=0, drop_last = True)
    valid_data = DataLoader(valid_data, batch_size=batch_size, num_workers=0, drop_last = True)
    
    model = Enformer_Model(batch_size, max_len)

    # add wandb logger 
    os.environ["WANDB_DIR"]=os.path.abspath("/data/rsg/chemistry/maggiejl/wandb_dir")
    wandb_logger = WandbLogger(log_model="all", project="training-1", dir = "/data/rsg/chemistry/maggiejl/P2/prelim_training")
    #wandb_logger.download(root = "/data/rsg/chemistry/maggiejl/P2/prelim_training", allow_missing_references = True, skip_cache = True)
    # cback = ModelCheckpoint(monitor="val_loss", mode="min", save_top_k=1)

    print("ready to train!")
    trainer = L.Trainer(max_epochs = num_epochs,
                        accelerator='gpu',
                        devices=num_dev,
                        logger=wandb_logger,
                        ) 
    trainer.fit(model, training_data, valid_data)


if __name__ == "__main__":
    train()