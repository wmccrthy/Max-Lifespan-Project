import torch, numpy as np
from enformer_pytorch import from_pretrained
from enformer_pytorch.finetune import HeadAdapterWrapper
from torch import nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset, random_split
from torch import Tensor
import pytorch_lightning as L
from pytorch_lightning.loggers import WandbLogger
import wandb, os
from torch.nn import Dropout
from pytorch_lightning.callbacks import EarlyStopping
import h5py
import csv

# dataset class for bigger dataset
# read from h5py file

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

class EnformerDataset(Dataset):
    """Dataset for supervised fine-tuning."""
    def __init__(self, seqs, labels):
        super().__init__()
        self.labels = labels #list of floats
        self.seq = seqs #numpy array of int64s

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        seq = self.seq[idx, :]
        label = self.labels[idx]
        return torch.tensor(seq, dtype=torch.int64), torch.Tensor([label])

class EnformerDataModule(L.LightningDataModule):
    def __init__(self, batch_size, max_data=829124):
        super().__init__()
        self.batch_size = int(batch_size)
        self.max_data = int(max_data)
        self.train_dataset = None
        self.valid_dataset = None
        
    #def prepare_data(self):
        h5py_path = "/data/rsg/chemistry/maggiejl/P2/prelim_training/data/tokenized_enformer_829124_training.h5"        
        meta_path = "/data/rsg/chemistry/maggiejl/P2/prelim_training/data/arbitrary_one2one_829124_training.csv"
        
        print("opening h5py file")
        hf = h5py.File(h5py_path, 'r')
        datapoints = hf.get('dataset')
        
        lifespans = []
        with open(meta_path) as read_from:
            for line in read_from:
                if len(lifespans) >= self.max_data: break
                line = line.split(",")
                lifespans.append(float(line[1])) 
        
        data_size = datapoints.shape[0]
        assert len(lifespans) == data_size
        print("size of data: ", data_size)
        train_size = int(0.8*data_size)
        val_size = data_size - train_size
        
        comb_dataset = EnformerDataset(seqs = datapoints, labels = lifespans)
        self.train_dataset, self.valid_dataset = random_split(comb_dataset, [train_size, val_size])
    
    def train_dataloader(self):
        return DataLoader(self.train_dataset, batch_size=self.batch_size, num_workers=80, drop_last = True)
    
    def val_dataloader(self):
        return DataLoader(self.valid_dataset, batch_size=self.batch_size, num_workers=80, drop_last = True)

if __name__ == "__main__":
    batch_size, num_epochs, num_dev = int(2), int(70), int(4)
    max_data = 829124
    
    print("initializing model")
    model = Enformer_Model(batch_size, max_len = 196608)
    data_module = EnformerDataModule(batch_size = batch_size, max_data = max_data)

    os.environ["WANDB_DIR"]=os.path.abspath("/data/rsg/chemistry/maggiejl/wandb_dir")
    wandb_logger = WandbLogger(log_model="all", project="training-1", dir = "/data/rsg/chemistry/maggiejl/P2/prelim_training")
    early_stopping = EarlyStopping(monitor='val_loss', min_delta=0.001, patience=5, verbose=True, mode='min')
    print("ready to train!")
    trainer = L.Trainer(max_epochs = num_epochs,
                        accelerator='gpu',
                        devices=num_dev,
                        logger=wandb_logger,
                        callbacks=[early_stopping],
                        strategy='ddp_find_unused_parameters_true',
                        ) 
    trainer.fit(model, datamodule=data_module)