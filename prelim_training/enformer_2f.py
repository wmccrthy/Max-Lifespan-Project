import torch
from enformer_pytorch import from_pretrained
from enformer_pytorch.finetune import HeadAdapterWrapper
from torch import nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset
from torch import Tensor
# import lightning as L
import pytorch_lightning as L
from pytorch_lightning.loggers import WandbLogger
import wandb, os
from torch.nn import Dropout
import numpy as np
from scipy.stats import kurtosis, skew
import re
import csv
import simplejson as json

#enformer_1b.py but on many different genes
# 5 fold cross validation

#frozen enformer layers
#learning rate scheduler 
#optional dropout/weight decay

def stats(data):
    data = np.array(data)
    mean = np.mean(data)
    std = np.sqrt(np.var(data, ddof=0))
    sk = skew(data, bias=False)
    kurt = kurtosis(data, fisher=True, bias=False)
    return mean, std, sk, kurt
    
def loadTokenized(num_data=829124):
    data_path = "/data/rsg/chemistry/maggiejl/P2/prelim_training/data/tokenized_enformer_829124_training.csv"
    data = []
    labels = []
    genes = []
    unique_genes = set()
    ind = 0
    non_line = 0
    pattern = re.compile(r'([^,]+),[^,]+,([^,]+),"(.*)"')
    with open(data_path) as to_read:
        for line in to_read:
            ind += 1
            if len(data) >= num_data: return data, labels, genes, unique_genes
            if len(data) % 1000 == 0: print("datapoint ", len(data), " done")
            match = pattern.match(line)
            if not match:
                non_line += 1
                continue
            gene, lifespan, seq = match.groups()
            seq = seq.strip()
            #seq = ast.literal_eval(seq)
            labels.append(float(lifespan))
            data.append(seq)
            genes.append(gene)
            unique_genes.add(gene)
            #print(seq[:30], type(seq[0]))
    print("finished loading in!")
    print("total data: ", ind)
    print("did not match: ", non_line)
    return data, labels, genes, unique_genes

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
        for param in self.model.enformer.parameters():
            param.requires_grad = False
        self.fc1 = nn.Linear(128, 1)
        self.dropout1 = Dropout(p=0.1)
        self.fc2 = nn.Linear(896, 1)
        self.droupout2 = Dropout(p=0.1)

    def forward(self, inputs):
        #print("inputs shape:", inputs.shape) #dim = (batch size, max_len) #usually (batch size, 1, max_len)
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
        optimizer = torch.optim.Adam(self.parameters(), lr=1e-4, weight_decay=1e-5)
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

def train(shuffled_training_inps, shuffled_training_labels, gene, fold):
    batch_size, num_epochs, num_dev = int(2), int(40), int(1)
    
    #split data into training and validation set
    fold_size = len(shuffled_training_labels) // 5
    val_start = fold*fold_size
    val_end = (fold+1)*fold_size
    print("indices for train: 0, ", val_start, ", and ", val_end, "to end")
    print("indices for val: ", val_start, ", ", val_end)
    train_dataset = shuffled_training_inps[:val_start] + shuffled_training_inps[val_end:]
    train_labels = shuffled_training_labels[:val_start] + shuffled_training_labels[val_end:]
    valid_dataset = shuffled_training_inps[val_start:val_end]
    valid_labels = shuffled_training_labels[val_start:val_end]
    
    train_dataset = EnformerDataset(train_dataset, train_labels)
    valid_dataset = EnformerDataset(valid_dataset, valid_labels)
    training_data = DataLoader(train_dataset, batch_size=batch_size, num_workers=0, drop_last = True)
    valid_data = DataLoader(valid_dataset, batch_size=batch_size, num_workers=0, drop_last = True)
    
    model = Enformer_Model(batch_size, 196608)

    os.environ["WANDB_DIR"]=os.path.abspath("/data/rsg/chemistry/maggiejl/wandb_dir")
    wandb_logger = WandbLogger(log_model="all", project="split-by-gene", name = f"{gene}_{len(training_labels)}_fold_{fold}", dir = "/data/rsg/chemistry/maggiejl/P2/prelim_training")
    
    print("ready to train!")
    trainer = L.Trainer(max_epochs = num_epochs,
                        accelerator='gpu',
                        devices=num_dev,
                        logger=wandb_logger,
                        )
    trainer.fit(model, training_data, valid_data)
    
    api = wandb.Api()
    wandb_id = wandb_logger.experiment.id
    wandb_run_link = f"split-by-gene/{wandb_id}"
    print(f"wandb run: {wandb_run_link}")
    run = api.run(f"split-by-gene/{wandb_id}")
    run = wandb.Api().run(f"{wandb_logger.experiment.project}/{wandb_id}")
    print(run.history().columns)
    try:
        train_loss = run.history(keys=["train_loss"]).tail(1).iloc[0]['train_loss']
    except:
        train_loss = "N/A"
    valid_loss = run.history(keys=["val_loss"]).tail(1).iloc[0]['val_loss']
    wandb.finish()
    return train_loss, valid_loss, wandb_run_link

if __name__ == "__main__":
    os.environ["TOKENIZERS_PARALLELISM"] = "false"  
    output_path = "/data/rsg/chemistry/maggiejl/P2/prelim_training/data/stats_enformer_bygene.csv"
    data, labels, genes, unique_genes = loadTokenized()
    
    #clear csv
    with open(output_path, mode='w', newline='') as file: pass
    #add first row
    with open(output_path, "a", newline='') as write_to:
                    writer = csv.writer(write_to)
                    writer.writerow(["gene", "datapoints",
                                     "mean", "std", "skew", "kurtosis",
                                     "wandb", "train_f1", "train_f2", "train_f3", "train_f4",
                                     "train_f5", "valid_f1", "valid_f2", "valid_f3", "valid_f4",
                                     "valid_f5", "train_avg", "valid_avg"])
                    
    for specific_gene in sorted(list(unique_genes)):
        torch.cuda.empty_cache()
        print("starting: ", specific_gene)
        # extract data and labels that match the gene
        ind = [i for i, gene in enumerate(genes) if gene==specific_gene]
        training_inps = [data[i] for i in ind]
        training_labels = [labels[i] for i in ind]
        print("length of data: ", len(training_inps), len(training_labels))
        mean, std, sk, kurt = stats(training_labels)
        
        #shuffle dataset
        indices = np.arange(len(training_labels))
        np.random.shuffle(indices)
        shuffled_training_inps = [training_inps[i] for i in indices]
        shuffled_training_labels = [training_labels[i] for i in indices]
        
        #train only if datapoints more than 100
        train_fold_results = []
        valid_fold_results = []
        if len(training_labels) < 100:
            wandbstr, train_fold_results, valid_fold_results, train_avg, valid_avg = "N/A", "N/A", "N/A", "N/A", "N/A"
        else:
            for fold in range(5):
                train_loss, valid_loss, wandbstr = train(shuffled_training_inps, shuffled_training_labels, specific_gene, fold)
                train_fold_results.append(train_loss)
                valid_fold_results.append(valid_loss)
            train_avg = sum(train_fold_results) / len(train_fold_results)
            valid_avg = sum(valid_fold_results) / len(valid_fold_results)
        
        print("writing line")
        
        writeLine = [specific_gene, len(training_labels), 
                     mean, std, sk, kurt,
                     wandbstr] + train_fold_results + valid_fold_results + [train_avg, valid_avg]
        with open(output_path, "a", newline='') as write_to:
                    writer = csv.writer(write_to)
                    writer.writerow(writeLine)
    
    