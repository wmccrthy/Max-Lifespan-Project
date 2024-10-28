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
from scipy.stats import kurtosis, skew
import csv

"""
This script is informed by EDA we ran to trim cumulative dataset s.t we only have one2one sequences of lengths between 104 and 31620
trimmed cumulative data was then organized into sets per gene
we then ran further analysis to determine the number of sequences, number of species, and max sequence length per gene set

Script for running enformer model on a per gene basis (5 fold validation for each gene, where we read in/load training data in separately)

PSUEDO/STRUCTURE:
enformer model code (repurposed from Maggie's work w/ tweaked params)

code to create map of gene:(num_seqs, num_species) so we only train if gene has num_species > 300

for each gene_one2one_file:
    tokenize as required for transformer (NEED METHOD FOR THIS -- can take from earlier enformer model)

    train

    save results


TO DO:
    - adjust paths from Maggie's dirs to mine
    - trace thru script make sure everything makes sense
"""

gene_one2one_datasets_path = "/data/rbg/users/wmccrthy/chemistry/Everything/gene_datasets/regulatory_one2one/"
gene_one2one_metadata_path = "/data/rbg/users/wmccrthy/chemistry/Everything/EDA/regulatory_one2one_sets_metadata.csv"


def map_gene_to_metadata():
    """
    iterates through regulatory_one2one_sets_metadata.csv to create dict of gene:(num_data, num_species)
    """
    genes = {}
    with open(gene_one2one_metadata_path) as read_from:
         for line in read_from:
              line = line.split(",")
              gene, num_data, num_species, max_len = line
              genes[gene] = (num_data, num_species, max_len)
    return genes

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
    GIVEN GENE SEQUENCE FILE, RETURNS TRAINING DATA (TOKENIZED SEQS) AND LABELS (LIFESPANS)
    """
    data, labels = [], []
    with open(file_path) as to_read:
        for line in to_read:
            line = line.split(",")
            if len(line) < 8 or line[0] == "organism": continue #skip corrupted entries and first line (header)

            lifespan, seq = line[1], line[-1].strip().replace('"',"").replace("\n", "") #make sure indices of seq and lifespan are modified according to what training data is passed 
            seq = seq.replace('"', "").replace(" ", "")
            if len(lifespan) > 5 or len(seq) == 0: 
                print("indicates lifespan is unknown")
                continue #indicates lifespan is 'Unknown' 
            seq = transform4Enformer(seq, max_len)
            data.append(seq)
            labels.append(float(lifespan))
    return data, torch.Tensor(labels)
    
def stats(data):
    data = np.array(data)
    mean = np.mean(data)
    std = np.sqrt(np.var(data, ddof=0))
    sk = skew(data, bias=False)
    kurt = kurtosis(data, fisher=True, bias=False)
    return mean, std, sk, kurt

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
        for p in self.model.parameters(): p.requires_grad = False #freeze enformer layer

        self.fc1 = nn.Linear(128, 1)
        # self.dropout1 = Dropout(p=0.1) REMOVING DROPOUT FOR NOW
        self.fc2 = nn.Linear(max_len, 1)
        # self.droupout2 = Dropout(p=0.1) REMOVING DROPOUT FOR NOW

    def forward(self, inputs):
        #print("inputs shape:", inputs.shape) #dim = (batch size, max_len) #usually (batch size, 1, max_len)
        outputs = self.model(inputs) #embed first seq in batch
        #print("outputs shape:", outputs.shape) # dim = (batch size, 896, 128) # usually (batch size, max_len, 768)

        output1 = self.fc1(outputs) #dim = (batch size, 896, 1)

        # output1 = self.dropout1(output1) REMOVING DROPOUTS FOR NOW

        output1 = output1.squeeze(-1)
        final_output = self.fc2(output1) #dim = (batch size, 1)

        # final_output = self.droupout2(final_output) REMOVING DROPOUTS FOR NOW

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
        # EXPERIMENT W ADAM and SGD
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
    def __init__(self, seqs, labels):
        super().__init__()
        self.labels = labels #list of floats
        self.seq = seqs #numpy array of int64s

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        seq = self.seq[idx]
        label = self.labels[idx]
        return torch.tensor(seq, dtype=torch.int64), torch.Tensor([label])

def train(shuffled_training_inps, shuffled_training_labels, gene, fold, gene_max_len):
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
    training_data = DataLoader(train_dataset, batch_size=batch_size, num_workers=8, drop_last = True)
    valid_data = DataLoader(valid_dataset, batch_size=batch_size, num_workers=8, drop_last = True)
    
    model = Enformer_Model(batch_size, gene_max_len)

    os.environ["WANDB_DIR"]=os.path.abspath("/data/rbg/users/wmccrthy/chemistry/Everything/fall_24_training/wandb_dir")
    wandb_logger = WandbLogger(log_model="all", project="split-by-gene", name = f"{gene}_{len(training_labels)}_fold_{fold}", dir = "/data/rbg/users/wmccrthy/chemistry/Everything/fall_24_training")
    
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
    output_path = "/data/rbg/users/wmccrthy/chemistry/Everything/fall_24_training/stats_enformer_by_gene.csv"
    
    # get mapping of gene:num_data, num_species, max_len
    gene_stats = map_gene_to_metadata()
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
                    
    for gene_path in os.listdir(gene_one2one_datasets_path):
        specific_gene = gene_path.split("_")[0] # extract specific gene from file path
        gene_num_species, gene_max_len = gene_stats[specific_gene][1:]
        print("starting: ", specific_gene)
        torch.cuda.empty_cache()

        # get tokenized data and labels for gene
        training_inps, training_labels = prepare4Enformer(gene_path, gene_max_len)
        print("length of data: ", len(training_inps), len(training_labels))
        mean, std, sk, kurt = stats(training_labels)

        #shuffle dataset
        indices = np.arange(len(training_labels))
        np.random.shuffle(indices)
        shuffled_training_inps = [training_inps[i] for i in indices]
        shuffled_training_labels = [training_labels[i] for i in indices]
        
        train_fold_results = []
        valid_fold_results = []
        #train only if gene set has representation from 300 or more species
        if gene_num_species < 300:
            writeLine = [specific_gene, len(training_labels), mean, std, sk, kurt,]
        else:
            for fold in range(5):
                train_loss, valid_loss, wandbstr = train(shuffled_training_inps, shuffled_training_labels, specific_gene, fold, gene_max_len)
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