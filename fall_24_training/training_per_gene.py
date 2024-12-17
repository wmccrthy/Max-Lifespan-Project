import torch
import numpy as np
from enformer_pytorch import Enformer, seq_indices_to_one_hot
from enformer_pytorch import from_pretrained
from enformer_pytorch.finetune import HeadAdapterWrapper, get_enformer_embeddings
from torch import nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset, random_split
from torch import Tensor
import pytorch_lightning as L
from pytorch_lightning.loggers import WandbLogger, CSVLogger
import wandb
import os
from torch.nn import Dropout
from pytorch_lightning.callbacks import EarlyStopping
from scipy.stats import kurtosis, skew
import csv
import sys, random

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

gene_one2one_datasets_path = (
    "/data/rbg/users/wmccrthy/chemistry/Everything/gene_datasets/regulatory_one2one/"
)
gene_one2one_metadata_path = "/data/rbg/users/wmccrthy/chemistry/Everything/EDA/regulatory_one2one_sets_metadata.csv"
training_dir = "/data/rbg/users/wmccrthy/chemistry/Everything/fall_24_training/"
UNIVERSAL_MAX_LEN = 31620


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
    if len(seq) > max_len:
        seq = seq[:max_len]
    newseq = []
    for letter in seq:
        if letter == "A":
            newseq.append(0)
        elif letter == "C":
            newseq.append(1)
        elif letter == "G":
            newseq.append(2)
        elif letter == "T":
            newseq.append(3)
        elif letter == "N":
            newseq.append(4)
    newseq = newseq + [-1] * (max_len - len(newseq))
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
            if len(line) < 8 or line[0] == "organism":
                continue  # skip corrupted entries and first line (header)
            # make sure indices of seq and lifespan are modified according to what training data is passed
            lifespan, seq = line[1], line[-1].strip().replace('"', "").replace(
                "\n", ""
            ).replace(" ", "")
            if len(lifespan) > 5 or len(seq) == 0:
                # print("indicates lifespan is unknown")
                continue  # indicates lifespan is 'Unknown'
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


def get_model_params(model):
    # Access and use only the embedding weights -> pass our inputs thru embedding weights
    print("====================== Named Params ======================")
    for name, param in model.named_parameters():
        # if "embedding" in name:
        print(name)
        # if "embedding" in name:  # Adjust based on the exact layer names in Enformer
        #     # Do something with `param`, e.g., copy or modify
        #     embedding_weights = param.detach().clone()
    print("====================== Modules ======================")
    for module in model.modules():
        print(module)


class MLP(nn.Module):
    def __init__(
        self, input_dim, hidden_dim1, hidden_dim2, output_dim, dropout_prob=0.1
    ):
        super(MLP, self).__init__()
        self.fc1 = nn.Linear(input_dim, hidden_dim1)
        self.dropout1 = nn.Dropout(dropout_prob)
        self.fc2 = nn.Linear(hidden_dim1, hidden_dim2)
        self.dropout2 = nn.Dropout(dropout_prob)
        self.fc3 = nn.Linear(hidden_dim2, output_dim)

    def forward(self, x):
        x = F.relu(self.fc1(x))  # First hidden layer
        x = self.dropout1(x)
        x = F.relu(self.fc2(x))  # Second hidden layer
        x = self.dropout2(x)
        x = self.fc3(x)  # Output layer
        return x


class Enformer_Model(L.LightningModule):
    def __init__(self, batch_size, max_len) -> None:
        super().__init__()

        dim = 1536
        target_length = 248  # target length for input seq length of 31620

        self.model = Enformer.from_hparams(
            dim=dim,
            depth=11,
            heads=8,
            output_heads=dict(human=5313, mouse=1643),
            target_length=target_length,
        )
        # for p in self.model.parameters(): p.requires_grad = False #freeze enformer layer

        self.fc1 = nn.Linear(dim * 2, 1)
        self.dropout1 = Dropout(p=0.05)

        # self.fc2 = nn.Linear(target_length, 1)
        # self.droupout2 = Dropout(p=0.1)

        self.mlp = MLP(
            input_dim=target_length,
            hidden_dim1=128,
            hidden_dim2=64,
            output_dim=1,
            dropout_prob=0.05,
        )

    def forward(self, inputs, target):
        # print("inputs shape pre-mod:", inputs.shape) #dim = (batch size, max_len) #usually (batch size, 1, max_len)
        # breakpoint()

        # ====================== PRE-TRAINED ======================
        # output, embeddings = self.model(inputs, target = target) # pass inputs

        # ====================== FROM_HPARAMS ======================
        embeddings = self.model(inputs, return_only_embeddings=True)
        # embeddings = get_enformer_embeddings(self.model, inputs)
        # print(embeddings.shape)

        # pass embeddings thru a linear layer and mlp layer
        output1 = self.fc1(embeddings)
        # print(output1.shape)
        output1 = self.dropout1(output1)
        output1 = torch.permute(output1, (0, 2, 1))
        # print("permuted shape:", output1.shape)

        # final_output = self.droupout2(self.fc2(output1))
        # final_output = self.fc2(output1)
        final_output = self.mlp(output1)
        # print(final_output.shape, final_output)
        return final_output

    def training_step(self, batch, batch_idx):
        # breakpoint()
        loss = self._shared_eval_step(batch, batch_idx)
        # self.log("train_loss", loss, sync_dist=True) ===================== FOR WANDB =====================
        # ===================== FOR CSV =====================
        self.log("train_loss", loss, on_step=False, on_epoch=True)
        return loss

    def validation_step(self, batch, batch_idx):
        loss = self._shared_eval_step(batch, batch_idx)
        # self.log("val_loss", loss, sync_dist=True) ===================== FOR WANDB =====================
        # ===================== FOR CSV =====================
        self.log("val_loss", loss, on_step=False, on_epoch=True)

    def _shared_eval_step(self, batch, batch_idx):
        inputs, target = batch
        # print("inputs: ", inputs, len(inputs[0]), len(inputs[1]), inputs.shape)
        output = self.forward(inputs, target)

        # print("eval input shape:",inputs.shape, "target shape :",target.shape,"output:",output.shape)
        # print("new target shape: ", reshaped_target.shape)

        criterion = F.mse_loss(output, target)

        loss = torch.sqrt(criterion)
        # print("loss:", loss)
        return loss

    def configure_optimizers(self):
        # EXPERIMENT W ADAM and SGD
        optimizer = torch.optim.Adam(self.parameters(), lr=1e-4, weight_decay=1e-5)
        return {
            "optimizer": optimizer,
            "lr_scheduler": {
                "scheduler": torch.optim.lr_scheduler.ReduceLROnPlateau(
                    optimizer, min_lr=1e-6
                ),
                "monitor": "val_loss",
                "frequency": 1,
            },
        }


class EnformerDataset(Dataset):
    """Dataset for supervised fine-tuning."""

    def __init__(self, seqs, labels):
        super().__init__()
        self.labels = labels  # list of floats
        self.seq = seqs  # numpy array of int64s

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        seq = self.seq[idx]
        label = self.labels[idx]
        # return torch.tensor(seq, dtype=torch.int64), torch.Tensor([label]) # seq should already be in tensor form, so try returning diff vals
        return seq, torch.Tensor([label])


def train(shuffled_training_inps, shuffled_training_labels, gene, fold, gene_max_len, tag=None, valid_inps = None):
    batch_size, num_epochs, num_dev = int(1), int(10), int(1)

    # ===================== Implements 5 cross fold validation splitting =====================
    # split data into training and validation set
    # fold_size = len(shuffled_training_labels) // 5
    # val_start = fold*fold_size
    # val_end = (fold+1)*fold_size
    # print("indices for train: 0, ", val_start, ", and ", val_end, "to end")
    # print("indices for val: ", val_start, ", ", val_end)
    # train_dataset = shuffled_training_inps[:val_start] + shuffled_training_inps[val_end:]
    # train_labels = shuffled_training_labels[:val_start] + shuffled_training_labels[val_end:]
    # valid_dataset = shuffled_training_inps[val_start:val_end]
    # valid_labels = shuffled_training_labels[val_start:val_end]

    # train_dataset = EnformerDataset(train_dataset, train_labels)
    # valid_dataset = EnformerDataset(valid_dataset, valid_labels)
    # training_data = DataLoader(train_dataset, batch_size=batch_size, num_workers=8, drop_last = True)
    # valid_data = DataLoader(valid_dataset, batch_size=batch_size, num_workers=8, drop_last = True)

    # ===================== Implements random split =====================
    # split data randomly
    # breakpoint()
    torch.manual_seed(42)
    train_ratio, val_ratio = 0.8, 0.2
    cumulative_dataset = EnformerDataset(
        shuffled_training_inps, shuffled_training_labels
    )
    train_data, val_data = random_split(cumulative_dataset, [train_ratio, val_ratio])
    print(len(train_data), len(val_data))
    training_data = DataLoader(
        train_data, batch_size=batch_size, num_workers=8, drop_last=True
    )
    valid_data = DataLoader(
        val_data, batch_size=batch_size, num_workers=8, drop_last=True
    )
    # re-compile training and validation sets if valid_inps was passed
    if valid_inps:
        training_data = DataLoader(
            cumulative_dataset, batch_size=batch_size, num_workers=8, drop_last=True
        )
        valid_data = DataLoader(
            EnformerDataset(valid_inps[0], valid_inps[1]), # valid inps will be passed as tuple of lists where valid_data[0] is seqs and valid_data[1] is labels
            batch_size=batch_size, num_workers=8, drop_last=True
        )

    model = Enformer_Model(batch_size, gene_max_len)

    os.environ["WANDB_DIR"] = os.path.abspath(
        "/data/rbg/users/wmccrthy/chemistry/Everything/fall_24_training/wandb_dir"
    )

    # ===================== FOR 5-CROSS FOLD =====================
    # wandb_logger = WandbLogger(log_model="all", project="split-by-gene", name = f"{gene}_{len(training_labels)}_fold_{fold}", dir = "/data/rbg/users/wmccrthy/chemistry/Everything/fall_24_training")
    # ===================== FOR RANDOM SPLIT =====================
    # wandb_logger = WandbLogger(log_model="all", project="split-by-gene", name = f"{gene}_{len(training_labels)}_random_split", dir = "/data/rbg/users/wmccrthy/chemistry/Everything/fall_24_training")

    if tag:
        tag = "_" + tag
    else:
        tag = ""
    log_dir = f"/data/rbg/users/wmccrthy/chemistry/Everything/fall_24_training/gene_training_metrics{tag}"
    # check if CSVLogger directory exists
    # Check if the directory exists
    if not os.path.exists(log_dir):
        # Create the directory if it doesn't exist
        os.makedirs(log_dir)
        print(f"Directory '{log_dir}' created.")
    else:
        print(f"Directory '{log_dir}' already exists.")

    # init CSVLogger
    csv_logger = CSVLogger(
        log_dir,
        name=f"{gene}_{len(cumulative_dataset)}_random_split",
    )

    # init EarlyStopping
    early_stopping = EarlyStopping(
        monitor="val_loss",  # metric to monitor
        # number of epochs with no improvement after which training will be stopped
        patience=4,
        verbose=True,
        mode="min",  # "min" because you want to stop when val_loss stops decreasing
    )

    print("ready to train!")
    trainer = L.Trainer(
        max_epochs=num_epochs,
        accelerator="gpu",
        devices=num_dev,
        callbacks=[early_stopping],
        # logger=wandb_logger,
        logger=csv_logger,
    )
    trainer.fit(model, training_data, valid_data)

    # # ===================== FOR WANDB =====================
    # api = wandb.Api()
    # wandb_id = wandb_logger.experiment.id
    # wandb_run_link = f"split-by-gene/{wandb_id}"
    # print(f"wandb run: {wandb_run_link}")
    # run = api.run(f"split-by-gene/{wandb_id}")
    # run = wandb.Api().run(f"{wandb_logger.experiment.project}/{wandb_id}")
    # print(run.history().columns)
    # try:
    #     train_loss = run.history(keys=["train_loss"]).tail(1).iloc[0]['train_loss']
    # except:
    #     train_loss = "N/A"
    # valid_loss = run.history(keys=["val_loss"]).tail(1).iloc[0]['val_loss']
    # wandb.finish()
    # return train_loss, valid_loss, wandb_run_link


"""
Method to train across all genes (individually)
"""
def train_all_genes(tag=None):
    os.environ["TOKENIZERS_PARALLELISM"] = "false"
    # output_path = "/data/rbg/users/wmccrthy/chemistry/Everything/fall_24_training/stats_enformer_by_gene.csv"

    # get mapping of gene:num_data, num_species, max_len
    gene_stats = map_gene_to_metadata()

    # clear csv
    # with open(output_path, mode='w', newline='') as file: pass
    # #add first row
    # with open(output_path, "a", newline='') as write_to:
    #     writer = csv.writer(write_to)
    #     # ===================== FOR 5-CROSS FOLD =====================
    #     # writer.writerow(["gene", "datapoints",
    #     #                 "mean", "std", "skew", "kurtosis",
    #     #                 "wandb", "train_f1", "train_f2", "train_f3", "train_f4",
    #     #                 "train_f5", "valid_f1", "valid_f2", "valid_f3", "valid_f4",
    #     #                 "valid_f5", "train_avg", "valid_avg"])
    #     # ===================== FOR RANDOM SPLIT =====================
    #     writer.writerow(["gene", "datapoints",
    #                     "mean", "std", "skew", "kurtosis",
    #                     "wandb", "train_loss", "valid_loss"])

    for gene_path in os.listdir(gene_one2one_datasets_path):
        # extract specific gene from file path
        specific_gene = gene_path.split("_")[0]
        gene_num_species, gene_max_len = gene_stats[specific_gene][1:]
        gene_num_species, gene_max_len = int(gene_num_species), int(gene_max_len)
        print(
            "starting: ",
            specific_gene,
            " num species: ",
            gene_num_species,
            " max seq len: ",
            gene_max_len,
        )

        # train only if gene set has representation from 300 or more species (ensure we perform this check bfore tokenizing and pre-processing data)
        if gene_num_species < 400:
            print("skipping", specific_gene, "")
            # writeLine = [specific_gene, len(training_labels), mean, std, sk, kurt,]
            continue

        # ====================== FOR CONVENIENCE, PADDING ALL GENES TO SAME LENGTH ======================
        gene_max_len = UNIVERSAL_MAX_LEN

        torch.cuda.empty_cache()

        # get tokenized data and labels for gene
        full_gene_path = gene_one2one_datasets_path + gene_path
        training_inps, training_labels = prepare4Enformer(full_gene_path, gene_max_len)
        print("length of data: ", len(training_inps), len(training_labels))
        # mean, std, sk, kurt = stats(training_labels)

        # shuffle dataset
        indices = np.arange(len(training_labels))
        np.random.shuffle(indices)
        shuffled_training_inps = [training_inps[i] for i in indices]
        shuffled_training_labels = [training_labels[i] for i in indices]

        # train_fold_results = []
        # valid_fold_results = []
        # ===================== FOR 5-CROSS FOLD =====================
        # for fold in range(5):
        #     train_loss, valid_loss, wandbstr = train(shuffled_training_inps, shuffled_training_labels, specific_gene, fold, gene_max_len)
        #     train_fold_results.append(train_loss)
        #     valid_fold_results.append(valid_loss)
        # train_avg = sum(train_fold_results) / len(train_fold_results)
        # valid_avg = sum(valid_fold_results) / len(valid_fold_results)
        # print("writing line")
        # writeLine = [specific_gene, len(training_labels),
        #             mean, std, sk, kurt,
        #             wandbstr] + train_fold_results + valid_fold_results + [train_avg, valid_avg]

        # ===================== FOR RANDOM SPLIT =====================
        # train_loss, valid_loss, wandbstr = train(shuffled_training_inps, shuffled_training_labels, specific_gene, -1, gene_max_len) # ===================== FOR WANDB =====================
        train(
            shuffled_training_inps,
            shuffled_training_labels,
            specific_gene,
            -1,
            gene_max_len,
            tag,
        )  # ===================== FOR CSV =====================
        # print("writing line")
        # writeLine = [specific_gene, len(training_labels),
        #             mean, std, sk, kurt,
        #             wandbstr] + [train_loss, valid_loss]

        # with open(output_path, "a", newline='') as write_to:
        #     writer = csv.writer(write_to)
        #     writer.writerow(writeLine)


def train_single_gene(gene):
    os.environ["TOKENIZERS_PARALLELISM"] = "false"

    # get mapping of gene:num_data, num_species, max_len
    gene_stats = map_gene_to_metadata()

    for gene_path in os.listdir(gene_one2one_datasets_path):
        # extract specific gene from file path
        specific_gene = gene_path.split("_")[0]

        if specific_gene != gene:
            print("skipping", specific_gene)
            continue  # iterate until we find gene of interest

        gene_num_species, gene_max_len = gene_stats[specific_gene][1:]
        gene_num_species, gene_max_len = int(gene_num_species), int(gene_max_len)
        print(
            "starting: ",
            specific_gene,
            " num species: ",
            gene_num_species,
            " max seq len: ",
            gene_max_len,
        )

        # ====================== FOR CONVENIENCE, PADDING ALL GENES TO SAME LENGTH ======================
        gene_max_len = UNIVERSAL_MAX_LEN

        torch.cuda.empty_cache()

        # get tokenized data and labels for gene
        full_gene_path = gene_one2one_datasets_path + gene_path
        training_inps, training_labels = prepare4Enformer(full_gene_path, gene_max_len)
        print("length of data: ", len(training_inps), len(training_labels))

        # shuffle dataset
        indices = np.arange(len(training_labels))
        np.random.shuffle(indices)
        shuffled_training_inps = [training_inps[i] for i in indices]
        shuffled_training_labels = [training_labels[i] for i in indices]

        train(
            shuffled_training_inps,
            shuffled_training_labels,
            specific_gene,
            -1,
            gene_max_len,
        )

"""
Method to train on an amalgamation of genes, then test on single left out gene
"""
def train_amalgamation_of_genes(train_genes, test_genes):
# what we could do is go through genes with decent val loss in prior tests (avg val loss < 11, min val loss < 7.5)
# sample random set of 5 of those genes for training
# sample one random for testing
    results_path = training_dir + "per_gene_results_pt6.csv"

    to_sample = []
    with open(results_path) as read_from:
        for line in read_from:
            line = line.split(",")
            if line[0] == "gene": continue #skip first line
            avg_val, min_val = line[2:]
            gene = line[0].strip()
            avg_val, min_val = int(avg_val), int(min_val)
            if avg_val < 11 and min_val < 7.5: to_sample.append(gene)

    train_genes = random.sample(to_sample, 10)
    test_gene = random.split(to_sample, 1)
    training_inps, training_labels = [], []
    # iterate thru and combine datasets into single list
    for gene in train_genes:
        gene_path = gene_one2one_datasets_path + gene + "_orthologs_trimmed_one2one.csv"
        gene_inps, gene_labels = prepare4Enformer(gene_path, UNIVERSAL_MAX_LEN)
        training_inps.extend(gene_inps)
        training_labels.extend(gene_labels)
    # shuffle training dataset
    indices = np.arange(len(training_labels))
    np.random.shuffle(indices)
    shuffled_training_inps = [training_inps[i] for i in indices]
    shuffled_training_labels = [training_labels[i] for i in indices]
    # create validation set
    val_gene_path = gene_one2one_datasets_path + test_gene[0] +  "_orthologs_trimmed_one2one.csv"
    val_inps, val_labels = prepare4Enformer(val_gene_path, UNIVERSAL_MAX_LEN)

    # call train on the combined data
    train(
        shuffled_training_inps,
        shuffled_training_labels,
        "gene_amalgamation",
        -1,
        UNIVERSAL_MAX_LEN,
        "_".join(train_genes) + f'_{test_gene[0]}',
        (val_inps, val_labels)
    )


if __name__ == "__main__":
    args = sys.argv
    print(args)
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])
