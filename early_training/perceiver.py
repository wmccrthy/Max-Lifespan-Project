from perceiver_pytorch import PerceiverIO, Perceiver
import torch, numpy as np, pandas as pd
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


DNABERT_path_local = "/Users/wyattmccarthy/Desktop/MIT Aging Project/Everything/early_training/DNABERT-S"
DNABERT_path = "/data/rsg/chemistry/wmccrthy/Everything/DNABERT-S/"

tokenizer = AutoTokenizer.from_pretrained(DNABERT_path, trust_remote_code=True)

"""
GIVEN SEQUENCE FILE, RETURNS MAX TOKENIZED LENGTH OF SEQ FROM THAT FILE S.T WE KNOW WHERE TO PAD TO
"""
def get_max_len(file_path):
    max_len = 0
    with open(file_path) as to_read:
        for line in to_read:
            line = line.split(",")
            if len(line) < 2: continue 
            lifespan, seq = line[1], line[-1].strip().replace('"',"").replace("\n", "") #make sure indices of seq and lifespan are modified according to what training data is passed 
            # print(seq[:8] + "....", lifespan)
            if len(lifespan) > 5 or len(seq) == 0: continue #indicates lifespan is 'Unknown' 
            seq = tokenizer(seq, return_tensors='pt')['input_ids'] #pad tokenized sequences up to max_len cutoff 
            max_len = max(max_len, len(seq[0]))
    print("Pad up to:", max_len)
    return max_len
            

def prepare4DNABERT_S(file_path, max_len):
    data, labels = [], []
    with open(file_path) as to_read:
        for line in to_read:
            line = line.split(",")
            if len(line) < 2: continue 
            lifespan, seq = line[1], line[-1].strip().replace('"',"").replace("\n", "") #make sure indices of seq and lifespan are modified according to what training data is passed 
            seq = seq.replace('"', "").replace(" ", "")
            # print(seq[:8] + "....", lifespan)
            if len(lifespan) > 5 or len(seq) == 0: continue #indicates lifespan is 'Unknown' 
            # print(len(seq), lifespan)
            seq = tokenizer(seq, return_tensors='pt', padding='max_length', max_length=max_len)['input_ids'] #pad tokenized sequences up to max_len cutoff 
  
            # print(seq, "\n")
            data.append(seq)
            labels.append(float(lifespan))
    return data, torch.Tensor(labels)


class Perceiver_with_DNABERT_S(L.LightningModule):
    def __init__(self, max_seq_len, pad_token, batch_size) -> None:
        super().__init__()
        self.max_len = max_seq_len
        self.pad_token = pad_token
        self.model = PerceiverIO(dim=768,
                                depth=6,
                                num_latents=256,
                                queries_dim=768,
                                latent_dim=256,
                                cross_heads=1,
                                latent_heads=8,
                                cross_dim_head=64,
                                latent_dim_head=64,
                                weight_tie_layers=False,
                                seq_dropout_prob=0)

        self.embedding = AutoModel.from_pretrained(DNABERT_path, trust_remote_code=True)
        for p in self.embedding.parameters(): p.requires_grad = False #freeze embedding layer 

        self.reduce = nn.Linear(512, 1) #shud b 2*latent_dim -> 1

        self.reduce_alt = nn.Linear(256, batch_size) #attempt to see if max pooling is harming algorithm

        self.mask_ignore_token_ids = set(pad_token)

    def mask_padding(self, tokenized_sequence):
        #return boolean tensor w True (or 1) down columns, c, for which tokenized_sequence[c] = PAD 
        #tensor dim needs to align with that of QK^T that is: dimension of (max_len, 768) * (768, max_len) = max_len * max_len 
        one_d = np.array([1 if i == self.pad_token else 0 for i in tokenized_sequence])
        # print(one_d)
        two_d = np.empty((self.max_len, self.max_len))
        for i in range(len(tokenized_sequence)): two_d[:,i] = one_d[i]
        two_d = torch.BoolTensor(two_d)
        # print(two_d.shape)
        return two_d

    def _mask_with_tokens(self, t, token_ids):
        init_no_mask = torch.full_like(t, False, dtype=torch.bool)
        mask = reduce(lambda acc, el: acc | (t == el), token_ids, init_no_mask)
        return mask


    def forward(self, inputs, pad_mask):
        # breakpoint()
        # print("pad mask shape:", pad_mask.shape) #dim = (batch size, 1, max_len)
        # print("inputs shape:", inputs.shape) #dim = (batch size, 1, max_len)

        embeddings = self.embedding(inputs[0])[0] #embed first seq in batch

        for i in range(1, len(inputs)): embeddings = torch.concat((embeddings, self.embedding(inputs[i])[0])) #concat additional seq embeddings to first in batch

        # print("embedding shape:", embeddings.shape) # dim = (batch size, max_len, 768)

        output = self.model(embeddings, mask=pad_mask) #spitting out dim = (batch_size, 256, 256)
        # print("perceiver io output shape:", output.shape)

        max_pool, mean_pool = torch.max(output, dim=0).values, torch.mean(output, dim=0)
        # print("max+mean pool shape:", max_pool.shape, mean_pool.shape) #each of dim [256, 256]

        max_pool, mean_pool = self.reduce_alt(max_pool), self.reduce_alt(mean_pool) #reduce each to [256, batch size]
        # print("max+mean pool reduced shape:", max_pool.shape, mean_pool.shape)

        to_collapse = torch.concat((max_pool, mean_pool)).transpose(0, 1) #concat mean and max pool, then tranpose for dim [batch size, 512]
        # print("concatenated+transposed mean+max pool shape:", to_collapse.shape)
        
        output = self.reduce(to_collapse) #reduce to dim = [batch size, 1]
        return output 

    def training_step(self, batch, batch_idx): 
        # breakpoint()
        # print(batch)
        loss = self._shared_eval_step(batch, batch_idx)
        self.log("train_loss", loss, sync_dist=True) 
        return loss 

    def validation_step(self, batch, batch_idx):
        loss = self._shared_eval_step(batch, batch_idx)
        self.log("val_loss", loss, sync_dist=True)
    

    def _shared_eval_step(self, batch, batch_idx):
        # breakpoint()
        inputs, target = batch
        # print("input:",inputs,"target:",target)

        no_mask = self._mask_with_tokens(inputs, self.mask_ignore_token_ids)
        pad_mask = ~no_mask
        # print(pad_mask)
        
        # target = target.view(-1, 1, 1)

        output = self.forward(inputs, pad_mask)  
        # print("input:",inputs,"target:",target,"prediction:",output)
        # print("target:",target,"prediction:", output)
        # print(target.shape, output.shape)

        criterion = F.mse_loss(output, target)
        # print("criterion:", criterion)

        loss = torch.sqrt(criterion) 
        # print("loss:", loss) 

        return loss

    def configure_optimizers(self):
        return torch.optim.Adam(self.parameters())


def train(training_data, batch_size, num_epochs, num_dev, is_local=False):
    batch_size, num_epochs, num_dev = int(batch_size), int(num_epochs), int(num_dev)

    #STEP 1: TOKENIZE TRAINING DATA 
    DNABERT_S_pad_token = [0, 1, 2, 3]
    
    seq_cutoff = get_max_len(training_data)

    training_inps, training_labels = prepare4DNABERT_S(training_data, int(seq_cutoff))
    training_inps = torch.stack(training_inps)
    
    #avoid hugging face error 
    os.environ["TOKENIZERS_PARALLELISM"] = "false"      

    # breakpoint()
    training_data = TensorDataset(training_inps, training_labels)
     #split data into training and validation set
    train_data_size = int(len(training_data) * .9)
    valid_data_size = len(training_data) - train_data_size
    seed = torch.Generator().manual_seed(42)

    if not is_local: training_data, valid_data = random_split(training_data, [train_data_size, valid_data_size], generator=seed)
    training_data = DataLoader(training_data, batch_size=batch_size, num_workers=19)
    if not is_local: valid_data = DataLoader(valid_data, batch_size=batch_size, num_workers=19)

    
    
    model = Perceiver_with_DNABERT_S(int(seq_cutoff), DNABERT_S_pad_token, batch_size)

    # breakpoint()
    #testing to see how we can freeze embeddings (want to see if this helps model overfit)
    # for p in model.parameters(): print(p.requires_grad)

    # add wandb logger 
    if not is_local: wandb_logger = WandbLogger(log_model="all", project="training-1")
    # cback = ModelCheckpoint(monitor="val_loss", mode="min", save_top_k=1)

    if is_local: trainer = L.Trainer(max_epochs=num_epochs, enable_checkpointing=False, log_every_n_steps=1)
    else:   trainer = L.Trainer(max_epochs = num_epochs,
                                #strategy="ddp_find_unused_parameters_true",
                                accelerator='gpu',
                                devices=num_dev,
                                logger=wandb_logger) 

    if is_local: trainer.fit(model, training_data)
    else: trainer.fit(model, training_data, valid_data)

    return 


if __name__ == "__main__":
    args = sys.argv
    print(args)
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])