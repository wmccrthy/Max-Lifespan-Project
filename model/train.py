import os
import torch, wandb, numpy as np, pandas as pd
from torch import nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset, TensorDataset
from torch import Tensor
import lightning as L
from pytorch_lightning.loggers import WandbLogger
from torchmetrics.functional import accuracy
from collections import Counter
import math


class Sequence_Tokenizer():
    """
    DATASET CLASS BUILT TO CONVERT (SEQUENCE, LIFESPAN) TABULAR DATA TO TOKENIZED FORMAT

    ARGS 
        - file_path: refers to file path containing csv data we want to tokenize
        - truncate_at: refers to the length at which we want to truncate sequences to ensure they are all of same length
        - token_len: desired token length (where token_len = 1 will have each base (A, T, G, C) as a token, token_len = 2 will have each pair of bases (AT, AC, AG, ..., TT))
    """
    def __init__(self, file_path, truncate_at, token_len):
        self.token_len = token_len
        self.truncate_at = truncate_at

        #find min length sequence s.t we know what length to truncate to 
        min_seq = float('inf')
        max_seq = 0

        #need to indicize text (associate words to indices)
        self.tokens = set() #n is special char
        #need to understand significance of incorporating (or not) 'N' 
        #N denotes that any nucleotide may be there 

        self.compile_tokens("", self.tokens)
        self.tokens_idx = {token:int(idx) for idx, token in enumerate(self.tokens)}
        self.idx_tokens = {val:key for key, val in self.tokens_idx.items()}

        self.vocab_size = len(self.tokens_idx)

        self.token_labels = {}

        print(self.tokens_idx)
        
        # vectorized_seqs = []
        seqs = []
        lifespan_labels = []
        with open(file_path) as read_from:
            for line in read_from:
                line = line.split(",")
                if len(line) < 2: 
                    # print(" ".join(line))
                    continue 
                seq, lifespan = line
                lifespan = lifespan.replace("\n", "")
                seq = seq.replace(" ", "").replace('"',"").replace("'", "")
                min_seq = min(len(seq), min_seq)
                max_seq = max(len(seq), max_seq)
                # print(len(seq))
                tokenized_seq = self.tokenize(seq)
                seqs.append(tokenized_seq)
                lifespan_labels.append(float(lifespan))
                self.token_labels[tuple(tokenized_seq)] = lifespan
            
        # print(self.decode_token(seqs[0])) #debugging for decode_token method 

        self.x = torch.LongTensor(seqs)
        self.y = torch.Tensor(lifespan_labels)
        print(self.x, self.y)
        # self.x = seqs
        # self.y = lifespan_labels

        self.tokenized_seq_len = len(self.x[0])

        self.emb_dim = int(truncate_at / token_len) #length of a tokenized sequence is the embedding dimension

        # print(self.x, self.y)       

        #i have raw sequences, raw lifespans, how do i go from that to what is accepted by model? 
        #if not transformer, i cud do BOW for each sequence, but how do we also encode positional information within embeddings? 
     
    def compile_tokens(self, cur_token, vocab):
        #dfs for token_len chars, taking path for each base (A, T, G, C) at each step
        if len(cur_token) == self.token_len: 
            vocab.add(cur_token)
            return
        self.compile_tokens(cur_token + 'A', vocab)
        self.compile_tokens(cur_token + 'G', vocab)
        self.compile_tokens(cur_token + 'T', vocab)
        self.compile_tokens(cur_token + 'C', vocab)
        self.compile_tokens(cur_token + 'N', vocab)

    def tokenize(self, seq):
        # print(Counter(seq)) #debugging
        return [self.tokens_idx[seq[i:i+self.token_len]] for i in range(0, self.truncate_at+1-self.token_len, self.token_len)]

    def decode_token(self, tokenized):
        return "".join([self.idx_tokens[i] for i in tokenized])
    
    def get_label(self, token):
        return self.token_labels[token]
        
"""
POSITIONAL ENCODING CLASS
"""
class PositionalEncoding(nn.Module):
    def __init__(self, d_model: int, dropout: float = 0.1, max_len: int = 5000):
        super().__init__()
        self.dropout = nn.Dropout(p=dropout)

        position = torch.arange(max_len).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2) * (-math.log(10000.0) / d_model))
        pe = torch.zeros(max_len, 1, d_model)
        pe[:, 0, 0::2] = torch.sin(position * div_term)
        pe[:, 0, 1::2] = torch.cos(position * div_term)
        self.register_buffer('pe', pe)

    def forward(self, x: Tensor) -> Tensor:
        """
        Arguments:
            x: Tensor, shape ``[seq_len, batch_size, embedding_dim]``
        """
        x = x + self.pe[:x.size(0)]
        return self.dropout(x)


"""
TRANSFORMER MODEL CLASS
"""
class Transformer(L.LightningModule): 
    def __init__(self, tokenized_seq_len, vocab_size, att_heads, num_layers, learning_rate) -> None:
        super().__init__()
        #cud also define individual layers 
        self.vocab_size = vocab_size

        self.d_model = tokenized_seq_len * 2 
        #come up w nice function for this that is dependent on vocab size and tokenized seq length; given there is an inverse relationship there... (as vocab size increase, tokenized_seq_len decrease)


        self.encoder_layer = nn.TransformerEncoderLayer(d_model=self.d_model, nhead=att_heads)
        self.transformer_encoder = nn.TransformerEncoder(self.encoder_layer, num_layers=num_layers)
        self.pos_encoder = PositionalEncoding(self.d_model, max_len=tokenized_seq_len)

        self.embedding = nn.Embedding(vocab_size, self.d_model)

        # print(self.embedding.weight)

        #for use in forward approach 1 
        self.linear1 = nn.Linear(self.d_model, self.d_model // 2)
        self.linear2 = nn.Linear(self.d_model // 2, 1)
        self.reduce = nn.Linear(tokenized_seq_len, 1)

        #for use in forward approach2 
        self.reduce2 = nn.Linear(tokenized_seq_len * 2, 1)


        self.testing_preds = []
        
        self.lr = learning_rate

    def forward(self, inputs):
        #we want to first embed inputs    

        # breakpoint()

        # print("Input Shape:", inputs.shape, inputs)
        inputs = self.embedding(inputs) #first embed input

        # print("Sequence Embedding Shape:", inputs.shape)
        #then, we need to add positional encodings to the embeddings -> sus about this 
        final_embeddings = self.pos_encoder(inputs) 
       
        # print("Pre-Encoder (Seq Embedding + Pos Encoding) Shape:", final_embeddings.shape)
        #then, given we are using a pre-implemented transformer model, we just pass the embedded input to our encoder 
        output = self.transformer_encoder(final_embeddings)
        # print("Post-Encoder Shape:", output.shape)
        
        #approach 1: put thru linear1 and linear2 layer (outshape of 128 x 1), than transpose output s.t it's shape 1x128, then put thru one more linear layer to reduce to 1x1
        # output = self.linear1(output)
        # print("Post Linear1 Shape:", output.shape)
        # output = self.linear2(output)
        # print("Post Linear2 Shape:", output.shape)
        # output = output.transpose(1, 2)
        # print("Post Tranpose Shape:", output.shape, output)
        # output = self.reduce(output)
        # print("Post Reduction Shape:", output.shape)


        #approach 2: instead of putting thru linear layers initially, concatenate maxpool2d and avgpool2d of encoder output providing output shape of 1x128, then put thru one linear layer to reduce to 1x1
        max_pool = nn.MaxPool2d((1, self.d_model))
        output_mp = max_pool(output)
        # print("Max Pool Shape:", output_mp.shape)


        avg_pool = nn.AvgPool2d((1, self.d_model))
        output_ap = avg_pool(output)
        # print("Avg Pool Shape:", output_ap.shape)

        output = torch.cat((output_mp, output_ap), 1).transpose(1, 2)
        # print("MaxPool Concat Avg Pool Shape:", output.shape)

        output = self.reduce2(output)
        # print("Post Reduction Shape:", output.shape)

        return output

    def training_step(self, batch, batch_idx):
        # breakpoint()

        inputs, target = batch

        target = target.view(-1, 1, 1)
        # print(target.shape, target)
        # print(inputs, target)

        output = self.forward(inputs)

        print("Target:", target, "Train Output:", output)

        criterion = F.mse_loss(output, target)
        loss = torch.sqrt(criterion)

        print("Loss: ", loss.item())
        return loss

    def test_step(self, batch, batch_idx):
        loss = self._shared_eval_step(batch, batch_idx)
        metrics = {'test_loss':loss}
        return metrics 

    def _shared_eval_step(self, batch, batch_idx):
        inputs, target = batch
        target = target.view(-1, 1, 1)
        # print(inputs, target)
        output = self.forward(inputs)
        # breakpoint()
        print("Target:", target, "Test Prediction:", output)

        self.testing_preds.append((target, output))

        criterion = F.mse_loss(output, target)
        loss = torch.sqrt(criterion)

        print("Loss: ", loss.item())
        # acc = accuracy(output, target, 'multilabel')
        # print("Testing Loss: ", loss.item(), "Testing Accuracy: ", acc)
        return loss

    def configure_optimizers(self):
        return torch.optim.Adam(self.parameters(), lr=self.lr)


"""
SCRIPT
"""
remote_training_set_path = "/Mounts/rbg-storage1/users/wmccrthy/early_training/arbitrary_100000_training.csv" 
remote_testing_set_path = "/Mounts/rbg-storage1/users/wmccrthy/early_training/lifespan_range_0to150_100_training.csv"

local_training_path = "/Users/wyattmccarthy/Desktop/MIT Aging Project/Everything/early_training/arbitrary_100_training.csv"
# local_training_path = "/Users/wyattmccarthy/Desktop/MIT Aging Project/Pre-Processing/early_training/lifespan_range_6to12_100_training.csv"
# local_training_path = "/Users/wyattmccarthy/Desktop/MIT Aging Project/Pre-Processing/early_training/lifespan_range_0to20_500_training.csv"

# local_testing_path = "/Users/wyattmccarthy/Desktop/MIT Aging Project/Pre-Processing/early_training/lifespan_range_0to20_unique_50_training.csv"
local_testing_path = "/Users/wyattmccarthy/Desktop/MIT Aging Project/Everything/early_training/lifespan_range_0to150_100_training.csv"

# breakpoint()
"""
BELOW IS SCRIPT FOR TESTING AND COMPARING VARIOUS APPROACHES TO TRAINING THE MODEL 
EACH USES A DIFF BASE TOKEN LENGTH WHERE: A TOKEN LENGTH OF 1 HAS A VOCABULARY OF 5 COMPRISED OF EACH NUCLEOTIDE (A,T,G,C,N) N can represent any nucleotide
                                          A TOKEN LENGTH OF 3 HAS A VOCABULARY OF 125 COMPRISED OF EACH CODON (AAA, ..., NNN) all length-3 permutations of 4 bases (A,T,G,C) and N
                                          ...
"""

results = {}

for i in range(1, 4):
    data = Sequence_Tokenizer(remote_training_set_path, 1024, i)
    test_d = Sequence_Tokenizer(remote_testing_set_path, 1024, i)

    dataset = TensorDataset(data.x, data.y)
    dataloader = DataLoader(dataset)

    testset = TensorDataset(test_d.x, test_d.y)
    testloader = DataLoader(testset)


    model = Transformer(data.tokenized_seq_len, data.vocab_size, 2,  4, .02)

    wandb_logger = WandbLogger(log_model="all", 
                                project="training-1")

    trainer = L.Trainer(max_epochs = 1,
                            accelerator='gpu',
                            devices=1,
                            logger=wandb_logger) 
    #remember to add configs for using GPUs within 


    trainer.fit(model, dataloader)
    print("Training Data Mean:", data.y.mean()) #testing to see if model is just getting stuck at average across all lifespans 

    trainer.test(model, testloader)

    results[i] = model.testing_preds
    print("RESULTS FROM TRAINING WITH TOKEN LENGTH:", i)
    print("Testing Data Mean:", test_d.y.mean())
    print("Training Data Mean:", data.y.mean())

for i in results:
    print("Token Length", i, "Testing Results:", results[i])


