#import pickle
#import json
#import re
import csv
import torch
import numpy as np
import h5py
import time
           
# tokenizing to a file
# enformer 829124 
# max length 22785 but enformer required length 196608

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
    if max_len - len(newseq) > 0:
        newseq = newseq + [-1]*(max_len - len(newseq))
    return newseq

def prepare4Enformer(max_len):
    """
    GIVEN SEQUENCE FILE, RETURNS MAX TOKENIZED LENGTH OF SEQ FROM THAT FILE S.T WE KNOW WHERE TO PAD TO
    """
    data_path = "/data/rsg/chemistry/maggiejl/P2/prelim_training/data/arbitrary_one2one_829124_training.csv"
    output_path = "/data/rsg/chemistry/maggiejl/P2/prelim_training/data/tokenized_enformer_829124_training.csv"
    data_batch = []
    index = 0

    with open(data_path) as to_read:
        for line in to_read:
            line = line.split(",")
            if len(line) < 2: continue 
            lifespan, seq = line[1], line[2].strip().replace('"',"").replace("\n", "") #make sure indices of seq and lifespan are modified according to what training data is passed 
            gene, species = line[3], line[4].strip().replace("/n", "")
            seq = seq.replace('"', "").replace(" ", "")
            if len(lifespan) > 5 or len(seq) == 0: 
                print("indicates lifespan is unknown")
                continue #indicates lifespan is 'Unknown' 
            data_batch.append([gene, species, float(lifespan), transform4Enformer(seq, max_len)])
            
            if len(data_batch) > 1000:
                index += 1
                with open(output_path, "a", newline='') as write_to:
                    writer = csv.writer(write_to)
                    writer.writerows(data_batch)
                data_batch = []
                print("batch ", index, " completed")
    if len(data_batch) > 0:
        with open(output_path, "a", newline='') as write_to:
            writer = csv.writer(write_to)
            writer.writerows(data_batch)

def loadTokenized(num_data=829124):
    start_time = time.time()
    data_path = "/data/rsg/chemistry/maggiejl/P2/prelim_training/data/tokenized_enformer_829124_training.csv"
    data = []
    labels = []
    genes = []
    ind = 0
    non_line = 0
    pattern = re.compile(r'([^,]+),[^,]+,([^,]+),"(.*)"')
    with open(data_path) as to_read:
        for line in to_read:
            ind += 1
            if len(data) >= num_data: 
                print("time for everything: ", time.time() - start_time)
                return data, labels, genes
            if len(data) % 1000 == 0: print("datapoint ", len(data), " done")
            match = pattern.match(line)
            if not match:
                non_line += 1
                continue
            gene, lifespan, seq = match.groups()
            seq = seq.strip()
            seq = json.loads(seq)
            labels.append(float(lifespan))
            data.append(torch.tensor(seq))
            genes.append(gene)
            #print(seq[:30], type(seq[0]))
    print("finished loading in!")
    print("total data: ", ind)
    print("did not match: ", non_line)
    
    print(genes)
    print(labels)
    print(data[0][:20])
    print(type(data[0]), type(data[0][0]))
    print("time for everything: ", time.time() - start_time)
    #for seq in data:
    #    seq = ast.literal_eval(seq)
    #print("haha")
    #print(data[0][:4])
    #print(type(data[0]), type(data[0][0]))
    return data, labels, genes

def prepare4EnformerPKL(max_len):
    """
    GIVEN SEQUENCE FILE, RETURNS MAX TOKENIZED LENGTH OF SEQ FROM THAT FILE S.T WE KNOW WHERE TO PAD TO
    """
    data_path = "/data/rsg/chemistry/maggiejl/P2/prelim_training/data/arbitrary_one2one_829124_training.csv"
    output_path = "/data/rsg/chemistry/maggiejl/P2/prelim_training/data/tokenized_enformer_829124_training.pkl"
    meta_data_path = "/data/rsg/chemistry/maggiejl/P2/prelim_training/data/metadata_enformer_829124_training.csv"
    data = []
    meta_data = []
    index = 0
    with open(meta_data_path, 'w', newline='') as file: pass

    with open(data_path) as to_read:
        for line in to_read:
            line = line.split(",")
            if len(line) < 2: continue 
            lifespan, seq = line[1], line[2].strip().replace('"',"").replace("\n", "") #make sure indices of seq and lifespan are modified according to what training data is passed 
            gene, species = line[3], line[4].strip().replace("/n", "")
            seq = seq.replace('"', "").replace(" ", "")
            if len(lifespan) > 5 or len(seq) == 0: 
                print("indicates lifespan is unknown")
                continue #indicates lifespan is 'Unknown' 
            meta_data.append([gene, species])
            data.append([float(lifespan), transform4Enformer(seq, max_len)])
            
            #if len(data) > 3000: break
            if len(data) % 1000 == 0: 
                index += 1
                print("batch ", index, " completed")
    print("saving to pkl")
    with open(output_path, "wb") as file:
        pickle.dump(data, file)
    print("saving to csv")
    with open(meta_data_path, "a", newline='') as write_to:
        writer = csv.writer(write_to)
        writer.writerows(meta_data)

def loadTokenizedPKL():
    start_time = time.time()
    pickle_path = "/data/rsg/chemistry/maggiejl/P2/prelim_training/data/tokenized_enformer_829124_training.pkl"
    #takes the longest
    with open(pickle_path, "rb") as file:
        data = pickle.load(file)
    print("time for pickle reading: ", time.time() - start_time)
    print(len(data))
    print(type(data[0][0]), type(data[0][1]), type(data[0][1][0]))
    print(data[0][0], data[0][1][:20])
    print('starting labels')
    labels = []
    seqs = []
    for entry in data:
        labels.append(entry[0])
        seqs.append(np.array(entry[1], dtype=np.int64))
    print("tensoring them")
    labels = torch.tensor(labels)
    print('tensoring again')
    for seq in seqs:
        if len(seq) != 196608: print("uh oh")
    print("length of all checked")
    print(type(seqs))
    print(len(seqs))
    print(len(seqs[0]))
    print(type(seqs[0]))
    print(type(seqs[0][0]))
    print("trying to array") #takes the longest
    seqs = np.array(seqs)
    print("numpied arrayed")
    seqs = torch.from_numpy(seqs)
    print("torched!!!")
    print(seqs.dtype)
    print("time for everything: ", time.time() - start_time)

def prepare4Enformer_h5py(max_len):
    """
    GIVEN SEQUENCE FILE, RETURNS MAX TOKENIZED LENGTH OF SEQ FROM THAT FILE S.T WE KNOW WHERE TO PAD TO
    """
    data_path = "/data/rsg/chemistry/maggiejl/P2/prelim_training/data/arbitrary_one2one_829124_training.csv"
    output_path = "/data/rsg/chemistry/maggiejl/P2/prelim_training/data/tokenized_enformer_829124_training.h5"
    meta_data_path = "/data/rsg/chemistry/maggiejl/P2/prelim_training/data/metadata_h5py_enformer_829124_training.csv"
    data_batch = []
    meta_data_batch = []
    index = 0

    #create empty h5py file with max_len + 1 spaces (lifespan, seq)
    with h5py.File(output_path, 'a') as file:
        if 'dataset' in file:
            del file['dataset']
        dataset = file.create_dataset('dataset', data=np.empty((0, max_len)), dtype=np.int64, maxshape=(None, max_len), chunks=True)

    with open(data_path) as to_read:
        for line in to_read:
            line = line.split(",")
            if len(line) < 2: continue 
            lifespan, seq = line[1], line[2].strip().replace('"',"").replace("\n", "") #make sure indices of seq and lifespan are modified according to what training data is passed 
            gene, species = line[3], line[4].strip().replace("/n", "")
            seq = seq.replace('"', "").replace(" ", "")
            if len(lifespan) > 5 or len(seq) == 0: 
                print("indicates lifespan is unknown")
                continue #indicates lifespan is 'Unknown' 
            data_batch.append(transform4Enformer(seq, max_len))
            meta_data_batch.append([float(lifespan), gene, species])
            
            if len(data_batch) > 1000:
                index += 1
                #if index > 2: return #early stop for testing
                with h5py.File(output_path, 'a') as file:
                    existing_dataset = file['dataset']
                    existing_dataset.resize((existing_dataset.shape[0] + len(data_batch)), axis=0)
                    existing_dataset[-len(data_batch):] = np.array(data_batch)
                with open(meta_data_path, "a", newline='') as write_to:
                    writer = csv.writer(write_to)
                    writer.writerows(meta_data_batch)
                data_batch = []
                meta_data_batch = []
                print("batch ", index, " completed")
    if len(data_batch) > 0:
        
        with h5py.File(output_path, 'a') as file:
            existing_dataset = file['dataset']
            existing_dataset.resize((existing_dataset.shape[0] + len(data_batch)), axis=0)
            existing_dataset[-len(data_batch):] = np.array(data_batch)
        with open(meta_data_path, "a", newline='') as write_to:
            writer = csv.writer(write_to)
            writer.writerows(meta_data_batch)

def loadTokenized_h5py(num_data=829124):
    output_path = "/data/rsg/chemistry/maggiejl/tokenized_enformer_829124_training.h5"
    
    start_time = time.time()
    with h5py.File(output_path, 'r') as file:
        data = file['dataset']
        print(data.shape)
        print(type(data[:1, :20]), data[:1, :20])
        print(data[0, 0], type(data[0, :]), type(data[0, 0]))
        test = data[7, :]
        print("lalala", test[:20])
        print("time to load in hf: ", time.time() - start_time)
        datapoints = data[:, :]
        print("time for loading in everything: ", time.time() - start_time)
    print("exiting loop")
    print(type(datapoints))
    print(datapoints.shape)
    
    
    print('starting new run')
    start_time = time.time()
    hf = h5py.File(output_path, 'r')
    datapoints = hf.get('dataset')
    print("time to load hf: ", time.time() - start_time)
    seq = datapoints[1, :]
    print("time to load one datapoint: ", time.time() - start_time)
    seq = datapoints[:, :]
    print("time to load all rows: ", time.time() - start_time)
    torch.tensor(seq, dtype=torch.int64)
    print("time to convert to tensor: ", time.time() - start_time)
    

if __name__ == "__main__":
    #tokenizing = False
    #if tokenizing: prepare4Enformer_h5py(196608)
    #else: 
        loadTokenized_h5py(3000)
    