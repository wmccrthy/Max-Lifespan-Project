import os, csv, sys, random
#import numpy as np
#import matplotlib.pyplot as plt # type: ignore

"""
SCRIPT FOR RETRIEVING AND OUTPUTTING TRAINING DATA CSV of (sequence, lifespan) format 

PSUEDO: 

    function takes in num_data parameter 
    (indicating how much training data we want to compile)
    we want to iterate thru cumulative data set according to this parameter
    want to select/use a data point every len(set)/num_data indices 
    interval = len(set)/num_data

    create species:lifespan mappings 
    
    iterate thru cumulative_data in intervals of 10,000 
    (starting from index i where i +=1 for each full iteration 
    of data where we don't scribe (num_data) points)

"""

def print_TOGA_dataset():
    """
    print the first few lines of any dataset for context
    """
    file_to_print = "/data/rsg/chemistry/wmccrthy/Everything/cumulativeTOGAset_trimmed.csv"
    #file_to_print = "/data/rsg/chemistry/maggiejl/P2/prelim_training/data/tokenized_enformer_829124_training.csv"
    #file_to_print = "/data/rsg/chemistry/maggiejl/P2/Max-Lifespan-Project/prelim_training/data/arbitrary_200000_training.csv"
    #file_to_print = "/data/rsg/chemistry/maggiejl/P2/prelim_training/data/metadata_h5py_enformer_829124_training.h5" 
    with open(file_to_print) as read_from:
        done = 0
        for line in read_from:
            if done == 4: return
            print(line[:1000])
            done +=1

def freq_counts(dict):
    frequency_counts = {}
    for freq in dict.values():
        if freq in frequency_counts:
            frequency_counts[freq] += 1
        else:
            frequency_counts[freq] = 1
    return {key: frequency_counts[key] for key in sorted(frequency_counts.keys())}
    
def get_overall_stats(path=None):
    """
    prints the statistics of a dataset
    """
    if path is None: path = '/data/rsg/chemistry/maggiejl/P2/prelim_training/data/top6genes_1000_training.csv'
    lifespans = []
    ids = set()
    seqs = set()
    genes = {}
    species = {}
    with open(path) as read_from:
        for line in read_from:
            line = line.strip().split(",")
            ids.add(line[0])
            lifespans.append(float(line[1]))
            seqs.add(line[2])
            if line[3] in genes: genes[line[3]] += 1
            else: genes[line[3]] = 1
            if line[4] in species: species[line[4]] += 1
            else: species[line[4]] = 1
    lifespans = np.array(lifespans)
    print("overall stats: ", np.mean(lifespans), np.std(lifespans))
    print("unique ids: ", len(ids))
    print("unique seqs: ", len(seqs))
    print("gene makeup: ", genes)
    print("species makeup: ", freq_counts(species))
    
def get_data(num_data, lifespan_min=None, lifespan_max=None, unique_lifespans=0):
    num_data = int(num_data)
    lifespan_path = '/data/rsg/chemistry/wmccrthy/Everything/lifespan_data.csv'
    lifespan_mappings = {}
    
    with open(f'{lifespan_path}') as file:
        for line in file:
            line = line.strip().split(",")
            species_name = line[0]
            lifespan = line[1]
            if lifespan != "Unknown" and lifespan != "Not Established": lifespan_mappings[species_name] = lifespan

    cumulative_set_MIT_path = "/data/rsg/chemistry/wmccrthy/Everything/cumulativeTOGAset_trimmed.csv"

    output_path = f"/data/rsg/chemistry/maggiejl/P2/prelim_training/data/arbitrary_{num_data}_training.csv"
    if lifespan_min and lifespan_max: output_path = f"/data/rsg/chemistry/maggiejl/P2/prelim_training/data/lifespan_range_{lifespan_min}to{lifespan_max}_{num_data}_training.csv"
    #prob: what if not max/min but unique
    if unique_lifespans != 0: output_path = f"/data/rsg/chemistry/maggiejl/P2prelim_training/data/lifespan_range_{lifespan_min}to{lifespan_max}_unique_{num_data}_training.csv"
    valid_data = []
    lifespans_included = set()

    index = 0
    print(len(lifespan_mappings.keys()))
    with open(cumulative_set_MIT_path) as read_from:
        for line in read_from:
            line = line.split(",")
            if len(line) < 9: continue  #avoid blank lines (just newlines for the most part)
            index += 1
            species, sequence = "_".join(line[0].split("_")[:2]), line[-1]
            if species in lifespan_mappings: lifespan = lifespan_mappings[species]
            else: continue
            sequence = sequence.replace("X", "").replace("\n", "").replace('"', "").replace(" ", "").upper()
            lifespan = lifespan.replace("\n", "")
            
            if "-" in sequence or len(sequence) < 1024 or len(sequence) > 12000: continue #limit length of sequence for sake of testing model  
            if unique_lifespans != 0 and lifespan in lifespans_included: continue 
                
            if lifespan_min == None or lifespan_max == None: 
                valid_data.append([index, lifespan, sequence]) 
                lifespans_included.add(lifespan)    
            elif float(lifespan) <= int(lifespan_max) and float(lifespan) >= int(lifespan_min): 
                valid_data.append([sequence, lifespan])  
                lifespans_included.add(lifespan)
    
    print("total valid data to sample from: ", len(valid_data))
    training_random = random.sample(valid_data, num_data)
    with open(output_path, "w") as write_to:
        writer = csv.writer(write_to)
        for d in training_random:
            writer.writerow(d)

def statistics_per_gene():
    """
    calculates the statistics of every gene
    """
    regulatory_set_path = f"/data/rsg/chemistry/wmccrthy/Everything/gene_datasets/regulatory/"
    valid_data = {}
    conf = 0
    genes_data = [["gene", "valid count", "mean", "std", "non121", "missing nucleotides", "unknown", "non_unique", "short"]]
    for file_name in os.listdir(regulatory_set_path):
        valid_species_for_gene = set()
        valid_lifespans_for_gene = []
        non_unique = unknown = missing = non121 = short = 0
        file_path = os.path.join(regulatory_set_path, file_name)
        if "_orthologs_trimmed" not in file_name: continue
        with open(file_path, "r") as read_from:
            for line in read_from:
                #do checks: missing fields, not one to one, unknown lifespan, unknown DNA gaps, short seq
                line = line.split(",")
                if len(line) < 10 or line[0] == 'organism': continue
                if line[3] != 'one2one': 
                    non121 += 1
                    continue
                lifespan, seq, species, gene_type = line[1], line[-1], line[0], line[2]
                gene_type = gene_type.split(".")[1]
                seq = seq.replace(" ", "").replace("\n", "").replace("X", "").replace('"', "").upper()
                lifespan = lifespan.replace("\n", "")
                if lifespan == 'Unknown': 
                    unknown +=1
                    continue 
                lifespan = float(lifespan)
                if "-" in seq: 
                    missing +=1
                    continue #exclude gaps for now, for simplicity sake
                if len(seq) < 10: 
                    short += 1
                    continue #exclude short/empty sequences
                if seq in valid_data: 
                    if lifespan - valid_data[seq][1] < 0.5: 
                        non_unique += 1
                        continue
                        #print("DIFFERENT SPECIES", lifespan, unique_seq[seq][1].replace("\n", ""))
                        #print(unedited_line, unique_seq[seq])
                    conf += 2
                    del valid_data[seq]
                    continue #exclude duplicates of (seq, lifespan)
                
                valid_lifespans_for_gene.append(lifespan)
                valid_species_for_gene.add(species)
                valid_data[seq] = (gene_type, lifespan, species)
        #plots under gene_hist 
        #plt.clf()
        #plt.hist(valid_lifespans_for_gene)
        #x = plt.xlim()[0] + plt.xlim()[1]
        #y = plt.ylim()[0] + plt.ylim()[1]
        #lt.text(3/4*x, 2/5*y, f"total datapoints: {len(valid_lifespans_for_gene)}")
        #lt.text(3/4*x, 3/5*y, f"unique species: {len(valid_species_for_gene)}")
        #plt.savefig(f'/data/rsg/chemistry/maggiejl/P2/prelim_training/data/gene_hist/{g}.png')
        if len(valid_species_for_gene) > 0:
            mean = np.mean(valid_lifespans_for_gene)
            std = np.std(valid_lifespans_for_gene)
            genes_data.append([gene_type, len(valid_lifespans_for_gene), mean, std, non121, missing, unknown, non_unique, short])
    
    valid_lifespans = [tup[1] for tup in valid_data.values()]
    plt.clf()
    plt.hist(valid_lifespans)
    print("mean and std of overall dataset: ", np.mean(valid_lifespans), np.std(valid_lifespans))
    print("duplicates: ", conf)
    print("valid data: ", len(valid_lifespans))
    plt.savefig('/data/rsg/chemistry/maggiejl/P2/prelim_training/data/stats_one2one_{len(valid_lifespans)}_lifespan_hist.png')
    
    output_path = f"/data/rsg/chemistry/maggiejl/P2/prelim_training/data/stats_one2one_{len(valid_lifespans)}_training.csv"
    with open(output_path, "w") as write_to:
        writer = csv.writer(write_to)
        writer.writerows(genes_data)
         
    output_path = f"/data/rsg/chemistry/maggiejl/P2/prelim_training/data/arbitrary_one2one_{len(valid_lifespans)}_training.csv"
    with open(output_path, "w") as write_to:
        writer = csv.writer(write_to)
        for seq in valid_data.keys():
            gene_type, lifespan, species = valid_data[seq]
            writer.writerow([0, lifespan, seq, gene_type, species])
            
def get_specific_data(num_data, gene_types=None):
    """
    DOES SAME AS ABOVE BUT GETS DATA FROM SPECIFIC GENE 
    """
    if gene_types is None: gene_types = ["ZNF860"]
    num_data = int(num_data)
    to_sample = []
    unique_seq = {}
    gene_counts = {key: 0 for key in gene_types}
    index = non_unique = unknown = missing = non121 = short = total = 0
    for gene_type in gene_types:
        regulatory_set_path = f"/data/rsg/chemistry/wmccrthy/Everything/gene_datasets/regulatory/{gene_type}_orthologs_trimmed.csv"
        with open(regulatory_set_path) as read_from:
            for line in read_from:
                #do checks: missing fields, not one to one, unknown lifespan, unknown DNA gaps, short seq
                line = line.split(",")
                if len(line) < 10 or line[0] == 'organism': continue
                total += 1
                if line[3] != 'one2one': 
                    non121 += 1
                    continue
                lifespan, seq, species = line[1], line[-1], line[0]
                seq = seq.replace(" ", "").replace("\n", "").replace("X", "").replace('"', "").upper()
                lifespan = lifespan.replace("\n", "")
                if lifespan == 'Unknown': 
                    unknown +=1
                    continue 
                if "-" in seq: 
                    missing +=1
                    continue #exclude gaps for now, for simplicity sake
                if len(seq) < 10: 
                    short += 1
                    continue #exclude short/empty sequences
                
                if seq in unique_seq: 
                    print(line)
                    print(unique_seq[seq])
                    non_unique += 1
                    continue
                
                index += 1 #adding index to keep format the same
                gene_counts[gene_type] += 1
                unique_seq[seq] = line
                to_sample.append([index, lifespan, seq, gene_type, species])
    
    #print overall statistics    
    print("invalid data stats ... non one2one:", non121, ", missing nucleotides:", missing, "unknown lifespan:", unknown, "repeat seq:", non_unique, "short:", short)
    print("total:", total)
    print("valid data to choose from: ", len(to_sample))
    print("valid data gene breakdown", gene_counts)
    
    #write to file
    num_data = min(num_data, len(to_sample))
    training = random.sample(to_sample, num_data)
    print("creating a length", num_data, "dataset")
    total = len(gene_types)
    output_path = f"/data/rsg/chemistry/maggiejl/P2/prelim_training/data/top{total}genes_{num_data}_training.csv"
    with open(output_path, "w") as write_to:
        writer = csv.writer(write_to)
        for line in training:
            writer.writerow(line)
    #print final statistics
    get_overall_stats(output_path)


def get_sig_data(num_data, top_k=None):
    """
    WORK IN PROGRESS METHOD: The idea here is to compile a training set by 
    sampling from gene sets with average cluster significance below a threshold, 
    for now .25 unweighted avg, .3 weighted avg

    step 1: build map of gene:[cluster stat sig unweighted avg, weighted avg]

    step 2: iterate thru all ortholog files in regulatory sets path 
        - if map[gene][0] <= .25 and map[gene][1] <= .35: add sequences from gene to training output 

    step 3: return random sample from all compiled seqs 
    """
    if top_k: top_k = int(top_k)
    num_data = int(num_data)

    reg_sets_path = "/data/rsg/chemistry/maggiejl/P2/Max-Lifespan-Project/EDA/gene_datasets/orthologs/"
    gene_ranks_path = "/data/rsg/chemistry/maggiejl/P2/Max-Lifespan-Project/prelim_training/gene_rankings.csv"
    gene_sigs = {}
    with open(gene_ranks_path) as read_from:
        for line in read_from:
            line = line.split(",")
            if len(line) < 6 or line[0] == 'gene': continue 
            gene, weighted, unweighted = line[:3]
            gene_sigs[gene] = [float(weighted), float(unweighted)]

    if top_k:
        rankings = sorted(gene_sigs.keys(), key = lambda x:(gene_sigs[x][0], gene_sigs[x][1]))
        print([(rankings[i], i, gene_sigs[rankings[i]]) for i in range(top_k)])
        to_consider = set(rankings[:top_k])
    
    all_data = []
    for file in os.listdir(reg_sets_path):
        gene = file.split("_")[0]
        if gene not in gene_sigs: continue
        if "embedding" not in file:
            if (top_k and gene in to_consider) or (top_k == None and gene_sigs[gene][1] <= .25 and gene_sigs[gene][0] < .35):
                print(file, gene, gene_sigs[gene])
                file_path = reg_sets_path + file
                with open(file_path) as read_from:
                    for line in read_from:
                        line = line.split(',')
                        if len(line) < 10 or len(line[1]) > 5: continue #avoids weird formatting issues and Unknown lifespans 
                        sequence = line[-1]
                        sequence = sequence.replace("X", "").replace("\n", "").replace('"', "").replace(" ", "").upper()
                        if len(sequence) == 0: continue #some outlying sequences 
                        line[-1] = sequence 
                        all_data.append(line)
    
    sampling = random.sample(all_data, min(num_data, len(all_data)))
    print(len(sampling))
    if top_k: new_training_path = f"arbitrary_{len(sampling)}_from_top_{top_k}_genes.csv"
    else: new_training_path = f"arbitrary_{len(sampling)}_from_good_genes.csv"
    with open(new_training_path, "w") as write_to:
        writer = csv.writer(write_to)
        for line in sampling:
            writer.writerow(line)
    
    return 

if __name__ == "__main__":
    args = sys.argv
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])




