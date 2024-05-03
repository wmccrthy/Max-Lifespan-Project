import torch 
from torch import nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset, TensorDataset
from torch import Tensor
from collections import Counter
import math, sys, os, scipy, csv, numpy
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from transformers import AutoTokenizer, AutoModel
import matplotlib.pyplot as plt
from scipy import stats


regulatory_sets_path = "/data/rsg/chemistry/wmccrthy/Everything/gene_datasets/regulatory/"
gene_datasets_path = "/data/rsg/chemistry/wmccrthy/Everything/gene_datasets/"
DNABERT_S_path_local = "/Users/wyattmccarthy/Desktop/MIT Aging Project/Everything/early_training/DNABERT-S"
DNABERT_S_path = "/data/rsg/chemistry/wmccrthy/Everything/DNABERT-S/"

embedder = AutoModel.from_pretrained(DNABERT_S_path, trust_remote_code=True)
tokenizer = AutoTokenizer.from_pretrained(DNABERT_S_path, trust_remote_code=True)


"""
TOKENIZES A TF SET OF SEQUENCES PRIOR TO DNABERT_S embeddings 
"""
def prepare4DNABERT_S(file_path):
    data, labels, species = [], [], []
    with open(file_path) as to_read:
        for line in to_read:
            line = line.split(",")
            if len(line) < 2: continue 
            spcies = " ".join(line[0].split("_")[:2])
            seq, lifespan = line[-1], line[1]
            seq = seq.replace('"', "").replace(" ", "")
            # print(seq[:8] + "....", lifespan)
            if len(lifespan) > 5: continue #indicates lifespan is 'Unknown' 

            seq = tokenizer(seq, return_tensors='pt')['input_ids']
  
            # print(seq, "\n")
            data.append(seq)
            labels.append(float(lifespan))
            species.append(spcies)
    return data, labels, species


"""
CREATES EMBEDDING FILE FOR GIVEN SEQUENCE SET 
    - embeddings are saved as rows of 769 numbers where the first 768 are the embedding and last number is the lifespan label
"""
def embed_tokenized(tokenized_seqs, lifespans, species, file, parent_path, isolated=False):
    embeddings = [] #maybe don't wanna store in memory; can just write to output file immediately (embeddings will line up index-wise with sequence set, that is: first line in embeddings file corresponds to sequence in first line of gene set)
    file_embeddings_path = file.split("_")[0] + "_embeddings.csv"
    if isolated: file_embeddings_path = "/".join(file.split("/")[:-1]) + "/" + file.split("/")[-1].split("_")[0] + "_embeddings.csv"
    file_embeddings_path = parent_path + file_embeddings_path
    
    #if embeddings path exists, don't want to re-write file so return 
    if os.system(f'test -e {file_embeddings_path}') == 0: 
        print('embedding exists for', file, '| exiting')
        return 
    print(file_embeddings_path)
    with open(file_embeddings_path, "w") as write_to:
        writer = csv.writer(write_to) 
        for i in range(len(tokenized_seqs)):
            s, l, spec = tokenized_seqs[i], lifespans[i], species[i]
            # print(s)
            s_embedding = embedder(s)[0]

            # get embeddings with mean/max pooling
            # s_embedding_mean, s_embedding_max = torch.mean(s_embedding[0], dim = 0), torch.max(s_embedding[0], dim = 0).values  #collapse w mean + max pooling where each embedding is simply of dim=768
            # s_embedding = torch.cat([s_embedding_mean, s_embedding_max]) #concatenate mean and max pooling for dim of 1536 
            # collapser = torch.nn.Linear(s_embedding.shape[0], s_embedding.shape[0] // 2) #linear layer to collapse concatenated mean/max embeddings to original dim of 768 
            # s_embedding = collapser(s_embedding)

            #get embeddings with just mean pooling (using this for now as mean/max approach seems to make all embeddings closer together in embedding space)
            s_embedding = torch.mean(s_embedding[0], dim = 0)


            #alt approach is to concatenate mean and max pool then collapse 
            s = [i.item() for i in s_embedding] + [l] + [spec]
            writer.writerow(s)
    write_to.close()
    return 
# test = prepare4DNABERT_S("/Users/wyattmccarthy/Desktop/MIT Aging Project/Everything/early_training/training_data/arbitrary_100_training.csv")[0]
# embed_tokenized(test)


"""
GETS EMBEDDINGS FOR EACH TRANSCRIPTION_FACTOR/GENE SET 
    - iterates through gene sets
    - tokenizes all sequences in a set (calling prepare4DNABERT_S method)
    - embeds tokenized sequences and writes out to csv file (calling embed_tokenized method)

Essentially, method will create DNABERT_S embeddings for all ortholog files in the passed directory 
(expects ortholog file csv format where DNA sequence is in last column (-1 index), species is in first column (0 index) and lifespan is in second column (1 index))


TO RUN: python3 get_regulatory_embeddings.py get_embeddings {dir_of_seqs_to_embed}
"""
def get_embeddings(dir=None):
    if not dir: to_search = regulatory_sets_path
    else: to_search = dir

    for file in os.listdir(to_search):
        file_path = to_search + file

        file_embeddings_path = regulatory_sets_path + file.split("_")[0] + "_embeddings.csv"
        #if embeddings path exists, don't want to re-write file so return or re-incur tokenizing
        if os.system(f'test -e {file_embeddings_path}') == 0: 
            print('embedding exists for', file, '| exiting')
            continue 

        print("Tokenizing", file)
        tokenized, lifespans, species = prepare4DNABERT_S(file_path)
        print("Embedding", file)
        embed_tokenized(tokenized, lifespans, species, file, to_search)
    return 

"""
PERFORMS ABOVE DESCRIBED PROCESS ON A SINGLE, GIVEN FILE (RATHER THAN ITERATING THRU DIRECTORY WITH ALL DATASETS)

(expects ortholog file csv format where DNA sequence is in last column (-1 index), species is in first column (0 index) and lifespan is in second column (1 index))

TO RUN: python3 get_regulatory_embeddings.py get_file_embeddings {file_path} 
"""
def get_file_embeddings(file):
    print("Tokenizing", file)
    tokenized, lifespans, species = prepare4DNABERT_S(file)
    print("Embedding", file)
    embed_tokenized(tokenized, lifespans, species, file, "", True)


"""
FOR EACH EMBEDDING FILE IN GENE_DATASETS or passed directory, PERFORM PCA AND K-MEANS CLUSTERING ON EMBEDDINGS 
THEN, FOR EACH CLUSTER COMPUTE MEAN, STD DEVIATION, Z TEST SCORE, CORRESPONDING STATISTICAL SIGNIFICANCE and MODIFIED Z SCORE 

This method can be applied to any directory that holds embedding files
it expects a csv file format where every value in the line up to index -2 is a value in the embedding | the value in index -2 is the lifespan and value in index -1 is the species 

TO RUN: python3 get_regulatory_embeddings.py cluster_all_embeddings {embedding_directory_path}
"""
def cluster_all_embeddings(dir = None):
    if dir: to_search = dir 
    else: to_search = regulatory_sets_path
    num_analyzed = 0
    error_files = []
    few_sample_files = []
    lifespan_map = get_lifespan_map()

    with open("regulatory_sets_cluster_data.csv", "w") as write_to:
        writer = csv.writer(write_to)
        writer.writerow(['gene type', 'set mean lifespan', 'set std dev', 'set iqr bounds', 'cluster', 'cluster mean lifespan', 'cluster std dev', 'cluster iqr', 'cluster iqr bounds', '#points in cluster', '#species in cluster', 'modified z-score', 'statistical significance'])
        for file in os.listdir(to_search):
            gene_type = file.split("_")[0]
            if 'embeddings_w_species' in file:
                embedding_file = to_search + file
                embeddings, labels, species = [],[], []
                print(f'analyzing {file} (#{num_analyzed}) \n')
                num_analyzed += 1
                with open(embedding_file) as to_read:
                    for line in to_read:
                        line = line.split(",")
                        if len(line) < 2: continue #avoid retarded formatting 
                        embedding, label, s = line[:-2], float(line[-2].replace("\n", "")), line[-1].replace("\n", "")
                        embeddings.append(embedding)
                        labels.append(label)
                        species.append(s)
                to_read.close()
                try:
                    embeddings = numpy.array(embeddings)
                    labels = numpy.array(labels)
                except:
                    print(f'error with embeddings for {file}\n')
                    error_files.append(file)
                    continue 

                #if less than 5 sequences in gene set, skip over, as statistical methods are not robust w so few samples (prob true for a greater threshold than 5 but for sake of avoiding errors we use 5...)
                if len(embeddings) < 5:
                    print(f'too few samples ({len(embeddings)}) in', file, '\n')
                    few_sample_files.append((file, len(embeddings)))
                    continue 

                # reduce to 3 dimensions 
                pca = PCA(n_components=3)
                reduced = numpy.ascontiguousarray(pca.fit_transform(embeddings))

                #approach: cluster w k ranging from 2-5, take max silhouette score across each k and then visualize k-means w highest scoring k 
                max_silhouette = 0
                optimal_k = 0
                for i in range(2,6):
                    try:
                        k_means = KMeans(i)
                        k_means.fit(reduced)
                        silhouette_avg = silhouette_score(reduced, k_means.labels_)
                        if silhouette_avg > max_silhouette: 
                            optimal_k = i
                            max_silhouette = silhouette_avg
                    except:
                        print(f'error with silhouette analysis for {file}\n')
                        optimal_k = 3
                        error_files.append(file)
                        break

                k_means = KMeans(optimal_k)
                k_means.fit(reduced)

                #compute mean, std dev, median, MAD, and IQR for entire set 
                cum_mean, cum_std  = labels.mean(), labels.std()
                cum_median, cum_MAD = numpy.median(labels), stats.median_abs_deviation(labels)
                q1, q3 = numpy.percentile(labels, 25), numpy.percentile(labels, 75)
                iqr = q3 - q1
                upper, lower = q3 + 1.5*iqr, q1 - 1.5*iqr 
    
                species_data_path = f'species_per_gene_cluster_data/{gene_type}'
                os.system(f"mkdir {species_data_path}")
                for label in range(optimal_k):
                    #for each cluster, get avg lifespan, std deviation, z score, modified z score, iqr, and p value  
                    cluster_lifespans = []
                    cluster_species = set()
                    for i in range(len(k_means.labels_)):
                        if k_means.labels_[i] == label: 
                            cluster_lifespans.append(labels[i])
                            cluster_species.add(species[i])
                    cluster_lifespans = numpy.array(cluster_lifespans)
                    cluster_mean, cluster_std, cluster_median = cluster_lifespans.mean(), cluster_lifespans.std(), numpy.median(cluster_lifespans)
                    cluster_q1, cluster_q3 = numpy.percentile(cluster_lifespans, 25), numpy.percentile(cluster_lifespans, 75)
                    cluster_iqr = cluster_q3 - cluster_q1
                    z_score = (cluster_mean - cum_mean)/(cum_std/math.sqrt(len(cluster_lifespans))) #z-score appears to be blowing up bc of non-normal distributions, so waive from analysis for now and instead use modified z score 

                    modified_z_score = .6745*(cluster_median - cum_median)/cum_MAD
                    p_value = stats.norm.sf(abs(modified_z_score))
                    #want to write out gene set, gene set mean, gene set IQR bounds, cluster label, cluster mean lifespan, cluster std dev, cluster iqr, cluster iqr boudns, # points in cluster, #species in cluster, modified z score, stat sig
                    writer.writerow([gene_type, round(cum_mean, 1), round(cum_std, 1), f'{round(q1)}-{round(q3)}', f'cluster{label}',round(cluster_mean, 1), round(cluster_std, 1), round(cluster_iqr, 1), f'{round(cluster_q1)}-{round(cluster_q3)}', len(cluster_lifespans), len(cluster_species), round(modified_z_score, 2), round(p_value, 4)])


                    # write additional output file to gene_datasets/species_per_gene_cluster_data/gene/clusterN_species.csv
                    with open(f"{species_data_path}/cluster{label}_species.csv", "w") as write_to_2:
                        writer_2 = csv.writer(write_to_2)
                        writer_2.writerow(['species', 'lifespan'])
                        for s in cluster_species:
                            s_key = "_".join(s.split(" "))
                            writer_2.writerow([s, lifespan_map[s_key]])
                            
    print("Finished clustering analysis \n")
    print("Errors had on the following:")
    for i in error_files: print(i)
    print("Too few samples in the following:")
    for i in few_sample_files: print(i)
    return 


"""
RETURNS DICTIONARY MAPPING SPECIES:MAX_LIFESPAN | Used in other methods 
"""
def get_lifespan_map():
    lifespan_path = '/data/rsg/chemistry/wmccrthy/Everything/lifespan_data.csv'
    lifespan_mappings = {}
    with open(f'{lifespan_path}') as file:
        for line in file:
            line = line.strip().split(",")
            organism = line[0]
            lifespan = line[1]
            if lifespan == 'Unknown': continue 
            lifespan_mappings[organism] = lifespan
    return lifespan_mappings


"""
ORIGINALLY SAVED EMBEDDINGS WITHOUT THEIR SPECIES LABEL (ONLY INCLUDED LIFESPAN LABEL)
THIS SCRIPT CIRCLES BACK AND ACCURATELY ASSIGNS SPECIES TO EACH EMBEDDING S.T WE CAN INCLUDE SPECIES IN CLUSTERING ANALYSIS 

Embedding scripts have since been updated to include species so this method should not need to be used again
"""
def add_species_embeddings():
    #embeddings shud align index wise w their gene sets 
    #iterate thru all gene set, assign species to each line # 
    #iterate thru corresponding embedding set, and add species assigned to corresponding lines
    files_updated = 0
    for file in os.listdir(regulatory_sets_path):
        if 'embedding' not in file:
            file_gene_type = file.split("_")[0]
            line_num_2_species = {}
            line_num = 0
            with open(f'{regulatory_sets_path}{file}') as read_from:
                for line in read_from:
                    line = line.split(",")
                    if len(line) < 2 or line[0] == 'organism' or len(line[1]) > 5: continue 
                    species = " ".join(line[0].split("_")[:2])
                    line_num_2_species[line_num] = species
                    line_num += 1
            expected_entries = line_num
            file_embeddings_path = f'{regulatory_sets_path}{file_gene_type}_embeddings.csv'
            file_embeddings_path_alt = f'{regulatory_sets_path}{file_gene_type}_embeddings_w_species.csv'
            files_updated += 1

            if os.system(f'test -e {file_embeddings_path_alt}') == 0: 
                print('embedding w species exists for', file, '| exiting')
                continue 

            print("creating", file_embeddings_path_alt, f'(#{files_updated})')
            with open(file_embeddings_path_alt, 'w') as write_to:
                writer = csv.writer(write_to)
                with open(file_embeddings_path) as read_from:
                    line_num = 0
                    for line in read_from:
                        line = line.split(",")
                        if len(line) < 2: continue #avoid faulty lines
                        line[-1] = line[-1].replace("\n", "")
                        # print(line[-1], line_num_2_species[line_num])
                        try:
                            line = line + [line_num_2_species[line_num]]
                        except: 
                            print("missing entry for line", line_num, "| cutting off")
                            break
                        line_num += 1
                        writer.writerow(line)


"""
GIVEN AN EMBEDDING FILE, PERFORM PCA AND K-MEANS CLUSTERING ON EMBEDDINGS THEN VISUALIZE OUTPUT w Matplot

it expects a csv file format where every value in the line up to index -2 is a value in the embedding | the value in index -2 is the lifespan and value in index -1 is the species 

TO RUN: python3 get_regulatory_embeddings.py visualize_embeddings {embedding_file_path}
"""
def visualize_embeddings(embedding_file, is2d = False):
    embeddings, labels = [],[]
    with open(embedding_file) as to_read:
        for line in to_read:
            line = line.split(",")
            if len(line) < 2: continue #avoid retarded formatting 
            embedding, label = line[:-2], float(line[-2].replace("\n", ""))
            embeddings.append(embedding)
            labels.append(label)
    embeddings = numpy.array(embeddings)

    # for testing purposes, reduce to 2 or 3 dimensions (play around w both)
    if not is2d: pca = PCA(n_components=3)
    else: pca = PCA(n_components=2)

    reduced = numpy.ascontiguousarray(pca.fit_transform(embeddings))

    #try k means clustering on reduced 
    #approach: cluster w k ranging from 2-5, take max silhouette score across each k and then visualize k-means w highest scoring k 
    max_silhouette = 0
    optimal_k = 0
    for i in range(2,6):
        k_means = KMeans(i)
        k_means.fit(reduced)
        silhouette_avg = silhouette_score(reduced, k_means.labels_)
        # print(i, silhouette_avg)
        if silhouette_avg > max_silhouette: 
            optimal_k = i
            max_silhouette = silhouette_avg

    # print("optimal k:", optimal_k)

    k_means = KMeans(optimal_k)
    k_means.fit(reduced)

    # print(reduced)
    # print(len(reduced), num_seqs)

    fig = plt.figure()

    
    if is2d: grph = fig.add_subplot() #for 2d
    else: grph = fig.add_subplot(projection='3d') #for 3d 

    if is2d: grph.scatter(reduced[:,0], reduced[:,1], c=k_means.labels_) #for 2d
    else: grph.scatter(reduced[:,0], reduced[:,1], reduced[:,2], c=k_means.labels_) #for 3d

    #for each cluster, display avg lifespan and std deviation
    for label in range(optimal_k):
        #get average lifespan of 
        lifespans = []
        for i in range(len(k_means.labels_)):
            if k_means.labels_[i] == label: lifespans.append(labels[i])
        lifespans = numpy.array(lifespans)

        if not is2d: grph.text(k_means.cluster_centers_[label][0], k_means.cluster_centers_[label][1], k_means.cluster_centers_[label][2], f'{round(lifespans.mean())}|{round(lifespans.std())}', bbox=dict(facecolor='white', alpha=0.4)) #for 3d
        else: grph.text(k_means.cluster_centers_[label][0], k_means.cluster_centers_[label][1], f'{round(lifespans.mean())}|{round(lifespans.std())}', bbox=dict(facecolor='white', alpha=0.4)) #for 2d

    # [grph.text(reduced[i][0], reduced[i][1], reduced[i][2], labels[i]) for i in range(len(labels))] #textually label 
    #try to textually label s.t only the mean is 
    grph.set_title(embedding_file.split("/")[-1][:-4] + f" Visualized (mean lifespan={round(numpy.mean(labels))})")
    plt.show()
    fig.savefig("embeddings_visuals/" + embedding_file.split("/")[-1][:-4] + "_visual.png")
    plt.close()


"""
Simple script to collect basic metrics; COUNTS NUMBER OF EMBEDDING FILES IN regulatory_sets_path
"""
def count_embedding_sets():
    embedded_cnt = 0
    total = 0
    for file in os.listdir(regulatory_sets_path):
        if 'embedding' in file: embedded_cnt += 1
        else: total += 1
    print(f'{embedded_cnt}/{total} regulatory sets embedded')

if __name__ == "__main__":
    args = sys.argv
    print(args)
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])

