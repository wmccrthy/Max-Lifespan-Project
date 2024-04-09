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
DNABERT_S_path_local = "/Users/wyattmccarthy/Desktop/MIT Aging Project/Everything/early_training/DNABERT-S"
DNABERT_S_path = "/data/rsg/chemistry/wmccrthy/Everything/DNABERT-S/"

embedder = AutoModel.from_pretrained(DNABERT_S_path_local, trust_remote_code=True)
tokenizer = AutoTokenizer.from_pretrained(DNABERT_S_path_local, trust_remote_code=True)


"""
TOKENIZES A TF SET OF SEQUENCES PRIOR TO DNABERT_S embeddings 
"""
def prepare4DNABERT_S(file_path):
    data, labels = [], []
    with open(file_path) as to_read:
        for line in to_read:
            line = line.split(",")
            if len(line) < 2: continue 
            seq, lifespan = line[-1], line[1]
            seq = seq.replace('"', "").replace(" ", "")
            # print(seq[:8] + "....", lifespan)
            if len(lifespan) > 5: continue #indicates lifespan is 'Unknown' 

            seq = tokenizer(seq, return_tensors='pt')['input_ids']
  
            # print(seq, "\n")
            data.append(seq)
            labels.append(float(lifespan))
    return data, labels


"""
CREATES EMBEDDING FILE FOR GIVEN SEQUENCE SET 
    - embeddings are saved as rows of 769 numbers where the first 768 are the embedding and last number is the lifespan label
"""
def embed_tokenized(tokenized_seqs, lifespans, file, parent_path, isolated=False):
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
            s, l = tokenized_seqs[i], lifespans[i]
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
            s = [i.item() for i in s_embedding] + [l]
            writer.writerow(s)

# test = prepare4DNABERT_S("/Users/wyattmccarthy/Desktop/MIT Aging Project/Everything/early_training/training_data/arbitrary_100_training.csv")[0]
# embed_tokenized(test)


"""
GETS EMBEDDINGS FOR EACH TRANSCRIPTION_FACTOR/GENE SET 
    - iterates through gene sets
    - tokenizes all sequences in a set (calling prepare4DNABERT_S method)
    - embeds tokenized sequences and writes out to csv file (calling embed_tokenized method)
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
        tokenized, lifespans = prepare4DNABERT_S(file_path)
        print("Embedding", file)
        embed_tokenized(tokenized, lifespans, file, to_search)
    return 

"""
PERFORMS ABOVE DESCRIBED PROCESS ON A SINGLE, GIVEN FILE (RATHER THAN ITERATING THRU DIRECTORY WITH ALL DATASETS)
"""
def get_file_embeddings(file):
    print("Tokenizing", file)
    tokenized, lifespans = prepare4DNABERT_S(file)
    print("Embedding", file)
    embed_tokenized(tokenized, lifespans, file, "", True)




"""
FOR EACH EMBEDDING FILE IN GENE_DATASETS, PERFORM PCA AND K-MEANS CLUSTERING ON EMBEDDINGS 
THEN, FOR EACH CLUSTER COMPUTE MEAN, STD DEVIATION, Z TEST SCORE AND CORRESPONDING STATISTICAL SIGNIFICANCE 
"""
def cluster_all_embeddings(dir = None):
    if dir: to_search = dir 
    else: to_search = regulatory_sets_path
    num_analyzed = 0
    error_files = []
    with open("regulatory_sets_cluster_data.csv", "w") as write_to:
        writer = csv.writer(write_to)
        writer.writerow(['gene type', 'set mean lifespan', 'cluster', 'cluster mean lifespan', 'cluster std dev', '# points in cluster', 'z-score', 'statistical significance'])
        for file in os.listdir(to_search):
            gene_type = file.split("_")[0]
            if 'embedding' in file:
                embedding_file = to_search + file
                embeddings, labels = [],[]
                print(f'analyzing {file} (#{num_analyzed}) \n')
                num_analyzed += 1
                with open(embedding_file) as to_read:
                    for line in to_read:
                        line = line.split(",")
                        if len(line) < 2: continue #avoid retarded formatting 
                        embedding, label = line[:-1], float(line[-1].replace("\n", ""))
                        embeddings.append(embedding)
                        labels.append(label)
                to_read.close()
                try:
                    embeddings = numpy.array(embeddings)
                except:
                    print(f'error with embeddings for {file}')
                    error_files.append(file)
                    continue 

                # reduce to 3 dimensions 
                pca = PCA(n_components=3)
                reduced = numpy.ascontiguousarray(pca.fit_transform(embeddings))

                #try k means clustering on reduced 
                #approach: cluster w k ranging from 2-5, take max silhouette score across each k and then visualize k-means w highest scoring k 
                max_silhouette = 0
                optimal_k = 0
                for i in range(2,6):
                    k_means = KMeans(i)
                    k_means.fit(reduced)
                    silhouette_avg = silhouette_score(reduced, k_means.labels_)
                    if silhouette_avg > max_silhouette: 
                        optimal_k = i
                        max_silhouette = silhouette_avg

                k_means = KMeans(optimal_k)
                k_means.fit(reduced)

                #for each cluster, get avg lifespan, std deviation and z score 
                for label in range(optimal_k):
                    cluster_lifespans, all_lifespans = [], []
                    for i in range(len(k_means.labels_)):
                        if k_means.labels_[i] == label: cluster_lifespans.append(labels[i])
                        all_lifespans.append(labels[i])
                    cluster_lifespans, all_lifespans = numpy.array(cluster_lifespans), numpy.array(all_lifespans)
                    cluster_mean, cluster_std = cluster_lifespans.mean(), cluster_lifespans.std()
                    cum_mean, cum_std  = all_lifespans.mean(), all_lifespans.std()
                    z_score = (cluster_mean - cum_mean)/(cum_std/math.sqrt(len(cluster_lifespans)))
                    p_value = stats.norm.sf(abs(z_score))
                    #want to write out gene set, cluster label (best way to denote this?), lifespans.max, lifespans.min, lifespans.mean, lifespans.median
                    writer.writerow([gene_type, cum_mean, f'cluster{label}',round(cluster_mean, 1), round(cluster_std), len(cluster_lifespans), round(z_score, 2), round(1-p_value, 3)])
    
    print("Errors had on the following:")
    for i in error_files: print(i)

    return 

"""
GIVEN AN EMBEDDING FILE, PERFORM PCA AND K-MEANS CLUSTERING ON EMBEDDINGS THEN VISUALIZE OUTPUT w Matplot
"""
def visualize_embeddings(embedding_file, is2d = False):
    embeddings, labels = [],[]
    with open(embedding_file) as to_read:
        for line in to_read:
            line = line.split(",")
            if len(line) < 2: continue #avoid retarded formatting 
            embedding, label = line[:-1], float(line[-1].replace("\n", ""))
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
        else: grph.text(k_means.cluster_centers_[label][0], k_means.cluster_centers_[label][1], f'{round(lifespans.mean())}|{round(numpy.std(lifespans))}', bbox=dict(facecolor='white', alpha=0.5)) #for 2d

    # [grph.text(reduced[i][0], reduced[i][1], reduced[i][2], labels[i]) for i in range(len(labels))] #textually label 
    #try to textually label s.t only the mean is 
    grph.set_title(embedding_file.split("/")[-1][:-4] + " Visualized")
    plt.show()
    fig.savefig("embeddings_visuals/" + embedding_file.split("/")[-1][:-4] + "_visual.png")
    plt.close()


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

