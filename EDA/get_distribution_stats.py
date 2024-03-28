"""
script for producing various stats on the similarity of TOGA orthologs across species (organized by orthologs gene type, if relevant)
motivation is that there is no point in trying to 'learn' on ortholog types for which there is not much variance across species 
METHODS: 
    - get_cumulative_length_distribution():
        - gets length distribution across all TOGA sequences 
            - length is considered as sequence (with all " ", \n, "X" and "-" chars removed)
    
    - get_length_distribution(gene_type_orthologs_path):
        - get length distribution across sequences of this gene type 

    - get_similarity(gene_type_orthologs_path):
        - get measure of similarity across sequences of this gene type (need to figure out measure of similarity that will be good for this)

    - get_pairwise_alignment_similarity(gene_type_orthologs_path):

    - get_full_alignment_similarity(gene_type_orthologs_path):


**remember to download matplot for plotting purposes 

"""
import matplotlib.pyplot as plt
import numpy as np
import sys, os, csv 
from collections import Counter
from Bio.Align import PairwiseAligner, Alignment



local_cumulative_set_path = "/Users/wyattmccarthy/Desktop/MIT Aging Project/Pre-Processing/cumulativeTOGAset.csv"
cumulative_set_untrimmed_path = "/Mounts/rbg-storage1/users/wmccrthy/cumulativeTOGAset.csv"
cumulative_set_trimmed_path = "/Mounts/rbg-storage1/users/wmccrthy/cumulativeTOGAset_trimmed.csv"

"""
GET LENGTH DISTRIBUTION ACROSS ALL TOGA ORTHOLOGS (REGARDLESS OF TYPE)
    - distinguish between max/min length for regulatory and non-regulatory orthologs 
"""
def get_cumulative_length_distribution():
    freqs = {}

    #build set of (regulatory) transcription factors s.t we can differentiate between max/min length for regulatory and non-regulatory orthologs 
    transciption_factor_path = "/Mounts/rbg-storage1/users/wmccrthy/transcription_factors.txt"
    tf_set = set()
    with open(transciption_factor_path) as read_from:
        for line in read_from:
            line = line.split("\t")
            gene_symbol = line[0]
            tf_set.add(gene_symbol)
            
    reg_max_len, reg_min_len, reg_avg = 0, float('inf'), 0
    non_max_len, non_min_len, non_avg = 0, float('inf'), 0
    total_reg = 0
    total_non = 0

    with open(cumulative_set_trimmed_path) as read_from:
        for line in read_from:
            line = line.split(",")

            if line[0] == 'organism' or len(line) < 9: continue 
            seq = line[-1] 
            trimmed_seq = seq.replace("\n", "").replace(" ", "").replace("X", "").replace("-", "")

            # if "-" in trimmed_seq: print(seq)
            gene_symbol = line[1].split(".")[1]

            if gene_symbol in tf_set: #is regulatory ortholog
                length = len(trimmed_seq)
                reg_max_len = max(length, reg_max_len)
                reg_min_len = min(length, reg_min_len)
                reg_avg += length 
                total_reg += 1
            else:
                length = len(trimmed_seq)
                non_max_len = max(length, non_max_len)
                non_min_len = min(length, non_min_len)
                non_avg += length 
                total_non += 1
            
            if length in freqs: freqs[length] += 1
            else: freqs[length] = 1
    reg_avg /= total_reg
    non_avg /= total_non

    #now we have length:# of sequences w length pairs and want to plot histogram 
    fig = plt.figure()
    grph = fig.add_subplot()
    grph.hist(freqs, range=(0, max(reg_max_len, non_max_len)))
    grph.set_title("TOGA Cumulative Sequence Length Distribution")
    grph.set_xlabel("Sequence Length")
    grph.set_ylabel("Frequency")
    plt.show()
    fig.savefig("cumulative_length_dist.png")


    with open("cumulative_seq_lengths.csv", "w") as write_to:
        writer = csv.writer(write_to)
        writer.writerow(['ortholog type', 'max seq length', 'min seq length', 'avg seq length'])
        writer.writerow(['regulatory', reg_max_len, reg_min_len, reg_avg])
        writer.writerow(['non-regulatory', non_max_len, non_min_len, non_avg])
    return 


"""
GET LENGTH DISTRIBUTION FOR A GIVEN GENE TYPE (AT PASSED FILE PATH)
"""
def get_length_distribution(gene_type_dataset_path):
    trimmed_name = gene_type_dataset_path.split("/")[-1][:-4]

    freqs = {}
    max_len = 0
    print("Collecting Length Distribution Stats on", trimmed_name)
    with open(gene_type_dataset_path) as read_from:
        for line in read_from:
            line = line.split(",")
            if line[0] == 'organism' or len(line) < 9: continue 
            seq = line[-1] 
            trimmed_seq = seq.replace("\n", "").replace(" ", "").replace("X", "").replace('"', "")
            length = len(trimmed_seq)
            max_len = max(length, max_len)
            if length in freqs: freqs[length] += 1
            else: freqs[length] = 1
    read_from.close()

    print("Max Sequence Length in", trimmed_name, " = ", max_len, "\n")

    fig = plt.figure()
    grph = fig.add_subplot()
    grph.hist(freqs, range=(0, max_len))
    grph.set_title(f"TOGA {trimmed_name} Sequence Length Distribution")
    grph.set_xlabel("Sequence Length")
    grph.set_ylabel("Frequency")
    # plt.show()
    fig.savefig(f"length_plots/{trimmed_name}_length_dist.png")
    plt.close(fig)
    return max_len

#test path: /Mounts/rbg-storage1/users/wmccrthy/gene_datasets/regulatory/TCF15_orthologs.csv

"""
GET SEQUENCE SIMILARITY (BY COSINE SIMILARITY) ACROSS SEQUENCES OF A CERTAIN GENE TYPE 
"""
def get_similarity(gene_type_dataset_path):
    """
    need to find a good metric of 'similarity' bw strings to use for this, that will also not be too slow 
    options:
        - BoW (and maybe tf-idf) for each sequence 
            - use cosine-similarity on the resultant vectors 
            - consider words as: all 3 codon combos
    """
    trimmed_name = gene_type_dataset_path.split("/")[-1][:-4]
    codons_base = {'---': 0}
    get_all_codons("", codons_base)

    vectorized_seqs = []
    print("Collecting Sequence Similarity Stats on", trimmed_name)

    with open(gene_type_dataset_path) as read_from:
        for line in read_from:
            codons_dist = codons_base.copy()

            line = line.split(",")
            if line[0] == 'organism' or len(line) < 9: continue 
            seq = line[-1] 
            trimmed_seq = seq.replace("\n", "").replace("X", "").replace('"', "").strip()
            trimmed_seq = trimmed_seq.split(" ")
            for codon in trimmed_seq:
                if codon in codons_dist: codons_dist[codon] += 1
            # print(codons_dist)
            vectorized_seqs.append([codons_dist[i] for i in codons_dist])
    read_from.close()

    avg_similarity = [0, 0]
    similarity_freqs = {}
    for i in range(len(vectorized_seqs)):
        for j in range(i+1, len(vectorized_seqs)):
            # print(i, j)
            vecA, vecB = vectorized_seqs[i], vectorized_seqs[j]
            similarity = cosine_similarity(vecA, vecB)
            avg_similarity[0] += similarity
            avg_similarity[1] += 1
            if similarity in similarity_freqs: similarity_freqs[similarity] += 1
            else: similarity_freqs[similarity] = 0

    print("avg similarity across all sequences in ", trimmed_name, " = ", avg_similarity[0]/avg_similarity[1])

    fig = plt.figure()
    grph = fig.add_subplot()
    grph.hist(similarity_freqs, range=(0, 1))
    grph.set_title(f"TOGA {trimmed_name} Sequence Similarity")
    grph.set_xlabel("Similarity")
    grph.set_ylabel("Frequency")
    # plt.show()
    fig.savefig(f"similarity_plots/{trimmed_name}_similarity.png")
    plt.close(fig)
    return avg_similarity[0]/avg_similarity[1]


"""
GET PAIRWISE SEQUENCE ALIGNMENT SIMILARITY ACROSS SEQUENCES OF A GIVEN GENE TYPE (RETURNS AVG ALIGNMENT SIMILARITY BW ALL PAIRWISE ALIGNMENTS)
"""
def get_pairwise_alignment(gene_type_dataset_path):
    trimmed_name = gene_type_dataset_path.split("/")[-1][:-4]
    #step 1: get all sequences from given path 
    seqs = []
    print("Collecting Pairwise Alignment Similarity Stats on", trimmed_name)
    with open(gene_type_dataset_path) as read_from:
        for line in read_from:
            line = line.split(",")
            if line[0] == 'organism' or len(line) < 9: continue 
            seq = line[-1] 
            trimmed_seq = seq.replace("\n", "").replace("X", "").replace('"', "").replace(" ", "")
            if len(trimmed_seq) > 0: seqs.append(trimmed_seq)
    read_from.close()
    
    alignment_scores = {}
    avg_alignment = [0, 0]
    aligner = PairwiseAligner()
    print("iterating thru:", len(seqs), "data points")
    for i in range(len(seqs)):
        print("comparing from:", i)
        for j in range(i+1, len(seqs)):
            seq1, seq2 = seqs[i], seqs[j]
            min_len = min(len(seq1), len(seq2)) #max possible alignment score is if all bases aligned, so divide whatever alignment score is by this value to get percentage of alignment
            sim_score = aligner.score(seq1, seq2)/min_len
            # print(sim_score)
            avg_alignment[0] += sim_score
            avg_alignment[1] += 1
            if sim_score in alignment_scores: alignment_scores[sim_score] += 1
            else: alignment_scores[sim_score] = 1

    print("avg pairwise alignment similarity across all sequences in ", trimmed_name, " = ", avg_alignment[0]/avg_alignment[1])

    fig = plt.figure()
    grph = fig.add_subplot()
    grph.hist(alignment_scores, range=(0, 1))
    grph.set_title(f"TOGA {trimmed_name} Sequence Alignment Similarity")
    grph.set_xlabel("Similarity")
    grph.set_ylabel("Frequency")
    # plt.show()
    fig.savefig(f"alignment_similarity_plots/{trimmed_name}_similarity.png")
    plt.close(fig)
    return avg_alignment[0]/avg_alignment[1]

"""
GET TOTAL SEQUENCE ALIGNMENT SIMILARITY (self-proclaimed 'wyatt similarity measure') ACROSS SEQUENCES OF A GIVEN GENE TYPE (RETURNS ALIGNMENT SCORE FOR ALL SEQUENCES ALIGNED)
"""
def get_full_alignment(gene_type_dataset_path):
    trimmed_name = gene_type_dataset_path.split("/")[-1][:-4]
    #step 1: get all sequences from given path 
    seqs = []
    max_len = 0
    print("Collecting Cumulative Alignment Similarity Stats on", trimmed_name, "\n")
    with open(gene_type_dataset_path) as read_from:
        for line in read_from:
            line = line.split(",")
            if line[0] == 'organism' or len(line) < 9: continue 
            seq = line[-1] 
            trimmed_seq = seq.replace("\n", "").replace("X", "").replace('"', "").replace(" ", "")
            if len(trimmed_seq) > 0: 
                seqs.append(trimmed_seq)
                max_len = max(max_len, len(trimmed_seq))
    read_from.close()
    
    #use max_len to pad all sequences s.t they are all of len == max_len (need equal length to align)
    #pad with '-' char, representing gap 
    for i in range(len(seqs)):
        s = seqs[i]
        if len(s) < max_len: 
            len_diff = max_len - len(s)
            seqs[i] += "-" * len_diff
    
    # total_alignment = Alignment(seqs)
    # print(total_alignment.counts())
    total_alignment = alignment_score(seqs)
    return total_alignment


regulatory_orthologs_path = "/Mounts/rbg-storage1/users/wmccrthy/gene_datasets/regulatory/"
def get_all_similarity():
    """
    FOR EACH FILE IN REGULATORY_PATH:
        get_similarity(file_path)
    """
    avg_sims = {}
    for file in os.listdir(regulatory_orthologs_path):
        if file[-11:] == 'trimmed.csv': #only consider data that's been trimemd (lost, missing sequences and sequences which have too many gaps are removed)
            file_path = regulatory_orthologs_path + file
            avg_file_sim = get_similarity(file_path)
            avg_sims[file] = avg_file_sim
    
    with open("avg_seq_similarity.csv", "w") as write_to:
        writer = csv.writer(write_to)
        writer.writerow(['gene type', 'cosine similarity'])
        for gene_type_file in avg_sims:
            gene_type = gene_type_file[:-4]
            writer.writerow[[gene_type, avg_sims[gene_type_file]]]
        

def get_all_len_dist():
    """
    FOR EACH FILE IN REGULATORY_PATH:
        get_length_distribution(file_path)
    """
    max_lens = {}
    for file in os.listdir(regulatory_orthologs_path):
        if file[-11:] == 'trimmed.csv': #only consider data that's been trimemd (lost, missing sequences and sequences which have too many gaps are removed)
            file_path = regulatory_orthologs_path + file
            max_len = get_length_distribution(file_path)
            max_lens[file] = max_len

    with open("max_seq_lengths.csv", "w") as write_to:
        writer = csv.writer(write_to)
        writer.writerow(['gene type', 'max seq length'])
        for gene_type_file in max_lens:
            gene_type = gene_type_file[:-4]
            max_len = max_lens[gene_type_file]
            writer.writerow([gene_type, max_len])


def get_all_pairwise_alignments():
    avg_alignments = {}
    for file in os.listdir(regulatory_orthologs_path):
        if file[-11:] == 'trimmed.csv': #only consider data that's been trimemd (lost, missing sequences and sequences which have too many gaps are removed)
            file_path = regulatory_orthologs_path + file
            avg_alignments[file] = get_pairwise_alignment(file_path)

    with open("average_pairwise_alignment_similarities.csv", "w") as write_to:
        writer = csv.writer(write_to)
        writer.writerow(['gene type', 'avg alignment similarity'])
        for gene_type_file in avg_alignments:
            gene_type = gene_type_file[:-4]
            avg_sim = avg_alignments[gene_type_file]
            writer.writerow([gene_type, avg_sim])

def get_all_full_alignments():
    w_alignments = {}
    for file in os.listdir(regulatory_orthologs_path):
        if file[-11:] == 'trimmed.csv': #only consider data that's been trimemd (lost, missing sequences and sequences which have too many gaps are removed)
            file_path = regulatory_orthologs_path + file
            w_alignments[file] =  get_full_alignment(file_path)

    with open("average_wyatt_alignment_similarities.csv", "w") as write_to:
        writer = csv.writer(write_to)
        writer.writerow(['gene type', 'avg similarity (where 0=max variability, .375=zero variability)'])
        for gene_type_file in w_alignments:
            gene_type = gene_type_file[:-4]
            avg_sim = w_alignments[gene_type_file]
            writer.writerow([gene_type, avg_sim])






"""
HELPER METHODS
"""
def euclidean_norm(vec):
    return np.linalg.norm(vec)

def cosine_similarity(vecA, vecB):
    # print(vecA, vecB)
    dot_prod = np.dot(vecA, vecB)
    euclidean_norm_product = (euclidean_norm(vecA) * euclidean_norm(vecB))

    if dot_prod == 0 or euclidean_norm_product == 0: return .5 #if error return unbiased .5 value 

    return dot_prod/euclidean_norm_product
    

def get_all_codons(cur_codon, vocab):
    if len(cur_codon) == 3: 
        vocab[cur_codon] = 0
        return 

    get_all_codons(cur_codon + 'A', vocab)
    get_all_codons(cur_codon + 'G', vocab)
    get_all_codons(cur_codon + 'T', vocab)
    get_all_codons(cur_codon + 'C', vocab)

def pad_bases_dict(dct):
    if 'N' not in dct: dct['N'] = 0
    if 'A' not in dct: dct['A'] = 0
    if 'T' not in dct: dct['T'] = 0
    if 'G' not in dct: dct['G'] = 0
    if 'C' not in dct: dct['C'] = 0

"""
METHOD THAT ITERATES OVER COLUMNS OF BASICALLY ALIGNED SEQUENCES (list of sequences of same gene type, all padded to same length and all starting w 'ATG') 
AND COMPUTES DISTRIBUTION AT EACH POSITION 
RETURNED SCORE IS AVERAGE OF DISTRIBUTION ACROSS ALL POSITIONS 
"""
def alignment_score(sequences):
    total_avg_dist = 0
    for i in range(len(sequences[0])):
        #for each column (aligned pos), get dist of bases 
        dist = {'A': 0, 'T':0, 'C':0, 'G': 0}
        total = 0
        for j in range(len(sequences)):
            if sequences[j][i] in dist: 
                dist[sequences[j][i]] += 1
                total += 1
        if total < 50: #too small a sample size from which to judge 
            # print(i, j, len(sequences))
            continue 
            
        # print(dist)
        dist = {i:dist[i]/total for i in dist}

        # if i < 3: print(dist) #testing to ensure that (for the most part) sequences begin w ATG, s.t aligning them from their start coordinates is valid 
        
        #perfectly dissimilar distribution wud see A, T, G, and C each having dist frequency of .25 here, thus, to garner an idea of similarity (where more similar indicates certain bases appear much more freq than others)
        #use a measure of avg distance from perfectly even distribution (.25)
        avg_dist = sum([abs(i - 0.25) for i in dist.values()])/4
        total_avg_dist += avg_dist
        # print("Average Distance from Perfectly Dissimilar Base Distribution:", avg_dist, dist, "Column", i, "\n")

    # print("total avg dist from dissimilarity:", total_avg_dist/len(sequences[0])) 

    #keep in mind range of 'score' here is from 0 - .375 where .375 corresponds to NO VARIABILITY  
    return total_avg_dist/len(sequences[0])
        


if __name__ == "__main__":
    args = sys.argv
    print(args)
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])