"""
SCRIPT FOR PERFORMING BASIC ANALYSIS ON ARBITRARILY SELECTED REGULATORY ORTHOLOG SETS; 
IDEA IS THAT CONFINING ANALYSIS TO SINGLE GENE TYPES MIGHT MAKE GLEANING FEATURES RELEVANT TO OUR MODEL A SIMPLER TASK / identify gene types which are relevant to look at 

COLLECT CUMULATIVE DATA ON: 
    - SIZE (MIN, MAX, MEDIAN, MEAN)
    - 'Wyatt' Alignment Score 

CREATE DICT OF LIFESPAN BUCKETS TO KEEP TRACK OF LIST OF SEQUENCES FOR EACH BUCKET
USE LIFESPAN BUCKETS DICT TO PERFORM ANALYSIS ON EACH BUCKET
IDEAS: 
    - IS THERE MORE SIMILARITY BETWEEN SEQUENCES IN SIMILAR BUCKETS? 
        - IN TERMS OF SIZE 
        - IN TERMS OF 'WYATT' ALIGNMENT SCORE 


"""
from get_distribution_stats import alignment_score
import statistics, sys, os, csv
from Bio import motifs 
from Bio.Seq import Seq 

def analyze(gene_type_dataset_path, get_bins=False):
    trimmed_name = gene_type_dataset_path.split("/")[-1][:-4]

    print("Performing Analysis on", trimmed_name)

    lifespan_buckets = {(0, 15):([], set()), (15, 30):([], set()), (30, 60):([], set()), (60, 100):([], set()), (100, 250):([], set())} #these were arbitrarily estimated 

    #CUMULATIVE SET ANALYSIS 
    seqs = []
    represented_species = set()
    with open(gene_type_dataset_path) as read_from:
        for line in read_from:
            line = line.split(",")
            if line[0] == 'organism' or len(line) < 9: continue 

            seq = line[-1] 
            lifespan = line[1]
            species = "_".join(line[0].split("_")[:2])
            represented_species.add(species)

            trimmed_seq = seq.replace("\n", "").replace("X", "").replace('"', "").replace(" ", "")

            if len(lifespan) > 5 or len(trimmed_seq) == 0: continue #indicates lifespan is either unknown or not established 
            else: lifespan = float(lifespan)

            for bucket in lifespan_buckets:
                if lifespan < bucket[1] and lifespan > bucket[0]:
                    lifespan_buckets[bucket][0].append(trimmed_seq)
                    lifespan_buckets[bucket][1].add(species)
                    continue 
            seqs.append(trimmed_seq)
    read_from.close()

    seq_lengths = sorted([len(s) for s in seqs])
    max_len, min_len, avg_len, median_len  = max(seq_lengths), min(seq_lengths), statistics.mean(seq_lengths), statistics.median(seq_lengths)
    total_species = len(represented_species)

    #use max_len to pad all sequences s.t they are all of len == max_len (need equal length to align)
    #pad with '-' char, representing gap 
    for i in range(len(seqs)):
        s = seqs[i]
        if len(s) < max_len: 
            len_diff = max_len - len(s)
            seqs[i] += "-" * len_diff

    total_alignment_score = alignment_score(seqs)
    # print(trimmed_name, " Cumulative Stats\n", "Max Len:", max_len, "Min Len:", min_len, "Mean Len:", avg_len, "Median Len:", median_len, "\nAlignment Similarity:", total_alignment_score, "\n")
    to_ret = [max_len, min_len, avg_len, median_len, round(total_alignment_score, 5), total_species]
    
    # return to_ret 

    #waive for now (until we can confirm validity of chosen buckets)
    #PER LIFESPAN BUCKET ANALYSIS
    #essentially performs the same analysis done on the cumulative set to each lifespan bin subset 

    bins_of_interest = [] #use list to store bins for which there is notable similarity difference between bin sequences and full set sequences 

    for bucket in lifespan_buckets:
        if len(lifespan_buckets[bucket][0]) == 0: continue 

        bucket_lengths = sorted([len(s) for s in lifespan_buckets[bucket][0]])
        
        max_len, min_len, avg_len, median_len = max(bucket_lengths), min(bucket_lengths), statistics.mean(bucket_lengths), statistics.median(bucket_lengths)
        #use max_len to pad all sequences s.t they are all of len == max_len (need equal length to align)
        #pad with '-' char, representing gap 
        for i in range(len(lifespan_buckets[bucket][0])):
            s = lifespan_buckets[bucket][0][i]
            if len(s) < max_len: 
                len_diff = max_len - len(s)
                lifespan_buckets[bucket][0][i] += "-" * len_diff
        
        bucket_alignment_score = alignment_score(lifespan_buckets[bucket][0])
        # print(len(lifespan_buckets[bucket]))
        # print(trimmed_name, bucket, "Stats\n", "Max Len:", max_len, "Min Len:", min_len, "Mean Len:", avg_len, "Median Len:", median_len, "\nAlignment Similarity:", bucket_alignment_score, "\n")


        if abs(total_alignment_score - bucket_alignment_score) > .05 and bucket_alignment_score > 0:

            # print("NOTABLE SIMILARITY DIFFERENCE")
            # print("Number of Seqs in Bucket:", len(lifespan_buckets[bucket][0]))
            # print("Number of Species in Bucket:", len(lifespan_buckets[bucket][1]))
            # print(trimmed_name, bucket, "Stats\n", "Max Len:", max_len, "Min Len:", min_len, "Mean Len:", round(avg_len, 3), "Median Len:", median_len, "\nAlignment Similarity:", round(bucket_alignment_score, 5), "vs. All Buckets Alignment Similarity", round(total_alignment_score, 5), "\n")
            
            bins_of_interest.append([bucket, bucket_alignment_score])
            #PERHAPS STORE THIS INFORMATION GIVEN ITS RELEVANCE (COULD HAVE SET OF TYPES/BUCKETS with Highest/Lowest Similarity)
    
    if get_bins: return [total_alignment_score, bins_of_interest]

    return to_ret

"""
CREATE WEBLOGOS FOR EACH GENE TYPE and lifespan buckets within each gene type TO VISUALIZE NUCLEOTIDE DISTRIBUTION AND POTENTIAL PRESENCE OF MOTIFS 
"""
def weblogo_analysis():
    for file in os.listdir(gene_datasets_path):
        gene_type = file.split("_")[0]
        if file[-11:] == "trimmed.csv": 
            print("generating motif plots for", gene_type, "orthologs")
            seqs = []
            max_len = 0
            lifespan_buckets = {(0, 15):[], (15, 30):[], (30, 60):[], (60, 100):[], (100, 250):[]} #these were arbitrarily estimated 


            with open(gene_datasets_path + '/' + file) as read_from:
                for line in read_from:
                    line = line.split(',')
                    if line[0] == 'organism' or len(line) < 9: continue 
                    seq = line[-1] 
                    lifespan = line[1]
                    trimmed_seq = seq.replace("\n", "").replace("X", "").replace('"', "").replace(" ", "")
                    if len(trimmed_seq) > 0: 
                        max_len = max(max_len, len(trimmed_seq))
                        seqs.append(trimmed_seq)

                    if len(lifespan) > 5 or len(trimmed_seq) == 0: continue #indicates lifespan is either unknown or not established 
                    else: lifespan = float(lifespan)
                    for bucket in lifespan_buckets:
                        if lifespan < bucket[1] and lifespan > bucket[0]:
                            lifespan_buckets[bucket].append(trimmed_seq)
                            continue 
            #pad sequences to b of same length 
            for i in range(len(seqs)):
                s = seqs[i]
                if len(s) < max_len: 
                    len_diff = max_len - len(s)
                    seqs[i] += "-" * len_diff

            #motif plot for entire gene-type set 
            os.system(f"mkdir gene_motif_plots/{gene_type}")
            motif_format = motifs.create([Seq(i) for i in seqs])
            motif_format.weblogo(f"gene_motif_plots/{gene_type}/general_motif_plot.png")
            #motif plots for each lifespan bin 
            for bucket in lifespan_buckets:
                bucket_max_len = max([len(i) for i in lifespan_buckets[bucket]])
                #pad sequences in bucket 
                for i in range(len(lifespan_buckets[bucket])):
                    s = lifespan_buckets[bucket][i]
                    if len(s) < bucket_max_len: 
                        len_diff = bucket_max_len - len(s)
                        lifespan_buckets[bucket][i] += "-" * len_diff
                bucket_motif = motifs.create([Seq(i) for i in lifespan_buckets[bucket]])
                bucket_motif.weblogo(f"gene_motif_plots/{gene_type}/{bucket}_motif_plot.png")
    return 

gene_datasets_path = "/Mounts/rbg-storage1/users/wmccrthy/gene_datasets/regulatory/"



"""
ITERATE THRU EACH REGULATORY GENE SET:
    - CALL ANALYSIS METHOD (RETURNS SEQUENCE LENGTH DISTRIBUTION DATA ACROSS SET)
        - ALSO PRINTS OUT INTERESTING SIMILARITY DATA PERTAINING TO LIFESPAN BUCKETS WITHIN ENTIRE SET 
    - WRITE OUT LENGTH DISTRIBUTION TO OUTPUT FILE 
"""
def analyze_multiple():
    with open("regulatory_sets_metadata.csv", "w") as write_to:
        writer = csv.writer(write_to)
        writer.writerow(['gene type', 'max seq len', 'min seq len', 'mean seq len', 'median seq len', 'wmc similarity score', 'species represented'])
        for file in os.listdir(gene_datasets_path):
            if file[-11:] == "trimmed.csv":
                gene_type = file.split("_")[0]
                gene_type_data = analyze(gene_datasets_path + "/" + file)
                writer.writerow([gene_type] + gene_type_data)
                # print()
    return 

"""
ITERATE THRU EACH REGULATORY GENE SET:
    - CALL ANALYSIS METHOD w True flag (RETURNS bins within set that have differing similarity to set as a whole)
    - WRITE OUT interesting bins to output csv 
"""
def analyze_multiple_bins():
    with open("regulatory_interesting_bins.csv", "w") as write_to:
        writer = csv.writer(write_to)
        writer.writerow(['gene type', 'lifespan bin', 'set similarity', 'bin similarity'])
        for file in os.listdir(gene_datasets_path):
            if file[-11:] == "trimmed.csv":
                gene_type = file.split("_")[0]
                gene_type_alignment, gene_type_bins = analyze(gene_datasets_path + "/" + file, True)
                for bin in gene_type_bins:
                    writer.writerow([gene_type, bin[0], gene_type_alignment, bin[1]])
    return 


if __name__ == "__main__":
    args = sys.argv
    print(args)
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])


