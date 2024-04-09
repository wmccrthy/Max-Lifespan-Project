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
from transformers import AutoTokenizer


"""
GIVEN A FILE PATH TO A GENE SET, ALIGNS SEQUENCES THEN COMPUTES CHI SQUARE PER POSITION IN SEQUENCE ALIGNMENT (DERIVED FROM DISTRIBUTION OF NUCLEOTIDES AT EACH POSITION)
SET_CHI_SQUARE = SUM(CHI SQUARE) FOR EACH POSITION IN ALIGNMENT 
COMPUTES CRAMER'S V FOR THE SET GIVEN SET_CHI_SQUARE 

ALSO PERFORMS ANALYSIS ON VARIOUS 'BINS' (ORGANIZED BY LIFESPAN RANGES) OF THE SET
"""
def analyze(gene_type_dataset_path, get_bins=False, chi_squared=False):
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

    
    if chi_squared and not get_bins: 
        return [alignment_score(seqs, True)]
    elif not chi_squared: total_alignment_score = alignment_score(seqs)
    else: parent_chi_data = [alignment_score(seqs, True)]

    # print(trimmed_name, " Cumulative Stats\n", "Max Len:", max_len, "Min Len:", min_len, "Mean Len:", avg_len, "Median Len:", median_len, "\nAlignment Similarity:", total_alignment_score, "\n")
    if not chi_squared: to_ret = [max_len, min_len, avg_len, median_len, round(total_alignment_score, 5), total_species]
    
    if not get_bins: return to_ret 

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
        
        if not chi_squared: bucket_alignment_score = alignment_score(lifespan_buckets[bucket][0])
        else: 
            bucket_chi_square_data = [alignment_score(lifespan_buckets[bucket][0], True)]
            #parent set cramer V, bucket, bucket cramer v, bucket cramer v - parent set cramer v, # species in bucket
            bins_of_interest.append([i for i in parent_chi_data] + [bucket] + [i for i in bucket_chi_square_data] + [bucket_chi_square_data[0] - parent_chi_data[0]] + [len(lifespan_buckets[bucket][1]), len(lifespan_buckets[bucket][0])])
            continue 
        # print(len(lifespan_buckets[bucket]))
        # print(trimmed_name, bucket, "Stats\n", "Max Len:", max_len, "Min Len:", min_len, "Mean Len:", avg_len, "Median Len:", median_len, "\nAlignment Similarity:", bucket_alignment_score, "\n")

        # if abs(total_alignment_score - bucket_alignment_score) > .05 and bucket_alignment_score > 0:
            # print("NOTABLE SIMILARITY DIFFERENCE")
            # print("Number of Seqs in Bucket:", len(lifespan_buckets[bucket][0]))
            # print("Number of Species in Bucket:", len(lifespan_buckets[bucket][1]))
            # print(trimmed_name, bucket, "Stats\n", "Max Len:", max_len, "Min Len:", min_len, "Mean Len:", round(avg_len, 3), "Median Len:", median_len, "\nAlignment Similarity:", round(bucket_alignment_score, 5), "vs. All Buckets Alignment Similarity", round(total_alignment_score, 5), "\n")
            
            # bins_of_interest.append([bucket, bucket_alignment_score, len(lifespan_buckets[bucket][1]), bucket_alignment_score - total_alignment_score])
            #PERHAPS STORE THIS INFORMATION GIVEN ITS RELEVANCE (COULD HAVE SET OF TYPES/BUCKETS with Highest/Lowest Similarity)

        bins_of_interest.append([bucket, bucket_alignment_score, len(lifespan_buckets[bucket][1]), len(lifespan_buckets[bucket][0]), bucket_alignment_score - total_alignment_score]) #append bucket to get all bucket data, regardless if different from that of parent set
    if not chi_squared and get_bins: return [total_alignment_score, bins_of_interest]
    elif get_bins: return bins_of_interest

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

gene_datasets_path = "/data/rsg/chemistry/wmccrthy/Everything/gene_datasets/regulatory/"



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
        writer.writerow(['gene type', 'lifespan bin', 'set similarity', 'bin similarity', '# species in bin', '# sequences in bin', 'similarity diff'])
        for file in os.listdir(gene_datasets_path):
            if file[-11:] == "trimmed.csv":
                gene_type = file.split("_")[0]
                gene_type_alignment, gene_type_bins = analyze(gene_datasets_path + "/" + file, get_bins=True)
                for bin in gene_type_bins:
                    writer.writerow([gene_type, bin[0], gene_type_alignment, bin[1], bin[2], bin[3], bin[4]])
    return 



"""
ITERATE THRU EACH REGULATORY GENE SET:
    - COMPUTE CHI SQUARED VALUE FOR THE SET 
    - OUTPUT [GENE TYPE, CHI SQUARED VALUE, DF] 

"""
def analyze_chi_squared():
    with open("regulatory_sets_chi_square.csv", "w") as write_to:
        writer = csv.writer(write_to)
        writer.writerow(['gene type', 'normalized chi squared', 'bonferroni-corrected 5% threshold', 'cramers v'])
        for file in os.listdir(gene_datasets_path):
            if file[-11:] == "trimmed.csv":
                gene_type = file.split("_")[0]
                gene_type_data = analyze(gene_datasets_path + "/" + file, chi_squared = True)
                writer.writerow([gene_type] + gene_type_data)

def analyze_chi_squared_bins():
    with open("regulatory_bins_chi_square.csv", "w") as write_to:
        writer = csv.writer(write_to)
        writer.writerow(['gene type', 'gene type cramer v', 'lifespan bin', 'bin cramer v', 'cramer v difference', '#species in bin', '#sequences in bin'])
        for file in os.listdir(gene_datasets_path):
            if file[-11:] == "trimmed.csv":
                gene_type = file.split("_")[0]
                gene_type_bins = analyze(gene_datasets_path + "/" + file, get_bins=True, chi_squared = True)
                for bin in gene_type_bins:
                    writer.writerow([gene_type] + [i for i in bin])


"""
TOKENIZES A TF SET OF SEQUENCES PRIOR TO DNABERT_S embeddings 
"""
def analyze_tokenized_seq_lengths():
    DNABERT_S_path = "/data/rsg/chemistry/wmccrthy/Everything/DNABERT-S/"
    tokenizer = AutoTokenizer.from_pretrained(DNABERT_S_path, trust_remote_code=True)
    with open("regulatory_sets_tokenized_lengths.csv", "w") as write_to:
        writer = csv.writer(write_to)
        writer.writerow(['gene type', 'max tokenized length', 'min tokenized length'])
        for file in os.listdir(gene_datasets_path):
            if file[-11:] == "trimmed.csv":
                gene_type = file.split("_")[0]
                file_path = gene_datasets_path + file
                print(file)
                with open(file_path) as to_read:
                    min_len, max_len = float('inf'), 0
                    for line in to_read:
                        line = line.split(",")
                        if len(line) < 9 or line[0] == 'organism': continue 
                        seq = line[-1]
                        seq = seq.replace('"', "").replace(" ", "").replace("\n", "").replace("X", "")

                        if len(seq) == 0: continue #odd outlying data 

                        seq_tokenized = tokenizer(seq, return_tensors='pt')['input_ids']
                        # print(seq.shape)
                        if len(seq_tokenized[0]) == 10: print(line)
                        min_len, max_len = min(len(seq_tokenized[0]), min_len), max(len(seq_tokenized[0]), max_len)
                print(min_len, max_len, "\n")
                writer.writerow([gene_type, max_len, min_len])

"""
PRINTS METADATA ON THE # OF GENE SETS W MAX (TOKENIZED) SEQUENCE LENGTHS AT VARIOUS THRESHOLDS 
"""
def count_tokenized_lengths():
    freqs = {0:0, 1:0, 2:0, 3:0}
    total_sets = 0
    with open("/Users/wyattmccarthy/Desktop/MIT Aging Project/Everything/EDA/regulatory_sets_tokenized_lengths.csv") as read_from:
        for line in read_from:
            line = line.split(',')
            if line[0] == 'gene type': continue 
            total_sets += 1
            line[1] = int(line[1])
            if line[1] > 500: freqs[3] += 1
            if line[1] > 1000: freqs[0] += 1
            if line[1] > 2000: freqs[1] += 1
            if line[1] > 5000: freqs[2] += 1
    print('# Gene Sets:', total_sets)
    print('# Sets w Max Seq Length > 500:', freqs[3])
    print('# Sets w Max Seq Length > 1000:', freqs[0])
    print('# Sets w Max Seq Length > 2000:', freqs[1])
    print('# Sets w Max Seq Length > 5000:', freqs[2])

if __name__ == "__main__":
    args = sys.argv
    print(args)
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])



"""
ATG CCG CGC TCC TTC TTG GTA AAG AAG ATC AAA GGG GAC GGC TTC CAG TGC AGC GGG GTG CCG GCC CCC ACC TAC CAC CCC TTG GAG (capensis ratel, 31.3)

ATG CCG CGC TCC TTC CTG GTA AAG AAG ATC AAA GGG GAC GGC TTC CAG TGC AGC GGG GTG CCG GCC CCC ACC TAC CAC CCC TTG GAG ACC GCC TAC GTG (harbor seal, 47.6)

ATG CCG CGC TCC TTC CTG GTA AAG AAG ATC AAA GGG GAC GGC TTC CAG TGC AGC GGG GTG CCG GCC CCC ACC TAC CAC CCT TTG GAG ACC GCC TAT (hyena, 41.1)
^all three above in range 30-60 lifespan


ATG CCC CGC TCC TTC TTG GTG AAG AAG ATC AAA GCA GAT GCC TTC CAG TGT AGC AGC GTC CCA GCT CCC AGC TAC CAC CCC CTG GGA TCC GCA TAT GTG CTG (yellow-footed antechinus, 3.9)

ATG CCC CGC TCC TTC CTG GTG AAG AAG ATC AAA GCA GAC GCC TTC CCG TGT AGC AGC GTC CCA GCT CCT AGC TAC CAC CCC CTG GGA TCC GCA TAT GTG TTG (agile gracile mouse opossum, 6.1)

ATG CCG CGC TCC TTC CTG GTC AAG AAG ATC AAA GGC GAC GGC TTC CAG TGC AGC GGG GTG CCG GCC CCC ACC TAC CAC CCC TTG GAG ACC GCC TAC GTG (pygmy shrew, 3.2)
(0-15 range)
"""




"""
DNAJC1


ATG ACG GCT CCT TGC TCC CAG CCG GCG CAG CTT CCT GGA CGC CGC CAG CTC AGG CTG GTG CCG TTC some typa gorilla

ATG ACG GCG CCC TGC TCC CGG CTG GCC CGG CTT CCT CCG CGC CGC CGG CTC CGG CTC GTG CCG TTC some typa whale 

ATG ACG GCG CCC TGC TGC CGG CTG GCC CGG CTT CCT CCG CGC CGC CGG CTC CGG CTA CGG CTC GTG another typa whale 
(60-100)

ATG ACG GCT CCC TCG TCC CCG GCC GCC CCG CGG CAG CCT CCG CCG CCG CTC CGC CTG CTG AGG TGG shrew 
(0-15)
"""