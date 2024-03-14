import os, csv, sys 
from Bio.Seq import Seq 
from Bio import motifs 

"""
SCRIPT FOR GETTING MOTIF ANNOTATIONS ACROSS ALL DATA 
ESSENTIALLY THE CUMULATIVE VERSION OF 'motif_annotation.py', JUST EXTENDS THE SAME PROCEDURE TO OPERATE ACROSS ALL FILES RATHER 
THAN SEARCH FOR SPECIFIC SPECIES WITHIN FILES 

PSUEDO / CONCEPT: 
    - takes in motif_instances parameter 
        - motif = motifs.create([Seq(i) for i in instances]) #create motif based on instances input (convert to Seq s.t input supports str or Seq type)
        - pwm = motif.counts.normalize(pseudocounts=0.5)
        - pssm = pwm.log_odds() 


    - USE dict TO STORE SPECIES:MOTIF_DATA PAIRS 
        - NEED TO THINK ABT WHAT DATA PERTAINING TO MOTIFS WE WANT 
    - motif_data = {}
        - motif_data[organism] = [{}, {}, [], set()]
        - motif_data[organism][0] = frequency dict for instances of motif found (only considering instances of motif found closest to assumed start-codon/transcript-start position)
        - motif_data[organism][1] = frequency dict for similar motifs found (only considering similar motifs found closer to assumed start-codon than closest exact motif)
        - motif_data[organism][2] = list of positions at which motif was found (relative to start codon)
            - from these we can derive average and mode 
        - motif_data[organism][3] = set of (chrom, start, end) checked as there r duplicates in data 

    - FOR FILE IN OS.LISTDIR(DATA_SUBSETS_path): 
        file_path = DATA_SUBSETS_path + file 
        with open(file_path) as read_from:
            for line in read_from:
                line = line.split(",)

                if line[0] == 'organism': continue #skip over header row
                
                organism, sequence = line[0], line[-1][:1005]
                only extract the portion of the sequence we care about 

                if organism not in motif_data: motif_data[organism] = [{}, {}, []]

                SEARCH FOR MOTIF IN SEQUENCE 
                    record frequencies of instance of motif found closest to start codon (motif_data[organism][0]) 
                    
                    record positions of instance of motif found closest to start codon (motif-data[organism][2])
                    

                SEARCH FOR SIMILAR MOTIF INSTANCES EVEN IF EXACT MOTIF ISN'T FOUND      
                    look for highest similarity-scoring pos/sequence that is closer than exact motif matches found 
                    
                    record frequencies of sequences most similar to motif (if closer than actual found motif) (motif_data[organism][1])

  HAVING FILLED DICTIONARY WITH ORGANISM MOTIF DATA, ITERATE THROUGH TOGA FILE STRUCTURE AND WRITE MOTIF DATA TO ORGANISM'S DIRECTORY 
    - what is a good way to store such data? 
        - frequency and distribution based 
        - CSV of form: motif_instance | frequency
            - have two csv's per organism: 1 for exact motif matches, 1 for similar motifs 
    - brainstorm what other data we might want to annotate and output 

"""       
def get_start_codon(organism = None):
    start_codon_instances = [Seq("ATG")]
    get_motif_annotations(start_codon_instances, organism, "start-codon")

def get_kozak_annotations(organism = None):
    canonical_kozak = "GCC(A/G)CCATGG"
    kozak_instances = [Seq("GCCACCATGG"), Seq("GCCGCCATGG")]
    get_motif_annotations(kozak_instances, organism, "kozak")

def get_tata_annotations(organism = None):
    tata_instances = []
    def get_TATA_variants(variant):
        if len(variant) >= len(canonical_tata): 
            tata_instances.append(variant)
            return 
        curInd = len(variant)
        if canonical_tata[curInd] == 'W': 
            get_TATA_variants(variant + "A")
            get_TATA_variants(variant + "T")
        elif canonical_tata[curInd] == 'N':
            get_TATA_variants(variant + "A")
            get_TATA_variants(variant + "T")
            get_TATA_variants(variant + "C")
            get_TATA_variants(variant + "G")
        else: get_TATA_variants(variant + canonical_tata[curInd])
    canonical_tata = "TATAWAWN"
    get_TATA_variants("TATA")
    tata_instances = [Seq(tata_variant) for tata_variant in tata_instances]
    get_motif_annotations(tata_instances, organism, "tata")  

def get_caat_annotations(organism = None):
    caat_instances = []
    canonical_caat = "GG(T/C)CAATCT"
    caat_instances = [Seq("GGCCAATCT"), Seq("GGTCAATCT")]
    get_motif_annotations(caat_instances, organism, "caat")

#what other conserved regions can we consider annotating? 

def get_motif_annotations(instances, organism_of_interest, motif_name):
    data_evaluated = 0 #variable used to track progress / collect stats to estimate total runtime 

    data_subsets_path = f"/Mounts/rbg-storage1/users/wmccrthy/dataSubSets/"

    motif = motifs.create([Seq(i) for i in instances]) #create motif based on instances input (convert to Seq s.t input supports str or Seq type)

    pwm = motif.counts.normalize(pseudocounts=0.5)
    # do we want to change pseudocounts based on general distribution of each sequence? 
    # if so, can update pwm (and pssm) for each iteration of inner-most loop based on each bases probability to appear in promoter sequence 

    pssm = pwm.log_odds() 

    motif_data = {}

    for file in os.listdir(data_subsets_path):
        file_path = data_subsets_path + file
        with open(file_path) as read_from:
            for line in read_from:
                line = line.split(",")

                if line[0] == 'organism': continue #skip over header line 

                organism, sequence = line[0].strip(), line[-1][:1005]
                #only extract the portion of the sequence we care about 

                if len(organism) < 10: continue #error data 

                #if organism given as input, only get data for that species 
                if organism_of_interest: 
                    if organism != organism_of_interest: continue

                if organism not in motif_data: motif_data[organism] = [{}, {}, [], set()] #ensure organism is in our motif_data dict s.t we can populate relevant values 

                chrom, start, end  = line[4:7]
                if (chrom, start, end) in motif_data[organism][3]: continue 
                else: motif_data[organism][3].add((chrom, start, end)) 

                #format sequence (get rid of spaces, and "" marks)
                sequence = sequence.replace('"', '').replace(" ", "").strip().upper()
                
                #preliminary fix for duplicate 'A' produced by querying 1000 BPs upstream from start rather than from start-1
                if sequence[-7:-2] == 'AAATG': sequence = sequence[:-7] + sequence[-6:]
                
                sequence = Seq(sequence)
                start_codon = 1000 
                #start codon generally assumed to be at position 1000; perhaps add functionality s.t if seq at 1000 is not ATG, we make note s.t data is less 'credible'
                #perhaps just have printed output if start codon not present that motif was found x bases from 'assumed' start codon position 

                data_evaluated += 1 
                print(data_evaluated)


                #SEARCH FOR MOTIF IN SEQUENCE 
                closest_motif = [None, None, 9999] #variable for storing instance of motif found closest to start codon
                if len(list(sequence.search(motif.alignment))) > 0:
                    for pos, seq in sequence.search(motif.alignment): 
                        distance_from_start_codon = start_codon - pos
                        if distance_from_start_codon < closest_motif[2]: 
                            closest_motif = [seq, pos, distance_from_start_codon]
                    
                    #record frequencies of instance of motif found closest to start codon (motif_data[organism][0]) 
                    if str(closest_motif[0]) not in motif_data[organism][0]: motif_data[organism][0][str(closest_motif[0])] = 1
                    else: motif_data[organism][0][str(closest_motif[0])] += 1
                    motif_data[organism][2].append(closest_motif[2])
                    #record positions of instance of motif found closest to start codon (motif-data[organism][2])
                    
                #SEARCH FOR SIMILAR MOTIF INSTANCES EVEN IF EXACT MOTIF ISN'T FOUND
                #little sus about capturing accurate data when pos is negative (indicating it pertains to sequence on the reverse strand)
                highest_scored_pos, most_similar_seq = 0, None 
                for pos, score in pssm.search(sequence, threshold= 7.5):
                    similar_seq = sequence[pos:pos+motif.length] # if pos < 0: similar_seq = sequence[pos:pos+motif.length].reverse_complement() ? 
                    distance_from_start_codon = start_codon - abs(pos)  
                    #look for highest similarity-scoring pos/sequence that is closer than exact motif matches found 
                    if score > highest_scored_pos and distance_from_start_codon < closest_motif[2]: 
                        most_similar_seq = [similar_seq, pos, distance_from_start_codon]
                        highest_scored_pos = score

                # record frequencies of sequences most similar to motif (if closer than actual found motif) (motif_data[organism][1])
                if most_similar_seq:
                    if str(most_similar_seq[0]) not in motif_data[organism][1]: motif_data[organism][1][str(most_similar_seq[0])] = 1
                    else: motif_data[organism][1][str(most_similar_seq[0])] += 1
                    
                # print(organism, "Closest Motif found: ", str(closest_motif[0]), "| at pos", str(closest_motif[1]), "|", str(closest_motif[2]), "base pairs from start codon")
                # if most_similar_seq: print(organism, "Most Similar (not exact) Sequence found: ", str(most_similar_seq[0]), "| at pos", str(most_similar_seq[1]), "|", str(most_similar_seq[2]), "base pairs from start codon")
                # print(sorted([(i, motif_data[organism][0][i]) for i in motif_data[organism][0]], key = lambda x:x[1]))
                # print(sorted([(i, motif_data[organism][1][i]) for i in motif_data[organism][1]], key = lambda x:x[1]))
                # print()


    #having iterated thru entirety of data, I want to now output csv files indicating motif instances found and their frequencies 
    
    #iterate thru TOGA file structure and at each organism, iterate thru motif_data[organism] writing data out as we go 

    hg_38_path = "/Mounts/rbg-storage1/users/wmccrthy/human_hg38_reference/"
    for order in os.listdir(hg_38_path):
        order_dir = f"{hg_38_path}{order}/"
        for organism in os.listdir(order_dir):

        
            exact_motif_output = f"{order_dir}{organism}/exact_{motif_name}_instances.csv"
            similar_motif_output = f"{order_dir}{organism}/similar_{motif_name}_instances.csv"

            #if files already exists, continue, no need to redundantly write data  
                # if os.system(f'test -e {exact_motif_output}') == 0 and os.system(f'test -e {similar_motif_output}') == 0: 
                #     print("file exists, skipping")
                #     continue
            
            
            if organism not in motif_data: 
                continue #avoid organisms for which data is missing (also thus only pertains to input organism, if passed) 
            else:
                print(organism, " not in motif data")

            print("outputting to:", organism)

            if os.system(f'test -e {exact_motif_output}') != 0:  #if file doesn't already exist, write it 
                with open(exact_motif_output, "w") as to_write:
                    writer = csv.writer(to_write)
                    writer.writerow(["motif instance", "frequency"])
                    for motif_instance in motif_data[organism][0]:
                        motif_freq = motif_data[organism][0][motif_instance]
                        writer.writerow([motif_instance, motif_freq])
            
            if os.system(f'test -e {similar_motif_output}') != 0: #if file doesn't already exist, write it 
                with open(similar_motif_output, "w") as to_write:
                    writer = csv.writer(to_write)
                    writer.writerow(["motif instance", "frequency"])
                    for motif_instance in motif_data[organism][1]:
                        motif_freq = motif_data[organism][1][motif_instance]
                        writer.writerow([motif_instance, motif_freq])
            
    return 


if __name__ == "__main__":
    args = sys.argv
    print(args)
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])

