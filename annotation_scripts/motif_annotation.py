from Bio import motifs
from Bio.Seq import Seq 
import csv, sys, os 


"""
WRITE A METHOD annotate_motif(), that given instances of a motif and a species:
    - motif = motifs.create(instances)
    - pwm = motif.counts.normalize(pseudocounts=0.5)
    - pssm = pwm.log_odds()

    - finds that species within our data sub sets 
        - for file in path/dataSubSets
            - iterate thru file until finding species (if don't find species, will simply proceed to next file)

    - iterates through species genes 
        - for each: 
            - if motif is found in current sequence: 
                - record data on position and instance of motif found 
            - also record positional similarity "scores" using pssm 
                - allows us to consider more flexible range of variants of given motif 

    CONCERNS
    - we might want to search different ranges of the sequence depending on the motif in question
    - we might care about different data depending on the motif in question 
        - BIG QUESTION: what data do we actually care about and want to record? 
"""

def get_start_codon(species, sample_size):
    start_codon_instances = [Seq("ATG")]
    annotate_motif(start_codon_instances, species, sample_size)
#honestly just testing to make sure start codon is at pos 0 for every gene 

def get_kozak_annotations(species, sample_size):
    canonical_kozak = "GCCRCCATGG"
    kozak_instances = [Seq("GCCACCATGG"), Seq("GCCGCCATGG")]
    annotate_motif(kozak_instances, species, sample_size)

def get_tata_annotations(species, sample_size):
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
    annotate_motif(tata_instances, species, sample_size)
    

data_subset_path = f"/Mounts/rbg-storage1/users/wmccrthy/dataSubSets/"

def annotate_motif(instances, species, sample_size):
    seen = set() #set used to ensure we only consider TATA occurences from same gene once (occasional duplicate in the data for some reason)

    motif = motifs.create([Seq(i) for i in instances]) #create motif based on instances input (convert to Seq s.t input supports str or Seq type)
    pwm = motif.counts.normalize(pseudocounts=0.5)
    pssm = pwm.log_odds()
    index = 0
    sample_size = int(sample_size)
    print(sample_size)
    motif_instance_frequency = {}
    similar_instance_frequency = {}

    for file in os.listdir(data_subset_path):
        print("checking in", file, "\n")
        with open(f'{data_subset_path}{file}') as cur_subset:
            for line in cur_subset:
                if index >= sample_size: break #exit when reached desired sample size

                line = line.split(",")
                cur_species = line[0]

                if cur_species != species: continue #proceed until we find inputted species 
                
                #avoid duplicates 
                chrom, start, end  = line[4:7]
                if (chrom, start, end) in seen: continue 
                else: seen.add((chrom, start, end)) 
                
                sequence, orientation = line[-1], line[-3]
                
                #format sequence (get rid of spaces, and "" marks)
                sequence = sequence.replace('"', '').replace(" ", "").strip().upper()

                
                #section off relevant portion of sequence (should always be first 1050 base pairs to accomdate for region before and slightly after start codon)
                sequence = sequence[:1005]
                print(sequence[min(998, len(sequence)-2):])

                #preliminary fix for duplicate 'A' produced by querying 1000 BPs upstream from start rather than from start-1
                if sequence[-7:-2] == 'AAATG': 
                    sequence = sequence[:-7] + sequence[-6:]
                    print("reduced to: ", sequence[min(998, len(sequence)-2):])

                sequence = Seq(sequence)
                start_codon = 1000 #start codon generally assumed to be at position 1000; perhaps add functionality s.t if seq at 1000 is not ATG, we make note s.t data is less 'credible'
                

                closest_motif = [None, None, float('inf')] #variable for storing instance of motif found closest to start codon
                if len(list(sequence.search(motif.alignment))) > 0:
                    index += 1 #used for containing method to inputted sample size 

                    for pos, seq in sequence.search(motif.alignment): 
                        distance_from_start_codon = start_codon - pos
                        if distance_from_start_codon < closest_motif[2]: 
                            closest_motif = [seq, pos, distance_from_start_codon]
                    
                    #record frequencies of instance of motic found closest to start codon 
                    if str(closest_motif[0]) not in motif_instance_frequency: motif_instance_frequency[str(closest_motif[0])] = 1
                    else: motif_instance_frequency[str(closest_motif[0])] += 1

                #search for similar motif instances even if exact match isn't found 
                #little sus about capturing the right data when pos is negative (indicating it pertains to sequence on the reverse strand)
                highest_scored_pos, most_similar_seq = 0, None 
                for pos, score in pssm.search(sequence, threshold= 7.5):
                    similar_seq = sequence[pos:pos+motif.length]
                    # if pos < 0: similar_seq = sequence[pos:pos+motif.length].reverse_complement()
                    distance_from_start_codon = start_codon - abs(pos)
                        
                    #look for highest similarity-scoring pos/sequence that is closer than exact motif matches found 
                    if score > highest_scored_pos and distance_from_start_codon < closest_motif[2]: 
                        most_similar_seq = [similar_seq, pos, distance_from_start_codon]
                        highest_scored_pos = score
                        # print(most_similar_seq)

                #record frequencies of sequences most similar to motif (if closer than actual found motif)
                if most_similar_seq: 
                    if str(most_similar_seq[0]) not in similar_instance_frequency: similar_instance_frequency[str(most_similar_seq[0])] = 1
                    else: similar_instance_frequency[str(most_similar_seq[0])] += 1

                        # print(similar_seq, pos, distance_from_start_codon, score)
                    
                print("Closest Motif found: ", str(closest_motif[0]), "| at pos", str(closest_motif[1]), "|", str(closest_motif[2]), "base pairs from start codon", orientation)
                if most_similar_seq: print("Most Similar Sequence found: ", str(most_similar_seq[0]), "| at pos", str(most_similar_seq[1]), "|", str(most_similar_seq[2]), "base pairs from start codon")
                print()

    print(sorted([(i, motif_instance_frequency[i]) for i in motif_instance_frequency], key = lambda x:x[1]))
    print(sorted([(i, similar_instance_frequency[i]) for i in similar_instance_frequency], key = lambda x:x[1]))
    return 

if __name__ == "__main__":
    args = sys.argv
    print(args)
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])
