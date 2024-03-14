from Bio import motifs
from Bio.Seq import Seq 
import csv, sys 

"""
PSUEDO CODE FOR ANNOTATING TATA BOX IN MINED-DATA

DEF FIND_TATA(SEQUENCE, ORIENTATION): 
    GIVEN A REFERENCE SEQUENCE, LOOKS FOR OCCURANCE OF CANONICAL TATA BOX WITHIN REFERENCE AND RETURNS LOCATION RELATIVE TO TSS 
        - TSS is computed as index 1000 if orientation is '+'
        - TSS computed as index (len(sequence) - 1000) if orientation is '-'

    SEE BIOPYTHON DOCUMENTATION: https://biopython-tutorial.readthedocs.io/en/latest/notebooks/14%20-%20Sequence%20motif%20analysis%20using%20Bio.motifs.html#Position-W[…]t-Matrices
        - command 'f': "searching for instances"

instances = [all variations of canonical tata box]
    - TATAWAWN, where W can be an A/T and N can be any nucleotide base.

species_sequences = [all 1000 bp sequences from species we want]
m = motifs.create(instances)

for seq in species_sequences:
    for pos, seq in test_seq.search(m.alignment):
        print("%i %s" % (pos, seq))

"""

#  CREATE TATA MOTIF 
canonical_tata = "TATAWAWN"
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

tata_instances = []
get_TATA_variants("TATA")
tata_instances = [Seq(tata_variant) for tata_variant in tata_instances]
tata_motif = motifs.create(tata_instances)
pwm = tata_motif.counts.normalize(pseudocounts=0.5)
pssm = pwm.log_odds()
print(pssm)

print(tata_motif.alignment)

#  MOCK FOR EXTRACTING RELEVANT PORTION OF SEQUENCE GIVEN HOW WE'VE FORMATTED DATA 

# data_subset_path = f"/Mounts/rbg-storage1/users/wmccrthy/dataSubSets/temp_API_subSet"

def generate_TATAfasta_for_weblogo(species, path):
    ind = 0
    logo_sequence_path = f"/Mounts/rbg-storage1/users/wmccrthy/logoTests/positive_tata_sequences_{species}.fa"

    tata_positions = {} #dict for holding average score for positions at which tata is found within species' sequences 
    seen = set()

    with open(logo_sequence_path, "w") as write_to:
        with open(path) as file:
            for line in file:

                if ind >= 1000: break #testing to see what logo is generated given 100 tata-containing sequences 

                line = line.split(",")
                chrom, strt, end = line[4:7]

                if (chrom, strt, end) not in seen: #ensures we don't consider duplicates 
                    seen.add((chrom, strt, end))
                    
                    if len(line) < 9 or line[0] != species: continue
                    #only want data pertaining to given species 

                    sequence, orientation, start_codon = line[-1], line[-3], None 

                    #format sequence (get rid of spaces, and "" marks)
                    sequence = sequence.replace('"', '').replace(" ", "").strip().upper()

                
                    #section off relevant portion of sequence according to orientation
                    if orientation == '+': 
                        sequence = sequence[:1000]
                        start_codon = min(1000, len(sequence))
                    else: 
                        sequence = sequence[len(sequence)-1000:]
                        start_codon = max(len(sequence) - 1000, 0)

                    sequence = Seq(sequence)

                    if len(sequence) != 1000: continue #only get sequences of length 1000, s.t we can align them and consider via logo chart (needs sequences of same length)

                    identified_tata = sequence.search(tata_motif.alignment)
        
                    if len(list(identified_tata)) > 0:
                        print(line[:-1])
                        ind += 1
                        closest_TATA = ["", float('inf')] #variable for storing tata box found closest to TSS 

                        for pos, seq in sequence.search(tata_motif.alignment):
                            distance_from_start_codon = None 
                            if orientation == '+': 
                                distance_from_start_codon = start_codon - pos
                                if distance_from_start_codon < closest_TATA[1]: closest_TATA = [seq, distance_from_start_codon]
                            else: 
                                distance_from_start_codon = pos   
                                if distance_from_start_codon < closest_TATA[1]: closest_TATA = [seq, distance_from_start_codon]
                            # print("%i %s" % (pos, seq), " found ", distance_from_tss, " base pairs from TSS")
                            # print(sequence[pos-15:pos+25])


                        """
                        THIS SECTION CONSIDERS "SIMILARITY SCORES", SO TATA-LIKE ELEMENTS THAT ARE FAR MISMATCHED 
                        """
                        for pos, score in pssm.search(sequence, threshold = 5.0):
                            similar_seq = sequence[pos:pos+8]
                            
                            if pos < 0: similar_seq = sequence[pos:pos+8].reverse_complement()

                            # print("Position %d: score = %5.3f" % (pos, score))
                            if pos not in tata_positions: tata_positions[pos] = [score, 1]
                            else: 
                                tata_positions[pos][0] += score
                                tata_positions[pos][1] += 1
                            if orientation == '+': 
                                distance_from_start_codon = start_codon - abs(pos)
                                if distance_from_start_codon < closest_TATA[1]: closest_TATA = [similar_seq, distance_from_start_codon]
                            else: 
                                distance_from_start_codon = abs(pos)   
                                if distance_from_start_codon < closest_TATA[1]: closest_TATA = [similar_seq, distance_from_start_codon]
                    
                        
                        print("Closest TATA Found: ", str(closest_TATA[0]), " ", closest_TATA[1], " base pairs from TSS")
                        write_to.write(">" + str(line[4:7] + [closest_TATA[1]]) + "\n" + str(closest_TATA[0]) + "\n")
                        #the line above writes out the TATA box in current sequence that is closest to the TSS 

                        # write_to.write(">" + str(line[4:7]) + "\n" + str(sequence) + "\n")
                        # the line above writes out sequence in which tata box was found 

                        print()

        print(sorted([(tata_positions[i][0]/tata_positions[i][1], i) for i in tata_positions], key=lambda x:x[0]))
        #the above line prints out the positions at which TATA was found, sorted by their 'score' (a measure of TATA frequency occurence at that position)
                    

if __name__ == "__main__":
    args = sys.argv
    print(args)
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])




"""
ARBITRARY SPECIES 


- Vicugna_pacos__alpaca__HLvicPac4 /Mounts/rbg-storage1/users/wmccrthy/dataSubSets/2bitSubSet_7804645_11706962.csv
- Myotis_lucifugus__little_brown_bat__HLmyoLuc1 /Mounts/rbg-storage1/users/wmccrthy/dataSubSets/2bitSubSet_15609283_19511604.csv

- Leptonychotes_weddellii__Weddell_seal__HLlepWed2 /Mounts/rbg-storage1/users/wmccrthy/dataSubSets/2bitSubSet_3902322_7804644.csv
- Hydrochoerus_hydrochaeris__capybara__HLhydHyd2 /Mounts/rbg-storage1/users/wmccrthy/dataSubSets/2bitSubSet_11706963_15609282.csv
- Elephas_maximus__Asiatic_elephant__HLeleMax1 /Mounts/rbg-storage1/users/wmccrthy/dataSubSets/2bitSubSet_
    - need to try all for him bc unsure 

command: python3 TATA_annotation.py generate_fasta_for_weblogo {species} {path}

"""

