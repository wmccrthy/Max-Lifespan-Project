import gzip
import os, sys, genome_browsing, time, csv
from Bio import motifs
from Bio.Seq import Seq 

hg_38_sequences = []

# The base path(s) upon which to start parsing TOGA files | change depending on if testing locally or running on servers 
MIT_human_ref_path = "/Mounts/rbg-storage1/users/wmccrthy/human_hg38_reference" 
human_reference_path = "/Users/wyattmccarthy/Desktop/MIT Aging Project/Pre-Processing/human_hg38_reference"

#  CREATE TATA MOTIF 
canonical_tata = "TATAWAWN"
def get_TATA_variants(variant):
    if len(variant) >= len(canonical_tata): 
        instances.append(variant)
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


instances = []
get_TATA_variants("TATA")
instances = [Seq(tata_variant) for tata_variant in instances]
tata_motif = motifs.create(instances)
pwm = tata_motif.counts.normalize(pseudocounts=0.5)
pssm = pwm.log_odds()
print(pssm)

print(tata_motif.alignment)

ind = 0

with open("/Users/wyattmccarthy/Desktop/MIT Aging Project/Pre-Processing/tata_logo_preliminary/hg_38_positive_tata.fa", "w") as write_to:
    for mammalian_order_dir in os.listdir(human_reference_path):
        if mammalian_order_dir == ".DS_Store": continue #bug fix | os.listdir always includes this file, just skip over

        print(f'{mammalian_order_dir}:') #output the current order being parsed for clarity 

        # for each organism/species in order directory 
        for organism_dir in os.listdir(human_reference_path + f'/{mammalian_order_dir}'):
            if organism_dir == ".DS_Store": continue #bug fix | os.listdir always includes this file, just skip over

            interior_dir_path = human_reference_path + f'/{mammalian_order_dir}/{organism_dir}' #interior_dir_path stores the path to this species (ex. human_hg38_reference/order/species/)
            
            print(interior_dir_path)

            #parse codonAlignments file (get actual sequences)
            with gzip.open(f'{interior_dir_path}/codonAlignments.fa.gz') as file:

                prevIdentifier = None #fasta format uses '>' to denote identifying lines | prevIdentifer will be used to check if found sequence pertains to reference or query genome 
                unAccounted = 0

                for line in file:
                    # ord() gets ascii value of passed char; bc gzip reads in bytes we need to check for symbolic chars like such 
                    if line[0] == ord('>'):
                        line = line.split(b'|') # formatting line s.t there are no blank spaces | line now of format [geneID, codon, query/reference]
                        for i in range(len(line)): line[i] = line[i].strip() 

                        line[0] = line[0][1:] #remove '>' from geneID of identifier line 
                        prevIdentifier = line #set prevIdentifier s.t we can determine if sequences found thereafter pertain to reference or query 
                    else: 
                        # line containing DNA sequence 
                        # we only care abt sequences corresponding to QUERY genome 
                        if prevIdentifier[2] == b'REFERENCE':
                            
                            sequence = line.strip().decode()
                            print("checking", sequence)
                            
                            #find TATA box in sequence 
                            #format sequence (get rid of spaces, and "" marks)
                            sequence = sequence.replace('"', '').replace(" ", "").strip().upper()

                            sequence = Seq(sequence)

                            identified_tata = sequence.search(tata_motif.alignment)
                
                            if len(list(identified_tata)) > 0:
                                for pos, seq in sequence.search(tata_motif.alignment):
                                    ind += 1
                                    write_to.write(">" + "hg38" + "\n" + str(seq) + "\n")
                                    if ind >= 1000: break 
                                    # print("%i %s" % (pos, seq), " found ", distance_from_tss, " base pairs from TSS")
                                    # print(sequence[pos-15:pos+25])

        break 
        print()
