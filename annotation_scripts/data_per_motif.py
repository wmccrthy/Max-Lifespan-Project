import os, csv, sys
from Bio import motifs
from Bio.Seq import Seq
from collections import Counter 

"""
with open(motif_output_path, "w") as write_to:
  writer = csv.writer(write_to)

    iterate thru all retrieved data (promoter regions)
      - for file in data_subsets:
          - for line in file... 
              - line = line.splt(",")
              - sequence = line[-1]
              - promoter_region = sequence[:1005]
              - if len(list(motifs.search(sequence))) > 0: writer.writerow(line)
"""

def get_kozak_orthologs(organism = None):
    canonical_kozak = "GCC(A/G)CCATGG"
    kozak_instances = [Seq("GCCACCATGG"), Seq("GCCGCCATGG")]
    get_orthologs_with_motif(kozak_instances, "kozak")

def get_tata_orthologs(organism = None):
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
    get_orthologs_with_motif(tata_instances, "tata")  

def get_caat_orthologs(organism = None):
    caat_instances = []
    canonical_caat = "GG(T/C)CAATCT"
    caat_instances = [Seq("GGCCAATCT"), Seq("GGTCAATCT")]
    get_orthologs_with_motif(caat_instances, "caat")

def get_orthologs_with_motif(motif_instances, motif_name):
    #create lifespan dict s.t we can add lifespan column to each row 
    lifespan_path = '/Mounts/rbg-storage1/users/wmccrthy/lifespan_data.csv'
    lifespan_mappings = {}
    with open(f'{lifespan_path}') as file:
        for line in file:
            line = line.strip().split(",")
            organism = line[0]
            lifespan = line[1]
            lifespan_mappings[organism] = lifespan

    motif = motifs.create([Seq(i) for i in motif_instances]) #create motif based on instances input (convert to Seq s.t input supports str or Seq type)    

    motif_output_path = f"/Mounts/rbg-storage1/users/wmccrthy/motif_datasets/orthologs_w_{motif_name}.csv"
    data_subsets_path = f"/Mounts/rbg-storage1/users/wmccrthy/dataSubSets/"
    
    with open(motif_output_path) as write_to:
        print(motif_output_path)
        writer = csv.writer(write_to)
        writer.writerow(['organism', 'max_lifespan', 'gene_id', 'orthologType', 'accession', 'chromosome', 'start', 'end', 'direction (+/-)', 'intactness', 'sequence'])
        for file in os.listdir(data_subsets_path):
            file_path = data_subsets_path + file
            with open(file_path) as read_from:
                for line in read_from:
                    line = line.split(",")
                    promoter_region = line[-1][:1005]
                    promoter_region = Seq(promoter_region)
                    if len(list(promoter_region.search(motif.alignment))) > 0: #if promoter region of sequence contains motif of interest, write line to output file 
                        species_name = "_".join(line[0].split("_")[:2])
                        species_lifespan = lifespan_mappings[species_name]
                        line = [line[0]] + [species_lifespan] + line[1:]
                        #incorporate lifespan to line -> write line to output 
                        writer.writerow(line)
        

"""
CONSIDER APPROACH THAT GETS ALL DATA FOR ALL MOTIFS AT ONCE (RATHER THAN 1 MOTIF AT A TIME, TAKING THE MOTIF OF INTEREST AS INPUT)

PSUEDOCODE: 
    - HAVE DICT OF MOTIFS WHERE MOTIF_DICT[MOTIF] = LIST OF PROMOTER REGION SEQUENCES THAT CONTAIN THIS MOTIF 
      ITERATE THRU ALL DATA SUBSETS:
        FOR EACH LINE:
            FOR EACH MOTIF:
                IF SEQUENCE.SEARCH(MOTIF.ALIGNMENT): MOTIF_DICT[MOTIF].ADD(LINE)

    FOR MOTIF IN MOTIF_DICT:
        OUTPUT FILE W ALL SEQUENCESI IN MOTIF_DICT[MOTIF]

*REMEMBER TO INCLUDE LIFESPAN DATA

*ONE BOTTLENECK FOR THIS APPROACH IS WE NEED A COMPREHENSIVE IDEA OF ALL THE MOTIFS WE WANT BEFOREHAND
THOUGH, WITH THAT BEING SAID, IF WE INITIALLY RUN THIS AND EVENTUALLY COME ACROSS MORE MOTIFS WE WANT TO LOOK FOR, WE CAN EMPLOY get_orthologs_with_motif (looks for 1 motif, passed as parameter)
"""

def get_all_motif_orthologs():
    #create lifespan dict s.t we can add lifespan column to each row 
    lifespan_path = '/Mounts/rbg-storage1/users/wmccrthy/lifespan_data.csv'
    lifespan_mappings = {}
    with open(f'{lifespan_path}') as file:
        for line in file:
            line = line.strip().split(",")
            organism = line[0]
            lifespan = line[1]
            lifespan_mappings[organism] = lifespan

    motif_set = get_all_motifs()

    motif_dict = {motif:set() for motif in motif_set}

    data_subsets_path = f"/Mounts/rbg-storage1/users/wmccrthy/dataSubSets/"
    
    index = 0
    #iterate thru all files that hold promoter region data 
    for file in os.listdir(data_subsets_path):
        if file[:3] != 'API':
            file_path = data_subsets_path + file
            with open(file_path) as read_from:
                for line in read_from:
                    print(index)
                    index += 1

                    line = line.split(",")
                    if line[0] == "organism" or len(line) != 10: continue 

                    promoter_region = line[-1][:1005]
                    promoter_region = Seq(promoter_region)

                    seq_to_eval = Seq(str(promoter_region).replace("-", "").replace(" ", "").replace("X", "").replace('"', "").strip().upper())
                    background = Counter(seq_to_eval)
                    background = {i:(background[i]/len(seq_to_eval)) for i in background}
                    print(background)

                    for motif in motif_set:
                        
                        # if len(list(promoter_region.search(motif.alignment))) > 0: 
                        #     species_name = "_".join(line[0].split("_")[:2])
                        #     species_lifespan = lifespan_mappings[species_name]
                        #     line = [line[0]] + [species_lifespan] + line[1:]
                        #     #incorporate lifespan to line -> write line to output 
                        #     motif_dict[motif].add(line)
                        # else:

                            #want to look on basis of motif-similarity if exact motif match not found (looking for just exact alignment will ignore interesting data)
                            #create pwm and pssm 
                            #create background distribution of nucleotides as such:
                            pwm = motif.counts.normalize(pseudocounts=background)
                            pssm = pwm.log_odds()

                            #search on basis of pssm (similarity "threshold" will matter here so explore means of finding informed threshold
                            #derive similarity threshold based on distribution of each nucleotide across entire promoter region and calculated threshold that will reap 
                            #false positive rate of 1% (this may need to be tweaked as we run scripts and learn more)

                            # distribution = pssm.distribution(background=background, precision=10**4)
                            # score_threshold = distribution.threshold_fpr(.05)
                            #waiving this step for now bc it's computationally very costly 

                            #if similar motif found w score threshold >= threshold required for a false positive rate of 1%:
                            #    add line to motif_dict[motif]
                            seq_similarity_search = pssm.search(seq_to_eval, threshold=8.0)

                            if len(list(seq_similarity_search)) > 0:
                                species_name = "_".join(line[0].split("_")[:2])
                                species_lifespan = lifespan_mappings[species_name]
                                line = [line[0]] + [species_lifespan] + line[1:]
                                #incorporate lifespan to line -> write line to output 
                                motif_dict[motif].add(tuple(line))

                            #for testing: print most similar sequence in similarity search 
                            continue 

                            most_similar = [0, None]

                            for pos, score in pssm.search(seq_to_eval, threshold=8.0):
                                # if pos > 0: print("Motif:", motif.name, "Consensus:", motif.consensus, "Similar Instance:", seq_to_eval[pos:pos+len(motif)], "Score: ", score)
                                #save position of motif w highest similarity 
                                if score > most_similar[0]: 
                                    if pos > 0: most_similar = [score, seq_to_eval[pos:pos+len(motif)]]
                                    else: most_similar = [score, seq_to_eval[pos:pos+len(motif)].reverse_complement()]
                                

                        
    data_output_path = "/Mounts/rbg-storage1/users/wmccrthy/motif_datasets/"
    for motif in motif_dict:
        motif_output_path = data_output_path + motif
        with open(motif_output_path, "w") as write_to:
            writer = csv.writer(write_to)
            writer.writerow(['organism', 'max_lifespan', 'gene_id', 'orthologType', 'accession', 'chromosome', 'start', 'end', 'direction (+/-)', 'intactness', 'sequence'])
            for line in motif_dict[motif]:
                writer.writerow(line)




def get_all_motifs():
    """
    use JASPAR data and Biopython JASPAR support to build and return comprehensive list of all motifs for transcription factors

    Psuedo:
        jaspar_data = open("jaspar_motifs.txt")
        for m in motifs.parse(jaspar_data, "jaspar"):        
    """

    jaspar_path = "/Mounts/rbg-storage1/users/wmccrthy/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
    motif_set= set()
    jaspar_data = open(jaspar_path)

    #for each motif in jaspar transcription factor binding sites data
    for m in motifs.parse(jaspar_data, "jaspar"):
        # print(m.name, m.matrix_id)
        # print(m.consensus, "\n")
        # print(m.alignment)
        #add motif to cumulative set
        motif_set.add(m)
    print(len(motif_set))

    #return cumulative set of motifs 
    return motif_set
   

if __name__ == "__main__":
    args = sys.argv
    print(args)
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])

