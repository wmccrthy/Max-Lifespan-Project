"""
DURING MEETING ON 4/16/24, DAVID SUGGESTED VALIDATING/CROSS-REFERENCING OUR DNA CLUSTERING ANAYLSIS WITH THE SAME ANALYSIS ON PROTEIN COUNTER PARTS
THAT IS, FOR EACH GENE DATASET:
    - TRANSLATE SEQUENCES TO PROTEIN 
    - RUN PROTEIN THRU SAME PIPELINE 
        - EMBED
        - PCA, K-MEANS CLUSTER
        - Z-TEST AND VISUALIZE 


THIS SCRIPT CONTAINS METHODS FOR PERFORMING THE PROCEDURE OUTLINED ABOVE; MAIN OBJECTIVES ARE
    - FOR EACH GENE DNA SET, CREATE CORRESPONDING GENE PROTEIN SET, CONSISTING OF EACH DNA SEQUENCE TRANSLATED TO PROTEIN 
    - FOR EACH GENE PROTEIN SET, COMPUTE AND SAVE EMBEDDINGS 
        - NEED TO FIND GOOD PRE-TRAINED PROTEIN EMBEDDING MODEL 
    - CLUSTER_ALL_EMBEDDINGS (THINK I CAN USE ALMOST IDENTICAL METHOD FROM GET_REGULATORY_EMBEDDINGS)
"""
import os, sys, csv
from Bio.Seq import Seq 

regulatory_dna_sets = "/data/rsg/chemistry/wmccrthy/Everything/gene_datasets/regulatory/"
regulatory_protein_sets = "/data/rsg/chemistry/wmccrthy/Everything/gene_datasets/regulatory_protein/"
#create alt directory regulatory_protein in which to store corresponding protein sequences 


"""
METHOD THAT ITERATES THROUGH GENE ORTHOLOG SETS, AND FOR EACH SET:  
    - FOR EACH LINE, L:
        - TRANSLATE DNA SEQUENCE TO PROTEIN
        - L_new = L with DNA replaced by PROTEIN
        - WRITE L_new to OUTPUT FILE GENE_PROTEINS
"""
def translate_gene_dna():
    num_translated = 0 
    for file in os.listdir(regulatory_dna_sets):
        if "embeddings" not in file: 
            gene = file.split("_")[0]
            file_path = regulatory_dna_sets + file
            protein_path = regulatory_protein_sets + f"{gene}_proteins.csv"
            print(f"translating {file} to proteins (#{num_translated})")
            with open(protein_path, "w") as write_to:
                writer = csv.writer(write_to)
                with open(file_path) as read_from:
                    for line in read_from:
                        line = line.split(",")
                        if len(line) < 9 or line[0] == 'organism': continue 
                        seq = line[-1].replace('"', "").replace("\n", "").replace("X", "").replace(" ", "")
                        # print(seq)
                        try:
                            protein_seq = Seq(seq).translate()
                        except: 
                            print("error translating seq:", seq, "\n proceeding to next seq")
                            continue 
                        # print(str(protein_seq), "\n")
                        writer.writerow(line[:-1] + [protein_seq])
            num_translated += 1
    return 



def get_protein_embeddings():


    return 




if __name__ == "__main__":
    args = sys.argv
    print(args)
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])
