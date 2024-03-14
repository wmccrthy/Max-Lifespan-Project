"""
BROAD SCRIPT FOR THE PURPOSE OF TRIMMING DATA 

    - def trim_cumulative()
        - ITERATE THRU CUMULATIVE TOGA SET, WRITING A LINE TO OUTPUT IF AND ONLY IF ITS SEQUENCE IS NOT:
            - Lost
            - Uncertain Loss   
            - Missing
            - Paralogous Projection
            - containing gap '-' char for more than 20% of its total length 

    - def trim_gene_sets():
        - ITERATE THRU gene_datasets directory
            - for each line in file:
                - write to output if and only if it's sequence is NOT:
                    - Lost
                    - Uncertain Loss   
                    - Missing
                    - Paralogous Projection
                    - containing gap '-' char for more than 20% of its total length 
"""
import os, csv, sys
from collections import Counter

cumulative_untrimmed_path = "/Mounts/rbg-storage1/users/wmccrthy/cumulativeTOGAset.csv"

gene_datasets_path = "/Mounts/rbg-storage1/users/wmccrthy/gene_datasets/regulatory/"


invalid_intactness = set(["lost", "missing", "paralogous projection"])

def trim_cumulative():
    trimmed_path = "/".join(cumulative_untrimmed_path.split("/")[:-1]) + "/cumulativeTOGAset_trimmed.csv"
    print(trimmed_path)
    with open(trimmed_path, "w") as write_to:
        writer = csv.writer(write_to)
        writer.writerow(['organism','gene_id','orthologType','chromosome','start','end','direction (+/-)','intactness','sequence'])
        with open(cumulative_untrimmed_path) as read_from:
            for line in read_from:
                line = line.split(",")
                intactness = line[-2]
                if len(line) < 9: continue 
                if intactness in invalid_intactness: continue #if missing or lost, don't output 
                elif Counter(line[-1])['-'] > len(line[-1]) / 5: continue #if too many gaps in sequence, don't output 
                writer.writerow(line)

def trim_gene_sets():
    for file in os.listdir(gene_datasets_path):
        file_path = gene_datasets_path + file 
        trimmed_file_path = "/".join(file_path.split("/")[:-1]) + "/" + "".join(file.split(".")[:-1]) + "_trimmed.csv"
        # print(trimmed_file_path)
        with open(trimmed_file_path, "w") as write_to:
            writer = csv.writer(write_to)
            writer.writerow(['organism','max_lifespan', 'gene_id','orthologType','chromosome','start','end','direction (+/-)','intactness','sequence'])
            with open(file_path) as read_from:
                for line in read_from:
                    line = line.split(",")
                    if len(line) < 10: continue 
                    intactness = line[-2]
                    if intactness in invalid_intactness: continue #if missing or lost, don't output 
                    elif Counter(line[-1])['-'] > len(line[-1]) / 5: continue #if too many gaps in sequence, don't output 
        
                    writer.writerow(line)
        write_to.close()

def ensure_gene():
    gene_types = set()
    for file in os.listdir(gene_datasets_path):
        type = file.split("_")[0]
        gene_types.add(type)
    print(len(gene_types))
    return gene_types

def remove_trimmings():
    num_rem = 0
    for file in os.listdir(gene_datasets_path):
        file_path = gene_datasets_path + file
        if file[-11:] == 'trimmed.csv': 
            num_rem += 1
            os.system(f"rm -rf {file_path}")
    print(num_rem)
    
if __name__ == "__main__":
    args = sys.argv
    print(args)
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])

