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

cumulative_untrimmed_path = "/data/rsg/chemistry/wmccrthy/Everything/cumulativeTOGAset.csv"

gene_datasets_path = "/data/rsg/chemistry/wmccrthy/Everything/gene_datasets/regulatory/"


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

        #don't execute on newly created 'trimmed' files 
        if 'trimmed' in file: continue

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
        os.system(f'rm -rf {file_path}')

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



"""
ADDS LIFESPAN DATA FOR 39 SPECIES WHO WERE PREVIOUSLY MISSING LIFESPANS TO TRIMMED REGULATORY GENE SETS 
one-time use 
"""
def update_missing_lifespans():
    #create lifespan dict s.t we can add lifespan column to each row 
    lifespan_path = '/data/rsg/chemistry/wmccrthy/Everything/lifespan_data.csv'
    lifespan_mappings = {}
    with open(f'{lifespan_path}') as file:
        for line in file:
            line = line.strip().split(",")
            organism = line[0]
            lifespan = line[1]
            lifespan_mappings[organism] = lifespan

    updated_species = set() #for testing purposes | to ensure we are updating expected lifespans 

    #iterate thru all un-updated sets
    for file in os.listdir(gene_datasets_path):

        if "updated" in file: continue #to ensure we don't infinitely loop as files are written to this directory throughout loop 

        file_path = gene_datasets_path + file 
        trimmed_file_path = "/".join(file_path.split("/")[:-1]) + "/" + "".join(file.split(".")[:-1]) + "_updated_trimmed.csv"
        # print(trimmed_file_path)
        with open(trimmed_file_path, "w") as write_to:
            writer = csv.writer(write_to)
            writer.writerow(['organism','max_lifespan', 'gene_id','orthologType','chromosome','start','end','direction (+/-)','intactness','sequence'])
            with open(file_path) as read_from:
                for line in read_from:
                    line = line.split(",")
                    if len(line) < 10: continue 
                    if line[0] == 'organism': continue 
                    species_name = "_".join(line[0].split("_")[:2])
                    species_lifespan = lifespan_mappings[species_name]
                    if line[1] == 'Unknown' or line[1] == 'Not Established': 
                        line[1] = species_lifespan 
                        if species_lifespan != 'Unknown': updated_species.add(species_name)
                    writer.writerow(line)
            os.system(f"rm -rf {file_path}")

    print(len(updated_species), updated_species)

"""
renames regulatory sets w updated lifespan data from "xxx_updated_trimmed.csv" to "xxx_trimmed.csv" for compatability w prior scripts 
"""
def rename_updated():
    for file in os.listdir(gene_datasets_path):
        file_path = gene_datasets_path + file
        if "trimmed_updated_trimmed" in file or "valid_trimmed" in file: #only files we care abt 
            new_name = "_".join(file.split("_")[:2])
            new_name += "_valid_trimmed.csv"
            print("renaming", file, "to", new_name, "\n")
            new_path = gene_datasets_path + new_name
            os.system(f"mv {file_path} {new_path}")
        else: 
            #get rid of redundant/old files 
            os.system(f"rm -rf {file_path}")
            print(f"removing {file}")
    

def remove_scrap_embeddings():
    for file in os.listdir(gene_datasets_path):
        file_path = gene_datasets_path + file
        if "embeddings" in file:
            print("removing:", file)
            os.system(f'rm -rf {file_path}')
        else: continue 

if __name__ == "__main__":
    args = sys.argv
    print(args)
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])

