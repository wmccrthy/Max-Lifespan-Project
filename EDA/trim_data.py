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

    - def create_one2one_set(file_path):
        - iterate thru given toga file path, writing a line to output if and only if it's sequence is:
            - one2one
            - between lengths of 104 and 31620

"""
import os, csv, sys, numpy as np
import matplotlib.pyplot as plt
from collections import Counter

cumulative_untrimmed_path = "/data/rsg/chemistry/wmccrthy/Everything/cumulativeTOGAset.csv"

gene_datasets_path = "/data/rbg/users/wmccrthy/chemistry/Everything/gene_datasets/regulatory/"
gene_one2one_datasets_path = "/data/rbg/users/wmccrthy/chemistry/Everything/gene_datasets/regulatory_one2one/"


invalid_intactness = set(["lost", "missing", "paralogous projection"])

"""
method to trim a given csv file path (expected format: ['organism','gene_id','orthologType','chromosome','start','end','direction (+/-)','intactness','sequence'])
generates new file only including one2one sequences btwn lengths of 104 and 31620 
"""
def create_one2one_set(file_path):
    new_path = file_path.split(".")[0] + "_one2one.csv"
    print(new_path)
    total_data = 0
    num_data = 0
    with open(new_path, "w") as write_to:
        writer = csv.writer(write_to)
        writer.writerow(['organism','gene_id','orthologType','chromosome','start','end','direction (+/-)','intactness','sequence'])
        with open(file_path) as read_from:
            for line in read_from:
                line = line.split(",")
                if len(line) < 9: continue
                total_data += 1
                seq = line[-1].strip()
                length = len(seq)
                if length < 104 or length > 31620: continue #only include seqs of length btwn this range
                ortholog_type = line[2]
                if ortholog_type != "one2one": continue #only include one2one seqs 
                num_data += 1
                writer.writerow(line)
    print("Total Data: ", total_data)
    print("Valid Data: ", num_data)


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

"""
method to trim a given gene's csv file path (expected format: ['organism','max_lifespan', 'gene_id','orthologType','chromosome','start','end','direction (+/-)','intactness','sequence'])
generates new file only including one2one sequences btwn lengths of 104 and 31620 (put in new directory for organization's sake)
also computes stats of # of sequences and # of species represented per gene (we can create histogram of these...)
"""
def create_one2one_gene_sets():
    gene_stats = {} #dict to hold gene:(num_data_entries, num_species) pairs
    num_genes = 0
    for file in os.listdir(gene_datasets_path):
        file_path = gene_datasets_path + file 
        new_file_path = "/".join(file_path.split("/")[:-2]) + "/regulatory_one2one/" + "".join(file.split(".")[:-1]) + "_one2one.csv"
        #only execute on most up-to-date 'trimmed' gene files 
        if file_path[-21:] != "orthologs_trimmed.csv": continue
        num_genes += 1
        cur_gene = file.split("_")[0]
        num_data = 0
        organisms = set()

        # print(new_file_path, cur_gene)

        # continue

        with open(new_file_path, "w") as write_to:
            writer = csv.writer(write_to)
            writer.writerow(['organism','max_lifespan', 'gene_id','orthologType','chromosome','start','end','direction (+/-)','intactness','sequence'])
            with open(file_path) as read_from:
                for line in read_from:
                    line = line.split(",")
                    if len(line) < 10: continue 
                    seq = line[-1].strip()
                    length = len(seq)
                    if length < 104 or length > 31620: continue #only include seqs of length in this range
                    ortholog_type = line[3]
                    if ortholog_type != "one2one": continue #only include one2one seqs 
                    num_data += 1
                    print(line[0])
                    organisms.add(line[0])
                    writer.writerow(line)

        gene_stats[cur_gene] = (num_data, len(organisms)) #fill dict entry for current gene

        write_to.close()
        os.system(f'rm -rf {file_path}')

    # after iterating thru all genes s.t we've generated one2one, trimmed, gene sets
    # create new file that outputs stats per gene (csv of format [gene | num_entries | num_species])
    # we can use this file to generate histograms thereafter
    
    with open("/data/rbg/users/wmccrthy/chemistry/Everything/EDA/regulatory_one2one_sets_metadata.csv", "w") as write_to:
        writer = csv.writer(write_to)
        writer.writerow(["gene", "# seqs", "# species"])
        for gene in gene_stats:
            num_seqs, num_species = gene_stats[gene]
            writer.writerow([gene, num_seqs, num_species])
    write_to.close()


"""
Method that iterates through one2one gene datasets and outputs csv file with [gene, num sequences, num species]
"""
def get_one2one_stats():
    gene_stats = {} #dict to hold gene:(num_data_entries, num_species) pairs
    # after iterating thru all genes s.t we've generated one2one, trimmed, gene sets
    # create new file that outputs stats per gene (csv of format [gene | num_entries | num_species])
    # we can use this file to generate histograms thereafter
    for file in os.listdir(gene_one2one_datasets_path):
        file_path = gene_one2one_datasets_path + file
        cur_gene = file.split("_")[0]
        organisms = set()
        num_data = 0
        with open(file_path) as read_from:
            for line in read_from:
                line = line.split(",")
                if len(line) < 10: continue #don't read empty lines
                org = "_".join(line[0].split("_")[:2])
                # print(org)
                organisms.add(org)
                num_data += 1
        gene_stats[cur_gene] = (num_data, len(organisms))
        # print(cur_gene, num_data, len(organisms))
        if len(organisms) > 450:
            print("============= WARNING =============")
            print(len(organisms), "> total species")

    with open("/data/rbg/users/wmccrthy/chemistry/Everything/EDA/regulatory_one2one_sets_metadata.csv", "w") as write_to:
        writer = csv.writer(write_to)
        writer.writerow(["gene", "# seqs", "# species"])
        for gene in gene_stats:
            num_seqs, num_species = gene_stats[gene]
            writer.writerow([gene, num_seqs, num_species])
    write_to.close()


"""
Iterate through regulatory_one2one_sets_metadata.csv and create histograms for each:
    - distribution of # sequences per gene
    - distribution of # species per gene
"""
def gene_one2one_histograms():
    data_path = "/data/rbg/users/wmccrthy/chemistry/Everything/EDA/regulatory_one2one_sets_metadata.csv"
    num_seqs = []
    num_species = []
    with open(data_path) as read_from:
        for line in read_from:
            line = line.split(",")
            if len(line) < 3: continue #avoid janky csv parsing stuff
            seqs, species = int(line[1].strip()), int(line[2].strip())
            num_seqs.append(seqs)
            num_species.append(species)
    
    num_seqs = np.array(num_seqs)
    num_species = np.array(num_species)
    # Generate the histogram
    plt.figure(figsize=(10, 6))
    plt.hist(num_seqs, bins=25, edgecolor='black')  # Adjust bins for better resolution
    plt.title('Histogram of # Seqs per Gene')
    plt.xlabel('# Seqs')
    plt.ylabel('# Genes')
    plt.grid(True)
    plt.show()
    plt.hist(num_species, bins=25, edgecolor='black')  # Adjust bins for better resolution
    plt.title('Histogram of # Species per Gene')
    plt.xlabel('# Species')
    plt.ylabel('# Genes')
    plt.grid(True)
    plt.show()



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

