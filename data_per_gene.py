import os, csv, sys 


def get_data_per_gene():
    #create lifespan dict s.t we can add lifespan column to each row 
    lifespan_path = '/Mounts/rbg-storage1/users/wmccrthy/lifespan_data.csv'
    lifespan_mappings = {}
    with open(f'{lifespan_path}') as file:
        for line in file:
            line = line.strip().split(",")
            organism = line[0]
            lifespan = line[1]
            lifespan_mappings[organism] = lifespan

    cumulative_set_path = f"/Mounts/rbg-storage1/users/wmccrthy/cumulativeTOGAset.csv"
    output_path = f"/Mounts/rbg-storage1/users/wmccrthy/gene_datasets/"

    genes = {}
    index = 0   

    with open(cumulative_set_path) as read_from:
        for line in read_from:
            line = line.split(",")
            if line[0] == "organism": continue #skip header

            # gene_id = ".".join(line[1].split(".")[:2])
            gene_symbol = line[1].split(".")[1]
            # print(gene_id)
            if gene_symbol not in genes: genes[gene_symbol] = [line]
            else: genes[gene_symbol].append(line)
            index += 1
            print(index, "\n")

    print(len(genes.keys()))
    print(genes.keys())

    for gene in genes:
        with open(f"{output_path}{gene}_orthologs.csv", "w") as write_to:
            writer = csv.writer(write_to)
            writer.writerow(['organism', 'max_lifespan', 'gene_id', 'orthologType', 'accession', 'chromosome', 'start', 'end', 'direction (+/-)', 'intactness', 'sequence'])     
            for line in genes[gene]:

                if line[0] == "organism": continue #skip header
                if len(line) < 9: 
                    continue #skip ed formatted data 

                species_name = "_".join(line[0].split("_")[:2])
                species_lifespan = lifespan_mappings[species_name]
                line = [line[0]] + [species_lifespan] + line[1:]
                #incorporate lifespan to line -> write line to output 
                writer.writerow(line)



#categorizes gene data sets into regulatory or non-regulatory (each gene data sets contains all TOGA orthologs pertaining to that gene)
def organize_data_per_gene():
    gene_data_path =  f"/Mounts/rbg-storage1/users/wmccrthy/gene_datasets/"
    transciption_factor_path = gene_data_path + "transcription_factors.txt"
    tf_set = set()
    with open(transciption_factor_path) as read_from:
        for line in read_from:
            line = line.split("\t")
            gene_symbol = line[0]
            tf_set.add(gene_symbol)
    print(tf_set)
    
    regulatory_path = gene_data_path + "regulatory/"
    non_regulatory_path = gene_data_path+ "non_regulatory/"

    for gene_file in os.listdir(gene_data_path):
        if gene_file[-3:] != "transcription_factors.txt" and gene_file != "regulatory" and gene_file != "non_regulatory":
            gene_symbol = gene_file.split("_")[0]
            print(gene_symbol)

            gene_file_path = gene_data_path + gene_file

            if gene_symbol in tf_set: os.system(f"mv {gene_file_path} {regulatory_path}")
            else: os.system(f"mv {gene_file_path} {non_regulatory_path}")

"""

#iterate thru all gene_ortholog files and add lifespan column to each entry 
def add_lifespan():
    #create lifespan dict s.t we can add lifespan column to each row 
    lifespan_path = '/Mounts/rbg-storage1/users/wmccrthy/lifespan_data.csv'
    lifespan_mappings = {}
    with open(f'{lifespan_path}') as file:
        for line in file:
            line = line.strip().split(",")
            organism = line[0]
            lifespan = line[1]
            lifespan_mappings[organism] = lifespan
        
    gene_data_path =  f"/Mounts/rbg-storage1/users/wmccrthy/gene_datasets/"

    #iterate thru every gene file, 
    #   write new output to file_name.csv 
    #   for each line in existing file 
    #       add lifespan data and write to new output
    #   once done iterating thru all lines remove prior file 
    file_num = 0
    for dir in os.listdir(gene_data_path):
        dir_path = gene_data_path + dir + "/"
        for file in os.listdir(dir_path):
            file_path = dir_path + file
            print(file_num, file)
            file_num += 1
            with open(f"{file_path}.csv", "w") as write_to: 
                writer = csv.writer(write_to)
                writer.writerow(['organism', 'max_lifespan', 'gene_id', 'orthologType', 'accession', 'chromosome', 'start', 'end', 'direction (+/-)', 'intactness', 'sequence'])
                with open(file_path) as read_from:
                    for line in read_from:
                        line = line.split(",")
                        if line[0] == "organism": continue #skip header
                        if len(line) < 9: 
                            continue #skip ed formatted data 
                        species_name = "_".join(line[0].split("_")[:2])
                        species_lifespan = lifespan_mappings[species_name]
                        line = [line[0]] + [species_lifespan] + line[1:]
                        # print(line)
                        # writer.writerow(line)
            os.system(f"rm -rf {file_path}")
"""


if __name__ == "__main__":
    args = sys.argv
    print(args)
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])



