import os, csv, sys, random

"""
SCRIPT FOR RETRIEVING AND OUTPUTTING TRAINING DATA CSV of (sequence, lifespan) format 

PSUEDO: 

    function takes in num_data parameter (indicating how much training data we want to compile)
    we want to iterate thru cumulative data set according to this parameter; that is, we want to select/use a data point every len(set)/num_data indices 
    interval = len(set)/num_data

    create species:lifespan mappings 
    
    iterate thru cumulative_data in intervals of 10,000 (starting from index i where i +=1 for each full iteration of data where we don't scribe (num_data) points)

"""
def get_data(num_data, lifespan_min=None, lifespan_max=None, unique_lifespans=0):
    num_data = int(num_data)

    lifespan_path = '/data/rsg/chemistry/wmccrthy/Everything/lifespan_data.csv'
    lifespan_mappings = {}
    with open(f'{lifespan_path}') as file:
        for line in file:
            line = line.strip().split(",")
            species_name = line[0]
            lifespan = line[1]
            if lifespan != "Unknown" and lifespan != "Not Established": lifespan_mappings[species_name] = lifespan

    cumulative_set_MIT_path = "/data/rsg/chemistry/wmccrthy/Everything/cumulativeTOGAset_trimmed.csv"

    output_path = f"/data/rsg/chemistry/wmccrthy/Everything/early_training/arbitrary_{num_data}_training.csv"
    if lifespan_min and lifespan_max: output_path = f"/data/rsg/chemistry/wmccrthy/Everything/early_training/lifespan_range_{lifespan_min}to{lifespan_max}_{num_data}_training.csv"
    if unique_lifespans != 0: output_path = f"/data/rsg/chemistry/wmccrthy/Everything/early_training/lifespan_range_{lifespan_min}to{lifespan_max}_unique_{num_data}_training.csv"
    valid_data = []
    lifespans_included = set()

    index = 0
    with open(cumulative_set_MIT_path) as read_from:
        for line in read_from:
            line = line.split(",")

            if len(line) < 9: continue  #avoid blank lines (just newlines for the most part)

            print(index)
            index += 1
            # print(index, next)
            # print(line, "\n")
            species, sequence = "_".join(line[0].split("_")[:2]), line[-1]
            sequence = sequence.replace("X", "").replace("\n", "").replace('"', "").replace(" ", "").upper()
            if species in lifespan_mappings: lifespan = lifespan_mappings[species]
            else: continue 

            if "-" not in sequence and len(sequence) >= 1024 and len(sequence) < 12000: #limit length of sequence for sake of testing model  
                lifespan = lifespan.replace("\n", "")
                if unique_lifespans != 0: 
                    if lifespan in lifespans_included: continue 
                    
                if lifespan_min == None or lifespan_max == None: 
                    valid_data.append([index, lifespan, sequence]) 
                    lifespans_included.add(lifespan)    
                else:
                    if float(lifespan) <= int(lifespan_max) and float(lifespan) >= int(lifespan_min): 
                        valid_data.append([sequence, lifespan])  
                        lifespans_included.add(lifespan)
    
    training_random = random.sample(valid_data, num_data)
    with open(output_path, "w") as write_to:
        writer = csv.writer(write_to)
        for d in training_random:
            writer.writerow(d)

     

def get_specific_data(num_data, gene_type):
    num_data = int(num_data)
    regulatory_set_path = f"/data/rsg/chemistry/wmccrthy/Everything/gene_datasets/regulatory/{gene_type}_orthologs_valid_trimmed.csv"
    to_sample = []
    with open(regulatory_set_path) as read_from:
        for line in read_from:
            #we just want sequence and lifespan
            line = line.split(",")
            if len(line) < 10 or line[0] == 'organism': continue 

            lifespan, seq = line[1], line[-1]
            if lifespan == 'Unknown': continue 
            if "-" in seq: continue #exclude gaps for now, for simplicity sake

            seq = seq.replace(" ", "").replace("\n", "").replace("X", "").replace('"', "")
            to_sample.append([seq, lifespan])

    num_data = min(num_data, len(to_sample))
    training = random.sample(to_sample, num_data)

    output_path = f"/data/rsg/chemistry/wmccrthy/Everything/early_training/training_data/{gene_type}_{num_data}_training.csv"
    with open(output_path, "w") as write_to:
        writer = csv.writer(write_to)
        for line in training:
            writer.writerow(line)

    return 

if __name__ == "__main__":
    args = sys.argv
    print(args)
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])




