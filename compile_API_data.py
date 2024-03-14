import os, csv 


"""
BECAUSE OF HOW AIDS THE RETRIEVAL HAS BEEN FOR API DATA, THERE ARE BOUND TO BE DUPLICATES 

THIS SCRIPT IS FOR COMPILING ALL COLLECTED API DATA IN A SET (S.T DUPLICATE PROBLEM IS RESOLVED)
AND WRITING THAT SET TO A CUMULATIVE_API_SET.CSV FILE 

PSUEDO: 
initialize set 
iterate thru all temp_API_subsets 
    for line in file:
        set.add(line[:-1])

with output open:
iterate thru all temp_API_subsets again:
    for line in file
        if line in set: write data to output cumulative_API_subset.csv -> remove line 
        else: continue 

doing it this way will be less memory intensive bc we don't save sequence to set | in the output stage, once we see a data point, we write it and remove it from set s.t if seen again it won't be outputted 

"""

total_data = 0
#variables for comparing total data points w len(set) to gauge how many duplicates were present 

data_path =  f"/Mounts/rbg-storage1/users/wmccrthy/dataSubSets/" 
cumulative_set = set()

#build set of all unique data points (identifed by all features except sequence to conserve memory)
for file in os.listdir(data_path):
    if file[:4] == "temp":
        file_path = data_path + file
        with open(file_path) as read_from:
            for line in read_from:
                line = line.split(",")
                if len(line) != 10: continue 
                cumulative_set.add(tuple(line[:-1]))
                total_data += 1

#debugging/stat output
print("total data read:", total_data)
print("unique data read:", len(cumulative_set))

total_data = 0
uniq_data = 0

#re-iterate thru API dataset files, for each line we haven't already outputted to cumulative API dataset -> write output and remove from duplicate checking set s.t we don't write same line again
with open(data_path + "cumulative_API_set.csv", "w") as write_to:
    writer = csv.writer(write_to)
    writer.writerow(['organism', 'gene_id', 'orthologType', 'accession', 'chromosome', 'start', 'end', 'direction (+/-)', 'intactness', 'sequence']) #write header row w orthologous sequence data fields 
    for file in os.listdir(data_path):
        if file[:4] == "temp":
            file_path = data_path + file
            with open(file_path) as read_from:
                for line in read_from:
                    line = line.split(",")
                    if len(line) != 10 or line[0] == 'organism': continue 
                    total_data += 1
                    #check for lines w slightly altered coords 
                    identifier = line[:-1]
                    start, end = identifier[5], identifier[6]
                    alt_1, alt_2, alt_3, alt_4 = identifier, identifier, identifier, identifier
                    alt_1[5], alt_2[5] = str(int(start)-1), str(int(start)+1)
                    alt_3[6], alt_4[6] = str(int(end)+1), str(int(end)-1)
                    alt_1, alt_2, alt_3, alt_4 = tuple(alt_1), tuple(alt_2), tuple(alt_3), tuple(alt_4)
                    #modify start/end coord of alt lines s.t we are checking for dupe lines w +1/-1 coordinate differences due to small script change   

                    if tuple(line[:-1]) in cumulative_set:
                        uniq_data += 1
                        writer.writerow(line)
                        cumulative_set.remove(tuple(line[:-1]))
                        #remove line and all functionally duplicate (but slightly diff identified) lines s.t we don't add duplicate data to final set
                        cumulative_set.discard(alt_1)
                        cumulative_set.discard(alt_2)
                        cumulative_set.discard(alt_3)
                        cumulative_set.discard(alt_4)
                    else: continue 

print("total data written:", total_data)
print("unique data written:", uniq_data)



         




