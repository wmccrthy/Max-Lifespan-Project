import csv, os, sys 
"""
GIVEN THAT THE UCSC GENOME BROWSER WENT DOWN OVER THE WEEKEND, MANY OF THE QUERIES MADE BY MY API-BASED SCRIPTS RETURNED NULL DATA (EVEN IF QUERY WAS VALID)
BECAUSE OF THAT, I HAVE TO RE-RUN API-BASED SCRIPTS S.T WE CAN RETRIEVE THE VALID DATA THAT AN INCAPACITATED GENOME BROWSER FAILED TO GIVE US 
WITH THAT BEING SAID, THE API-BASED SCRIPTS DID RETRIEVE SOME VALID DATA (FROM ALL QUERIES MADE BETWEEN RUN-START-TIME AND THE MOMENT GENOME BROWSER WENT DOWN)
I WOULD STILL LIKE TO HAVE ACCESS TO/USE THIS DATA WHILE WE WAIT FOR THE RE-RUN OF API-BASED SCRIPTS TO FINISH (MIGHT TAKE 7-8 DAYS)

AS SUCH, THIS IS A SCRIPT FOR ITERATING THROUGH ALL RETRIEVED API DATA AND EXTRACTING THAT WHICH IS NOT NULL (and compiling it into a dataset)

PSUEDO: 

with open({output_path}, 'w') as write_to_file:
    writer = csv.writer(write_to_file)

    writer.writerow(['organism', 'gene_id', 'orthologType', 'accession', 'chromosome', 'start', 'end', 'direction (+/-)', 'intactness', 'sequence']) 
    #write header row w orthologous sequence data fields 

    for each file in os.listdir({dataSubSet_path}):
        if file[:3] = "API":
            with open(file) as data_file:
                for line in data_file:
                    line = line.split(",")
                    if line[-1] != "Not Queryable": writer.writerow(line)

"""
data_path_MIT = f"/Mounts/rbg-storage1/users/wmccrthy/dataSubSets/" 
API_subset_num = 1

for file in os.listdir(data_path_MIT): #ensure subset file created is tagged correctly (so it doesn't overwrite prior dataset files)
    if file[:4] == "temp": API_subset_num += 1

not_null_data_cnt = 0
total_data_cnt = 0 
not_queryable_data_cnt = 0
coordinate_errors = 0

output_path_MIT = f"/Mounts/rbg-storage1/users/wmccrthy/dataSubSets/temp_API_subSet_{API_subset_num}.csv" 
print(output_path_MIT)

with open(f'{output_path_MIT}', 'w') as write_to_file:
    writer = csv.writer(write_to_file)
    writer.writerow(['organism', 'gene_id', 'orthologType', 'accession', 'chromosome', 'start', 'end', 'direction (+/-)', 'intactness', 'sequence']) 
    #write header row w orthologous sequence data fields 
    for file in os.listdir(data_path_MIT):
        if file[:3] == 'API':
            print(file)
            with open(f'{data_path_MIT}/{file}') as data_file:
                for line in data_file:
                    line = line.split(",")
                    cur_org = line[0]
                    if cur_org == "organism": continue 

                    total_data_cnt += 1
                    if line[-1] != "Not Queryable" and "Not Queryable" not in line[-1] and len(line) == 10: 
                        not_null_data_cnt += 1
                        line[-1] = line[-1].replace("\n", "").replace('"', '')
                        # print(line)
                        writer.writerow(line)
                        if len(line[-1]) < 50: 
                            if len(line) < 10: print(line)
                            coordinate_errors += 1
                    elif len(line) == 10 and "Not Queryable" in line[-1]: 
                        not_queryable_data_cnt += 1

#debugging/stat output 
print("Total Data Points:", total_data_cnt)
print("Data Points w Full Sequences (not null):", not_null_data_cnt)
print("Not Queryable Data Points:", not_queryable_data_cnt)
print("Coordinate Errors", coordinate_errors)




