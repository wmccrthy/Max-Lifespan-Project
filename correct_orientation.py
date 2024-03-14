"""
SOMEHOW NEED TO PERFORM A TEST TO SEE HOW DNA ZOO 2BIT FILES STORE SEQUENCES S.T I CAN FIGURE OUT HOW TOGA DERIVED NEGATIVE (-) STRANDS FROM SUCH GENOMES 

    - PICK DATA POINT MANUALLY (I KNOW SEQUENCE, I KNOW SPECIES, ETC)
    - query for sequence from corresponding 2bit file
    - see if sequence aligns with that recorded by TOGA 

    chosen data point: Phalanger_gymnotis__-__HLphaGym1,ENST00000382102.SLC45A2.7,one2one,HiC_scaffold_1,323568230,323622859,-,intact, ...

    query for sequence: Desktop/twoBitToFa /Users/wyattmccarthy/Desktop/'MIT Aging Project'/Pre-Processing/'aids test'/HLphaGym1.2bit:HiC_scaffold_1:323568230-323622859 /Users/wyattmccarthy/Desktop/'MIT Aging Project'/Pre-Processing/2bit_Sequences/HLphaGym1.fa

    CONCLUSION: DNA ZOO 2bit files and UCSC genomes are all in + orientation, when TOGA has a - strand from such a file/database with coordinates x - y, 
                that corresponds to the reverse complement of the + strand from x - y; since our mining scripts got + sequence data
                from y - y + 1000 and appended it to end of TOGA-given - strand, we can correct data collected as such by taking the reverse complement of 
                the last 1000 bp of - sequences and inserting it at front of sequence (remove from end, insert at beginning)
"""

import os, sys, Bio, csv 
from Bio.Seq import Seq

"""
SCRIPT FOR PERFORMING PROCEDURE DESCRIBED ABOVE IN "CONCLUSION"

PSUEDO: 
for each file in os.listdir(dataSubSets_path):
    with open(dataSubSets_path/file, "r") as read_from:
        with open(dataSubSets_path/file+"corrected") as write_to:
            for line in read_from:
                line = line.split(,)
                if line[orientation] != '-':        data is correct so write it to replacement file 
                    write_to.write(line)
                else: 
                    TOGA_provided = line[-1][:-1000]
                    sub_sequence_to_modify = line[-1][-1000:]
                    sub_sequence_to_modify = str(Seq(sub_sequence_to_modify).reverse_complement())
                    corrected_sequence = sub_sequence_to_modify + TOGA_provided 
                    corrected_line = line[:-1] + [corrected_sequence]
                    print(line)
                    print(corrected_line)
                    write_to.write(corrected_line)
"""
data_path_MIT = f"/Mounts/rbg-storage1/users/wmccrthy/dataSubSets/" 
corrupted_path_MIT = f"/Mounts/rbg-storage1/users/wmccrthy/corrupted_data/"
total_data = 0
corrected_data = 0

for file in os.listdir(corrupted_path_MIT):
    if file[:3] == 'API': continue #all valid API data at this point is aggregated into temp_API_subset so just correct data in that file 

    file_path = corrupted_path_MIT + file
    new_file_path = data_path_MIT + file[:-4] + "_corrected.csv"

    # if "corrected" in file: continue #if file already exists (we have already collected data), don't incur redundant time 

    # if os.system(f'test -e {new_file_path}') == 0 or "corrected" in file: continue #if file already exists (we have already collected data), don't incur redundant time 

    print(new_file_path, " for ", file_path)
    with open(file_path) as read_from:
        with open(new_file_path, "w") as write_to_file:
            write_to  = csv.writer(write_to_file)
            for line in read_from:
                total_data += 1
                line = line.split(",")
                if len(line) < 3: 
                    print(line)
                    continue #avoid ed data (weird)
                if line[-3] != '-':        #data is correct so write it to replacement file 
                    write_to.writerow(line)
                else:
                    corrected_data += 1
                    TOGA_provided = line[-1][:-1000]
                    sub_sequence_to_modify = line[-1][-1000:]
                    sub_sequence_to_modify = str(Seq(sub_sequence_to_modify).reverse_complement())
                    corrected_sequence = sub_sequence_to_modify + TOGA_provided
                    corrected_sequence = corrected_sequence.replace("\n", "") #get rid of fucking new line chars at beginning, totally fucking up data parsing 
                    corrected_line = line[:-1] + [corrected_sequence]
                    # print(line)
                    # print(corrected_line)
                    write_to.writerow(corrected_line)

print("total data rescribed: ", total_data)
print("total data corrected: ", corrected_data)



