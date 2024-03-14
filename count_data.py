import os, sys, Bio, csv 
from Bio.Seq import Seq

"""
SIMPLE SCRIPT FOR COUNTING THE DATA COLLECTED SO FAR (USEFUL UTILITY)
"""

data_path_MIT = f"/Mounts/rbg-storage1/users/wmccrthy/dataSubSets/" 
total_data = 0
total_data_w_seq = 0 

for file in os.listdir(data_path_MIT):
    file_path = data_path_MIT + file
    if file[:3] == "API": continue 
    with open(file_path) as read_from:
        for line in read_from:
            line = line.split(",")
            if len(line) == 10: 
                total_data += 1
                if len(line[-1]) > 50: total_data_w_seq += 1

print("Total Data Points:", total_data)
print("Total Data Points w (not null) Sequences:",total_data_w_seq)
