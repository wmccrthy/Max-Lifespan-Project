import os, sys, csv 

"""
SIMPLE SCRIPT FOR COUNTING THE DATA COLLECTED SO FAR (USEFUL UTILITY)
"""
def count_seqs(file_name):
    data_path_MIT = f"/data/rbg/users/wmccrthy/rbgquanta1/wmccrthy/" 
    total_data = 0
    total_data_w_seq = 0
    file_path = data_path_MIT + file_name
    print(file_path)
    with open(file_path) as read_from:
        for line in read_from:
            line = line.split(",")
            if len(line) == 9: 
                total_data += 1
                if len(line[-1]) > 50: total_data_w_seq += 1

    print("Total Data Points:", total_data)
    print("Total Data Points w (not null) Sequences:",total_data_w_seq)

if __name__ == "__main__":
    args = sys.argv
    print(args)
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])