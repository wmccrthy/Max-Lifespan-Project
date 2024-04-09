"""
script w method that takes in csv file path and prompts for user input query
supports basic queries (return all rows corresponding to certain conditions)
call: python3 query_data.py prompt {file_path}
"""
import csv, os, sys 

"""
given a csv file allows user to query for entries from that csv that conform to given value constraints 
"""
def prompt(file):
    fields2ind = set()
    ind2fields = {}
    with open(file) as read_from:
        for line in read_from:
            line = line.split(',')
            fields2ind = {i[1].replace("\n", ""):i[0] for i in enumerate(line)}
            break
    ind2fields = {fields2ind[i]:i for i in fields2ind}
    print(ind2fields)
    query = input("Give query below in format col#:lower_bound, upper_bound ex 2:.2,.3 4:50,9999 \n")
    conds = query.split(" ")
    conds = [i.split(":") for i in conds]
    for i in range(len(conds)):
        conds[i][0] = int(conds[i][0])
        conds[i][1] = [float(i) for i in conds[i][1].split(",")]
    response = handle_query(file, conds, len(fields2ind.keys()))
    for i in response: print(i)
    return

def handle_query(file, conds, fields_num):
    to_ret = []
    line_num = 0
    with open(file) as read_from:
        for line in read_from:
            line = line.split(",")
            if len(line) != fields_num: continue 

            if line_num == 0:
                line_num += 1
                continue 
            valid = True
            for c in conds:
                col, constraint = c
                if float(line[col]) < constraint[0] or float(line[col]) > constraint[1]: 
                    valid = False
                    break
            if valid: to_ret.append(line)
    return to_ret


"""
used to slightly modify data in regulatory bins data files s.t bins are denoted like (x-y) as opposed to (x,y) since the latter makes interferes with parsing
"""
def fix_bins_data(file):
    ind = 0
    with open(file[:-4] + '_fixed.csv', 'w') as write_to:
        writer = csv.writer(write_to)
        with open(file) as read_from:
            for line in read_from:
                    line = line.split(',')
                    if ind == 0: writer.writerow(line)
                    else:
                        #columns 1 and 2 for regulatory_interesting
                        # line = line[:1] + ["-".join(line[1:3]).replace('"', '').replace(" ", "")] + [i.replace("\n", "") for i in line[3:]]
                        # print(line)

                        #cols 2 and 3 for bins_chi_square 
                        line = line[:2] + ["-".join(line[2:4]).replace('"', '').replace(" ", "")] + [i.replace("\n", "") for i in line[4:]]

                        writer.writerow(line)
                    ind += 1
    print(ind)


if __name__ == "__main__":
    args = sys.argv
    print(args)
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])
