"""
script w method that takes in csv file path and prompts for user input query
supports basic queries (return all rows corresponding to certain conditions)
call: python3 query_data.py prompt {file_path}
"""
import csv, os, sys, numpy as np, math
from scipy import stats
import matplotlib.pyplot as plt 

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


"""
GET AVERAGE SEQUENCE LENGTH (RAW AND TOKENIZED) ACROSS ALL GENE SETS BY ITERATING THRU METADATA CSV'S 
"""
def get_avg_lengths():
    avg_raw, avg_tokenized = 0, 0
    num_raw, num_tokenized = 0, 0
    avg_species_rep = 0
    with open("/Users/wyattmccarthy/Desktop/MIT Aging Project/Everything/EDA/gene_sets_data/regulatory_sets_metadata.csv") as read_from:
        for line in read_from:
            line = line.split(",")
            if num_raw > 0: 
                avg_raw += float(line[-3])
                avg_species_rep += float(line[-1])
            num_raw += 1
            
    avg_raw /= (num_raw-1)
    avg_species_rep /= (num_raw-1)

    with open("/Users/wyattmccarthy/Desktop/MIT Aging Project/Everything/EDA/gene_sets_data/regulatory_sets_tokenized_lengths.csv") as read_from:
        for line in read_from:
            line = line.split(",")
            if num_tokenized > 0:
                avg_tokenized += float(line[1])
            num_tokenized += 1
    avg_tokenized /= num_tokenized

    print("Average Sequence (raw) length:", avg_raw)
    print("Average Species Represented per set:", avg_species_rep)
    print("Average Max Sequence Length (tokenized)", avg_tokenized)


regulatory_sets_path = "/data/rsg/chemistry/wmccrthy/Everything/gene_datasets/regulatory/"

"""
GET AVG # OF SEQUENCES PER GENE SET 
"""
def average_per_set():
    file_seq_cnts = []
    for file in os.listdir(regulatory_sets_path):
        file_path = regulatory_sets_path + file
        if 'embedding' not in file:
            with open(file_path) as read_from:
                file_cnt = 0
                for line in read_from:
                    if len(line.split(",")) >= 9: file_cnt += 1
            file_seq_cnts.append(file_cnt)

    print("Sequence Count per file data:", max(file_seq_cnts), min(file_seq_cnts), sum(file_seq_cnts)/len(file_seq_cnts))


"""
CREATE DICTS: 
    - GENUS:[ANIMALS IN GENUS]
    - FAMILY:[ANIMALS IN FAMILY]
WRITES OUT GENUSES/FAMILIES WITH LIFESPAN DISPARITY 

PLOTS HISTOGRAM OF LIFESPANS PER FAMILY, GENUS (perhaps add order too)
"""
def find_lifespan_anomalies():
    genera, families, orders = {}, {}, {}

    with open("/Users/wyattmccarthy/Desktop/MIT Aging Project/Everything/EDA/lifespans_by_order.csv") as read_from:
        for line in read_from:
            line = line.split(",")
            if line[0] == "order" or line[-1].replace("\n", "") == 'Unknown': continue 
            order, family, genus = line[:3]
            if len(genus.strip()) == 0: continue 
            species, lifespan = line[-2], float(line[-1])
            if order not in orders: orders[order] = [(species, lifespan)]
            else: orders[order].append((species, lifespan))
            if family not in families: families[family] = [(species, lifespan)]
            else: families[family].append((species, lifespan))
            if genus not in genera: genera[genus] = [(species, lifespan)]
            else: genera[genus].append((species, lifespan))
    
    with open('orders_of_interest_IQR.csv', 'w') as write_to:
        writer = csv.writer(write_to)
        # writer.writerow(['order', 'median max lifespan', 'MAD', 'outlying species', 'outlying max lifespan'])
        writer.writerow(['order', 'median lifespan', 'IQR', 'outlying species', 'outlying max lifespan'])
        for order in orders:
            lifespan2species = {i[1]:i[0] for i in orders[order]}
            order_lifespans = np.array([i[1] for i in orders[order]])
            median_l = np.median(order_lifespans)
            MAD = stats.median_abs_deviation(order_lifespans)
            q1, q3 = np.percentile(order_lifespans, 25), np.percentile(order_lifespans, 75)
            iqr = q3 - q1
            lower, upper = q1 - 1.5*iqr, q3 + 1.5*iqr 
            if len(order_lifespans) > 2:
                for i in range(len(order_lifespans)):
                    lifespan = order_lifespans[i]
                    if lifespan < lower or lifespan > upper: writer.writerow([order, median_l, iqr, lifespan2species[lifespan], lifespan])
                # tmp = 0.6745*(order_lifespans - median_l)/MAD
                # for i in range(len(tmp)):
                #     if tmp[i] >= 2: writer.writerow([order, median_l, MAD, lifespan2species[order_lifespans[i]], order_lifespans[i]])

                fig = plt.figure()
                grph = fig.add_subplot()
                grph.hist(order_lifespans)
                grph.set_title(f"{order} lifespans (n={len(order_lifespans)})")
                # [grph.text(order_lifespans[i], i/100, s=lifespan2species[order_lifespans[i]]) for i in range(len(order_lifespans))] hard to label clearly given how lifespans are 'binned' 
                fig.savefig(f"/Users/wyattmccarthy/Desktop/MIT Aging Project/Everything/EDA/lifespan_histograms/orders/{order}_hist.png")
                plt.close()

            # else:
            #     avg_l = order_lifespans.mean()
            #     std_l = order_lifespans.std()
            #     for l in order_lifespans:
            #       if l - avg_l > 1.5  * std_l: writer.writerow([order, avg_l, std_l, lifespan2species[l], l])

    with open("genera_of_interest_IQR.csv", "w") as write_to:
        writer = csv.writer(write_to)
        # writer.writerow(['genus', 'median max lifespan', 'MAD', 'outlying species', 'outlying max lifespan'])
        writer.writerow(['genus', 'median lifespan', 'IQR', 'outlying species', 'outlying max lifespan'])
        for genus in genera:
            lifespan2species = {i[1]:i[0] for i in genera[genus]}
            genus_lifespans = np.array([i[1] for i in genera[genus]])
            median_l = np.median(genus_lifespans)
            MAD = stats.median_abs_deviation(genus_lifespans)
            q1, q3 = np.percentile(genus_lifespans, 25), np.percentile(genus_lifespans, 75)
            iqr = q3 - q1
            lower, upper = q1 - 1.5*iqr, q3 + 1.5*iqr
            if len(genus_lifespans) > 2:
                for i in range(len(genus_lifespans)):
                    lifespan = genus_lifespans[i]
                    if lifespan < lower or lifespan > upper: writer.writerow([genus, median_l, iqr, lifespan2species[lifespan], lifespan])
                fig = plt.figure()
                grph = fig.add_subplot()
                grph.hist(genus_lifespans)
                grph.set_title(f"{genus} lifespans (n={len(genus_lifespans)})")
                fig.savefig(f"/Users/wyattmccarthy/Desktop/MIT Aging Project/Everything/EDA/lifespan_histograms/genera/{genus}_hist.png")
                plt.close()
            #MAD BASED OUTLIER COMP
            #     tmp = 0.6745* (genus_lifespans - median_l)/MAD 
            #     # print(genus, genus_lifespans, tmp)
            #     for i in range(len(tmp)):
            #         if tmp[i] >= 2: 
            #             writer.writerow([genus, median_l, MAD, lifespan2species[genus_lifespans[i]], genus_lifespans[i]])
            #             continue 
            # else:
            #     avg_l = genus_lifespans.mean()
            #     std_l = genus_lifespans.std()
            #     for l in genus_lifespans:
            #         if l - avg_l > 1.05 * std_l: writer.writerow([genus, avg_l, std_l, lifespan2species[l], l])
                

    with open("families_of_interest_IQR.csv", "w") as write_to:
        writer = csv.writer(write_to)
        # writer.writerow(['family', 'median max lifespan', 'MAD', 'outlying species', 'outlying lifespan'])
        writer.writerow(['family', 'median lifespan', 'IQR', 'outlying species', 'outlying max lifespan'])
        for family in families:
            lifespan2species = {i[1]:i[0] for i in families[family]}
            family_lifespans = np.array([i[1] for i in families[family]])
            median_l = np.median(family_lifespans)
            MAD = stats.median_abs_deviation(family_lifespans)
            q1, q3 = np.percentile(family_lifespans, 25), np.percentile(family_lifespans, 75)
            iqr = q3 - q1
            lower, upper = q1 - 1.5*iqr, q3 + 1.5*iqr
            if len(family_lifespans) > 2:
                for i in range(len(family_lifespans)):
                    lifespan = family_lifespans[i]
                    if lifespan < lower or lifespan > upper: writer.writerow([family, median_l, iqr, lifespan2species[lifespan], lifespan])
                fig = plt.figure()
                grph = fig.add_subplot()
                grph.set_title(f"{family} lifespans (n={len(family_lifespans)})")
                fig.savefig(f"/Users/wyattmccarthy/Desktop/MIT Aging Project/Everything/EDA/lifespan_histograms/families/{family}_hist.png")
                plt.close()
            #MAD BASED OUTLIERS
            #     tmp = 0.6745* (family_lifespans - median_l)/MAD 
            #     # print(family, family_lifespans, tmp)
            #     for i in range(len(tmp)):
            #         if tmp[i] >= 2: writer.writerow([family, median_l, MAD, lifespan2species[family_lifespans[i]], family_lifespans[i]])
            # else:
            #     avg_l = family_lifespans.mean()
            #     std_l = family_lifespans.std()
            #     for l in family_lifespans:
            #         if l - avg_l > 1.25 * std_l: writer.writerow([family, avg_l, std_l, lifespan2species[l], l])

    return 


"""
SCRIPT FOR RANKING GENES BASED ON avg stat sig of clusters
ITERATES THRU GIVEN CLUSTER DATA, COMPUTES AVG stat sig of clusters, OUTPUTS FOLLOWING FORMAT:
[gene, avg cluster stat sig, most stat sig, least stat sig, ranking score]

genes_dict will hold gene:[total stat sig, total clusters, max stat sig, min stat sig] | recall that max stat sig equates to least significant and min stat sig equates to most significant

ranking score is computed as 1 / (avg_stat_sig * # of species in clusters)
"""
def rank_genes(cluster_data_path):
    col_to_avg = -1
    genes = {}
    with open(cluster_data_path) as read_from:
        for line in read_from:
            line = line.split(",")
            if line[0] == "gene type": continue 
            # if not col_to_avg:
            #     for i in range(len(line)): 
            #         if 'statistical' in line[i]: 
            #             col_to_avg = i
            #             continue 

            gene, stat_sig = line[0], float(line[col_to_avg])
            if stat_sig == 0: stat_sig += .0001
            num_species_in_cluster = int(line[-3])
            if gene in genes: 
                genes[gene][0] += stat_sig
                genes[gene][1] += 1
                genes[gene][2], genes[gene][3] = max(genes[gene][2], stat_sig), min(genes[gene][3], stat_sig)
                genes[gene][4] += 1/(stat_sig/math.sqrt(num_species_in_cluster))
                genes[gene][5] += num_species_in_cluster
            else: genes[gene] = [stat_sig, 1, stat_sig, stat_sig, 1/(stat_sig / math.sqrt(num_species_in_cluster)), num_species_in_cluster]
    
    with open("gene_rankings.csv", "w") as write_to:
        writer = csv.writer(write_to)
        writer.writerow(['gene', 'avg cluster stat sig', 'most stat sig cluster', 'least stat sig cluster', 'avg # species per cluster'])
        for g in genes:
            stat_sig_avg = genes[g][0] / genes[g][1]
            ranking_score_avg = genes[g][4] / genes[g][1]
            avg_species_per_cluster = genes[g][5] / genes[g][1]
            writer.writerow([g, stat_sig_avg, genes[g][3], genes[g][2], avg_species_per_cluster])


if __name__ == "__main__":
    args = sys.argv
    print(args)
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])
