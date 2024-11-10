import os, csv, sys


"""
Script for iterating through training results for per-gene enformer approach
Loops through directories (corresponding to genes in: 
/data/rbg/users/wmccrthy/chemistry/Everything/fall_24_training/gene_training_metrics/


"""

RESULTS_DIR = "/data/rbg/users/wmccrthy/chemistry/Everything/fall_24_training/gene_training_metrics/"
TRAINING_DIR = "/data/rbg/users/wmccrthy/chemistry/Everything/fall_24_training/"

"""
Given a gene's training results directory, scrapes the metric file to return:
    avg train loss, avg val loss, min val loss
"""
def scrape_metrics(path):
    train_losses, val_losses = [], []
    i = 0
    val_min = float("inf")
    with open(path + "/version_0/metrics.csv") as read_from:
        for line in read_from:
            line = line.split(",")
            # for now we don't care abt epochs, we care abt averaging train_loss and val_loss
            train_loss, val_loss = line[2:]
            if i > 0:
                if len(train_loss) > 0: train_losses.append(float(train_loss))
                elif len(val_loss) > 0: val_losses.append(float(val_loss))
            i += 1
    train_loss = sum(train_losses)/len(train_losses)
    val_loss = sum(val_losses)/len(val_losses)
    val_min = min(val_losses)
    return train_loss, val_loss, val_min

"""
Iterates through gene result directories and outputs organized CSV with all results.
"""
def compile_results():
    with open(TRAINING_DIR + "per_gene_results.csv", "w") as write_to:
        writer = csv.writer(write_to)
        writer.writerow(["gene", "avg train loss", "avg val loss", "min val loss"])
        for sub_dir in os.listdir(RESULTS_DIR):
            gene = sub_dir.split("_")[0]
            avg_train, avg_val, val_min = scrape_metrics(RESULTS_DIR + "/" + sub_dir)
            print(gene, avg_train, avg_val, val_min)
            writer.writerow([gene, avg_train, avg_val, val_min])

"""
Returns list of tuples of results by parsing "per_gene_results.csv" 
Format is: List[(gene, train_avg_loss, val_avg_loss, val_min_loss)]
"""
def get_results():
    results = []
    # iterate thru organized results CSV to create dataset
    with open(TRAINING_DIR + "per_gene_results.csv") as read_from:
        for line in read_from:
            line = line.split(",")
            if line[0] == "gene": continue #skip first row
            results.append(line) # add line to results list (format: gene, train_avg, val_avg, val_min)
    return results
        







            

