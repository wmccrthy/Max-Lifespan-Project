import os, csv, sys, numpy as np
import matplotlib.pyplot as plt


"""
Script for iterating through training results for per-gene enformer approach
Loops through directories (corresponding to genes in: 
/data/rbg/users/wmccrthy/chemistry/Everything/fall_24_training/gene_training_metrics/


"""
RESULTS_DIR = "/data/rbg/users/wmccrthy/chemistry/Everything/fall_24_training/gene_training_metrics"
TRAINING_DIR = "/data/rbg/users/wmccrthy/chemistry/Everything/fall_24_training/"
GENE_LOSS_DIR =  "/data/rbg/users/wmccrthy/chemistry/Everything/fall_24_training/loss_plots/"

"""
Given a gene's training results directory, scrapes the metric file to return:
    avg train loss, avg val loss, min val loss
"""
def scrape_metrics(path, per_epoch = False):
    train_losses, val_losses = [], []
    i = 0
    val_min = float("inf")
    # check if exists
    if not os.path.exists(path + "/version_0/metrics.csv"): return 999, 999, 999, 999

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
    train_min = min(train_losses)
    if per_epoch: return train_losses, val_losses
    return train_loss, train_min, val_loss, val_min

"""
Iterates through gene result directories and outputs organized CSV with all results.
"""
def compile_results(tag = None):
    if tag: tag = "_" + tag
    else: tag = ""
    results_dir = TRAINING_DIR + f"per_gene_results{tag}.csv"
    with open(results_dir, "w") as write_to:
        writer = csv.writer(write_to)
        writer.writerow(["gene", "avg train loss", "min train loss", "avg val loss", "min val loss"])
        for sub_dir in os.listdir(RESULTS_DIR + tag + "/"):
            gene = sub_dir.split("_")[0]
            avg_train, train_min, avg_val, val_min = scrape_metrics(RESULTS_DIR + tag + "/" + sub_dir)
            print(gene, avg_train, train_min, avg_val, val_min)
            writer.writerow([gene, avg_train, train_min, avg_val, val_min])

"""
Returns list of tuples of results by parsing "per_gene_results.csv" 
Format is: List[(gene, train_avg_loss, val_avg_loss, val_min_loss)]
"""
def get_results(is_dict = False, tag = None):
    results = []
    if is_dict: results = {}
    if tag: tag = "_" + tag
    else: tag = ""
    results_dir = TRAINING_DIR + f"per_gene_results{tag}.csv"
    # iterate thru organized results CSV to create dataset
    with open(results_dir) as read_from:
        for line in read_from:
            line = line.split(",")
            if line[0] == "gene": continue #skip first row
            if not is_dict: results.append(line) # add line to results list (format: gene, train_avg, train_min, val_avg, val_min)
            else: results[line[0]] = (line[1:]) # add line to results list (format: gene, train_avg, train_min, val_avg, val_min)
    return results

def results_by_train_loss():
    arr = sorted(get_results(), key = lambda x:float(x[2])) # sorts results by min train loss
    for i in arr: print(i)


"""
Plot min val loss by gene
"""
def plot_min_val(tag = None):
    training_results = get_results(tag = tag)

    # trim training results so we're only plotting best genes (those w val loss < 10)
    training_results = [i for i in training_results if float(i[-1]) <= 10]

    labels, values = [i[0] for i in training_results], [float(i[-1]) for i in training_results]
    plt.figure(figsize=(14, 8))
    plt.bar(labels, values)

    # Set labels and title
    plt.xlabel("Gene")
    plt.ylabel("Min Validation Loss")
    plt.title("Minimum Validation Loss per Gene")

    # Rotate gene names if needed for readability
    plt.xticks(rotation=45, ha="right")

    # Add grid for better readability
    plt.grid(visible=True, axis='y', linestyle='--', alpha=0.7)

    # Display the plot
    # plt.tight_layout()
    if tag: tag = "_" + tag
    else: tag = ""
    plt.show()
    plt.savefig(f"per_gene_min_val_loss{tag}.png")

"""
Plot min training loss by gene
"""
def plot_min_train(tag = None):
    training_results = get_results(tag = tag)

    # trim training results so we're only plotting best genes (those w min training loss <= 7)
    training_results = [i for i in training_results if float(i[2]) <= 7]

    labels, values = [i[0] for i in training_results], [float(i[2]) for i in training_results]
    plt.figure(figsize=(14, 8))
    plt.bar(labels, values)

    # Set labels and title
    plt.xlabel("Gene")
    plt.ylabel("Min Training Loss")
    plt.title("Minimum Training Loss per Gene")

    # Rotate gene names if needed for readability
    plt.xticks(rotation=45, ha="right")

    # Add grid for better readability
    plt.grid(visible=True, axis='y', linestyle='--', alpha=0.7)

    # Display the plot
    # plt.tight_layout()
    if tag: tag = "_" + tag
    else: tag = ""
    plt.show()
    plt.savefig(f"per_gene_min_train_loss{tag}.png")

"""
method to get results for a given gene and tag
"""
def get_results(gene, tag = None):
    train_losses, val_losses = [], []
    if tag: tag = "_" + tag
    else: tag = ""
    for sub_dir in os.listdir(RESULTS_DIR + tag + "/"):
        cur_gene = sub_dir.split("_")[0]
        if cur_gene == gene:
            train_losses, val_losses = scrape_metrics(RESULTS_DIR + tag + "/" + sub_dir, per_epoch = True)
            break
    return train_losses, val_losses

"""
Method to compare gene results btwn two different runs/tags
"""
def compare_gene_loss(gene, tag1, tag2 = None):
    gene_tag1_train, gene_tag1_val = get_results(gene, tag1)
    gene_tag2_train, gene_tag2_val = get_results(gene, tag2)
    
    if tag1: tag1 = "_" + tag1
    else: tag1 = ""
    if tag2: tag2 = "_" + tag2
    else: tag2 = "_default"
    plt.figure(figsize=(10, 6))  # Set figure size
    plt.plot([i for i in range(len(gene_tag1_train))], gene_tag1_train, marker='o', linestyle='-', color='blue', label = f"training{tag1}")  # Plot with markers
    plt.plot([i for i in range(len(gene_tag1_val))], gene_tag1_val, marker='o', linestyle='-', color='green', label = f"validation{tag1}")  # Plot with markers
    plt.plot([i for i in range(len(gene_tag2_train))], gene_tag2_train, marker='o', linestyle='-', color='red', label = f"training{tag2}")  # Plot with markers
    plt.plot([i for i in range(len(gene_tag2_val))], gene_tag2_val, marker='o', linestyle='-', color='orange', label = f"validation{tag2}")  # Plot with markers
   
    plt.legend(loc = "upper right")

    # Set labels and title
    plt.xlabel("Epoch")
    plt.ylabel("Loss")
    plt.title("Training and Validation Loss for " + gene)

    # Add grid for better readability
    plt.grid(visible=True)

    # Display the plot
    plt.tight_layout()
    plt.show()
    plt.savefig(GENE_LOSS_DIR + gene + f"_loss{tag1}_vs{tag2}.png")
    plt.close()

"""
Method to plot train loss and val loss for given gene (one plot per gene)
"""
def plot_gene_loss(gene, path = None, tag = None):
    if not tag: tag = ""
    # get results for a given gene
    train_losses, val_losses = [], []
    if not path:
        for sub_dir in os.listdir(RESULTS_DIR + tag + "/"):
                cur_gene = sub_dir.split("_")[0]
                if cur_gene == gene:
                    train_losses, val_losses = scrape_metrics(RESULTS_DIR + tag + "/" + sub_dir, per_epoch = True)
                    break
    else:
        train_losses, val_losses = scrape_metrics(path, per_epoch = True)

    plt.figure(figsize=(10, 6))  # Set figure size
    plt.plot([i for i in range(len(train_losses))], train_losses, marker='o', linestyle='-', color='blue', label = "training")  # Plot with markers
    plt.plot([i for i in range(len(val_losses))], val_losses, marker='o', linestyle='-', color='green', label = "validation")  # Plot with markers
    plt.legend(loc = "upper right")
    # Set labels and title
    plt.xlabel("Epoch")
    plt.ylabel("Loss")
    plt.title("Training and Validation Loss for " + gene + f" ({tag})")

    # Add grid for better readability
    plt.grid(visible=True)

    # Display the plot
    plt.tight_layout()
    plt.show()
    plt.savefig(GENE_LOSS_DIR + gene + f"_loss{tag}.png")
    plt.close()

"""
Iterate through results dir and plot train/val loss for each gene
"""
def plot_all_gene_losses(tag = None):
    if tag: tag = "_" + tag
    else: tag = ""
    for sub_dir in os.listdir(RESULTS_DIR + tag + "/"):
            gene = sub_dir.split("_")[0]
            plot_gene_loss(gene, path = RESULTS_DIR + tag + "/" + sub_dir, tag = tag)

if __name__ == "__main__":
    args = sys.argv
    print(args)
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])







            

