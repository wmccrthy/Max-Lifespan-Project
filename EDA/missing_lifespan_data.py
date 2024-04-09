"""
SCRIPT FOR QUICKLY COMBING THRU LIFESPAN DATA TO GET IDEA OF HOW MANY ENTRIES ARE MISSING AND FOR WHAT SPECIES 
    - as of now best bet is to manually find missing lifespan data 
"""
import csv, os, sys 
from lifespanScraping import scrape_family

lifespan_data = "/Users/wyattmccarthy/Desktop/MIT Aging Project/Everything/EDA/lifespan_data.csv"


def find_missing():
    total_species = 0
    missing_data = 0
    missing = set()
    present = set()
    with open(lifespan_data) as read_from:
        for line in read_from:
            line = line.split(",")
            if line[0] == "Species": continue 
            
            total_species += 1
            lifespan = line[1]
            if lifespan == "Unknown" or lifespan == "Not Established":
                missing_data += 1
                missing.add(line[0])
    print("Total Species:", total_species)
    print("Missing Lifespan for:", missing_data)
    print("Have Lifespan for:", total_species - missing_data)

    # with open("missing_lifespans.csv", "w") as write_to:
    #     writer = csv.writer(write_to)
    #     writer.writerow(["Species", "Max Lifespan", "Source"])
    #     for species in missing:
    #         writer.writerow([species, 'N/A', 'N/A'])

    return 


"""
WRITES OUT MAX LIFESPANS W ORDER FIELD S.T WE HAVE MEANS OF FINDING 'SIMILAR' ANIMALS W VARYING LIFESPANS 
"""
def get_lifespans_per_order():
    #build species->max lifespan map 
    lifespan_path = '/data/rsg/chemistry/wmccrthy/Everything/lifespan_data.csv'
    lifespan_mappings = {}
    
    with open(f'{lifespan_path}') as file:
        for line in file:
            line = line.strip().split(",")
            organism = line[0]
            lifespan = line[1]
            lifespan_mappings[organism] = lifespan

    #iterate thru TOGA directory structure, for each order write out species, order, lifespan
    toga_home_dir = "/data/rsg/chemistry/wmccrthy/Everything/human_hg38_reference"
    written = set()
    with open("lifespans_by_order.csv", "w") as write_to:
        writer = csv.writer(write_to)
        writer.writerow(['order', 'species', 'max lifespan'])
        for order in os.listdir(toga_home_dir):
            order_path = toga_home_dir + f"/{order}"
            for organism in os.listdir(order_path):
                #write out line w order | species | max lifespan
                species = "_".join(organism.split("_")[:2])
                if species not in written and species in lifespan_mappings:
                    written.add(species)
                    family, genus = scrape_family(species)
                    writer.writerow([order, family, genus, species, lifespan_mappings[species]])
    
    return 

#NEDD4L in top ssh 



if __name__ == "__main__":
    args = sys.argv
    print(args)
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])






