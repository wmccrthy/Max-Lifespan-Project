"""
SCRIPT FOR QUICKLY COMBING THRU LIFESPAN DATA TO GET IDEA OF HOW MANY ENTRIES ARE MISSING AND FOR WHAT SPECIES 
    - as of now best bet is to manually find missing lifespan data 
"""
import csv 

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

    with open("missing_lifespans.csv", "w") as write_to:
        writer = csv.writer(write_to)
        writer.writerow(["Species", "Max Lifespan", "Source"])
        for species in missing:
            writer.writerow([species, 'N/A', 'N/A'])

    return 


find_missing()