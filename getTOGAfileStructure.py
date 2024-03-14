from urllib.request import urlopen 
import ssl 
import re 
from lifespanScraping import scrape_max_lifespan
import time 
import csv 
import os 


# SOME CHEESE TO AVOID SSL CERTIFICATE ERROR THAT WAS OCCURING WHEN TRYING TO WEB-SCRAPE
try:
    _create_unverified_https_context = ssl._create_unverified_context
except AttributeError:
    pass
else:
    ssl._create_default_https_context = _create_unverified_https_context


# METHOD THAT DISSECTS A TOGA DATA URL, RETURNING A LIST OF THE RELEVANT DIRECTORIES FROM THAT PAGE
# USED TO BUILD/MOCK THE FILE STRUCTURE AS WGET DOESN'T DOWNLOAD FROM THE TOP DOWN 
def dissect_directory(url): 
    test_page =  urlopen(url)
    html_bytes = test_page.read()
    parseable_html = html_bytes.decode('utf-8')
    # get parseable html from passed TOGA url (parseable_html is in string format)

    find_rows = [row.start() for row in re.finditer("<tr>", parseable_html)]
    # TOGA data page lists directories/files as table rows (<tr> in html) in html
    # As such, find_rows is a list of the starting index (in html string) of each of those rows 

    rows = []
    for i in range(len(find_rows)-1): # iterate through the starting index of each row 
        start, end = find_rows[i], find_rows[i+1]  
        sub_row = parseable_html[start:end]#the starting index of current row and starting index of next row describe the range of the current row

        entry_start, entry_end = sub_row.find("<a href="), sub_row.find('/">') #find regexes within the current row that indicate where the directory/file name is 
        entry = sub_row[entry_start+9:entry_end] #entry is set = to the directory/file name of current row 

        rows.append(entry) #add entry to returned list 
    
    
    return rows 

# STEP 1: WEBSCRAPE ALL THE ORDERS FROM HG38 DIRECTORY 
#   - BUILD DICTIONARY OF ORDERS THAT WILL HAVE ORDER:{dict of species} PAIRS 
#   - species dict within orders dict will map species to their TOGA URL 
test_url = "http://genome.senckenberg.de/download/TOGA/human_hg38_reference"
orders = dissect_directory(test_url)
orders = orders[3:-2]
orders = {order:{} for order in orders}

del orders["Matrix"]
del orders["MultipleCodonAlignments"]
# these two are not mammalian orders, have separate data that we do not need 

print(orders)
print(len(orders))


# NOW WE HAVE DICTIONARY OF ORDERS | STRUCTURE WILL FOLLOW S.T WE HAVE ORDER:[LIST OF SPECIES IN THAT ORDER]

# STEP 2: ITERATE THROUGH ALL ORDERS; AS WE GO, CREATE FILE STRUCTURE, WEBSCRAPE ALL SPECIES FROM EACH ORDER AND FILL ORDER DICTIONARY ,
# CAN ALSO CREATE CUMULATIVE LIFESPAN DATA SHEET HERE 
#   - also want to incorporate shell commands for building directory structure here 
#   - species dict within orders dict will map species to their TOGA URL 

# BC WGET IS BEING STUPID, I AM GOING TO WEBSCRAPE THE DIRECTORY STRUCTURE, 
# BY THAT, I WILL FIRST RETRIEVE ALL ORDERS FROM: https://genome.senckenberg.de/download/TOGA/human_hg38_reference/
    # THEN FOR EACH ORDER, I WILL RETRIEVE ALL SPECIES FROM: https://genome.senckenberg.de/download/TOGA/human_hg38_reference/{ORDER}/
        # THEN FOR EACH SPECIES, WE call wget on urls:
        #                   1. https://genome.senckenberg.de/download/TOGA/human_hg38_reference/{order}/{species}/geneAnnotation.bed.gz
        #                   2. https://genome.senckenberg.de/download/TOGA/human_hg38_reference/{order}/{species}/codonAlignment.fa.gz
        #                   3. https://genome.senckenberg.de/download/TOGA/human_hg38_reference/{order}/{species}/orthologsClassification.tsv.gz
        #                   4. https://genome.senckenberg.de/download/TOGA/human_hg38_reference/{order}/{species}/geneAnnotation.gtf.gz 
    

# human_ref_path = "/Users/wyattmccarthy/Desktop/'MIT Aging Project'/fileTest"
MIT_ref_path = "/Mounts/rbg-storage1/users/wmccrthy/human_hg38_reference"
# change this tomorrow when i move to diff cluster 

lifespan_data = {}

# for each order in the built dictionary of TOGA's mammalian orders 
for order in orders:  
    order_url = f'http://genome.senckenberg.de/download/TOGA/human_hg38_reference/{order}' #follow TOGA url structure to get URL of the current order

    species = dissect_directory(order_url)[3:] #retrieve all species directories from the current order's TOGA page 

    # make directory for this order at path: /.../human_hg38_reference/order (following TOGA file structure)
    order_path = MIT_ref_path + f'/{order}'
    os.system(f'mkdir {order_path}')

    # for each species in this order 
    for s in species:
        species_url = order_url + f'/{s}' #follow TOGA url structure of get url of current species 

        orders[order][s] = species_url

        # make directory for species at path: .../human_hg38_reference/order/species (following TOGA file structure)
        species_path = order_path + f'/{s}'
        os.system(f'mkdir {species_path}')

        # download relevant TOGA files to species directory with wget; relevant files are retrieved via species_url and downloaded at species_path 
        if os.system(f'test -f {species_path}/geneAnnotation.bed.gz') == 256: os.system(f'wget --no-check-certificate -P {species_path} {species_url}/geneAnnotation.bed.gz')
        if os.system(f'test -f {species_path}/codonAlignments.fa.gz') == 256: os.system(f'wget --no-check-certificate -P {species_path} {species_url}/codonAlignments.fa.gz')
        if os.system(f'test -f {species_path}/orthologsClassification.tsv.gz') == 256: os.system(f'wget --no-check-certificate -P {species_path} {species_url}/orthologsClassification.tsv.gz')
        if os.system(f'test -f {species_path}/geneAnnotation.gtf.gz') == 256: os.system(f'wget --no-check-certificate -P {species_path} {species_url}/geneAnnotation.gtf.gz')
       
        """
        # AnAge database only needs scientific name which is extracted from TOGA full species identifier 
        species_scientific_name = "_".join(s.split("_")[:2]) #extract species scientific name from full TOGA species id (currently stored in variable 's')
        species_lifespan = scrape_max_lifespan(species_scientific_name) #call method to scrape max lifespan data for current species from AnAge database 
        lifespan_data[species_scientific_name] = species_lifespan #update species entry in lifespan_data dictonary 
        """

        # WAIT S.T WEB-SCRAPING DOESN'T OVERLOAD THE DB'S SERVER W QUERIES (no waiting will do so)
        time.sleep(.001)
    


# print(len(lifespan_data.keys()))
# print(len([species for species in lifespan_data if lifespan_data[species][0] != None and lifespan_data[species][0] != "Not Established"]))
#indicates how many species we have and how many have lifespan data available 


# STEP 3: write lifespan data to sheet (already done, just copy sheet over)
"""
with open("lifespan_data.csv", 'w') as csvfile: 
    writer = csv.writer(csvfile)
    writer.writerow(["Species", "Max Lifespan", "Source"]) #write header row w lifespan data fields 

    # for each species, write it's lifespan data to a row 
    for key in lifespan_data:   
        max_lifespan, source = lifespan_data[key]
        writer.writerow([key, max_lifespan if max_lifespan != None else "Unknown", source if source != None else "N/A"])
"""





# STEP 5: RUN broaderTesting.py (ensure it has correct path), it will parse every single TOGA file 
#   - question is how we want to output data?
#   - maybe hold off on this for now and modify broaderTesting s.t it does the next step we need (look 1000 base pairs from TSS)






    

    



