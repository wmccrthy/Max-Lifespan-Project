"""
SCRIPT FOR QUICKLY COMBING THRU LIFESPAN DATA TO GET IDEA OF HOW MANY ENTRIES ARE MISSING AND FOR WHAT SPECIES 
    - as of now best bet is to manually find missing lifespan data 
"""
import csv, os, sys, time 
import requests as req
from lifespanScraping import scrape_family
from openai import OpenAI
import google.generativeai as genai


lifespan_data = "/Users/wyattmccarthy/Desktop/MIT Aging Project/Everything/EDA/lifespan_data/lifespan_data.csv"


"""
METHOD USED LOCALLY TO FIND AND RETURN A LIST OF ALL SPECIES WITH UNKNOWN MAX LIFESPANS
"""
def find_missing():
    total_species = 0
    missing_data = 0
    missing = set()
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
    
    # for i in missing: print(i)

    # with open("missing_lifespans.csv", "w") as write_to:
    #     writer = csv.writer(write_to)
    #     writer.writerow(["Species", "Max Lifespan", "Source"])
    #     for species in missing:
    #         writer.writerow([species, 'N/A', 'N/A'])

    return missing 


"""
QUERIES CHATGPT/OPENAI API for max lifespan data for missing species

query format: What is the max lifespan of {species}? Please cite your reference.
"""
gpt_client = OpenAI(api_key="sk-proj-D3X7XWYZ2zzTBNaudWHuT3BlbkFJHmerUcXBaCaJN8dw2xPz")

def query_gpt():
    #for species in missing lifespan: 
    missing = find_missing() #retrieve list of all species w missing lifespan data 
    lifespan_info = []
    for species in missing:
        species = " ".join(species.split("_"))
        req = gpt_client.chat.completions.create(
                model="gpt-3.5-turbo",
                messages=[
                    {"role": "system", "content": "You are a helpful assistant."},
                    {"role": "user", "content": f"What is the max lifespan of {species}? Please cite your reference."}
                ]
                )
        lifespan_info.append(req.choices[0].message)
        print(lifespan_info[-1])
    return lifespan_info


"""
QUERIES GOOGLE GEMINI API FOR MAX LIFESPAN DATA ON MISSING SPECIES 


query format: What is the max lifespan of {species}? Please cite your reference.

TO RUN: python3 missing_lifespan_data.py query_gemini
"""
gemini_api_key = "AIzaSyDbYafelS05UpJY63Q9PT5-tEmxJzwotcA"
genai.configure(api_key=gemini_api_key)
model = genai.GenerativeModel("gemini-pro") 
def query_gemini():
    missing = find_missing() #retrieve list of all species w missing lifespan data 
    # lifespans = {}
    with open("lifespan_data/missing_lifespans_info.csv", "w") as write_to:
        writer = csv.writer(write_to)
        writer.writerow(['species', 'lifespan info'])
        
        for species in missing:
            species = " ".join(species.split("_"))
            resp = model.generate_content(f"What is the max lifespan of {species}? Find and cite academic reference")
            time.sleep(5) #15 request-per-minute limit so wait 5 seconds btwn each to ensure not overdoing it 
            try: 
                writer.writerow(["_".join(species.split(" ")), resp.text])
            except:
                print(f'error on response for {species}', resp)
    return 
"""
APPEARS THAT Gemini API output is providing some false/fugaze citations... perhaps there is a better way to prompt it s.t citations are more valid 
perhaps ask for link as well? so at very least it has to send us somewhere 
"""



"""
WRITES OUT MAX LIFESPANS W ORDER FIELD S.T WE HAVE MEANS OF FINDING 'SIMILAR' ANIMALS W VARYING LIFESPANS 
Essentially spits out same information w added taxonomic info in following format: order | family | genus | species | species max lifespan 

TO RUN: python3 missing_lifespan_data.py get_lifespans_per_order
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






