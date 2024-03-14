"""
SCRIPT FOR DEDUCING HOW MANY SPECIES WE WERE ABLE TO GET TATA ANNOTATIONS FOR (GIVEN LIMITATIONS OF API QUERYING OR DATA QUALITY IN GENERAL)

RESULTS: 
    TOTAL ORGANISMS:  536
    ORGANISMS WITH TATA ANNOTATIONS:  390

KEEP IN MIND, SOME ORGANISMS WITH ANNOTATIONS HAD NO INSTANCES OF TATA FOUND 
ALSO STILL MISSING DATA FOR MANY SPECIES 
"""
import os 


hg_38_path = "/Mounts/rbg-storage1/users/wmccrthy/human_hg38_reference/"
annotated_orgs = 0
total_orgs = 0
for order in os.listdir(hg_38_path):
    order_dir = f"{hg_38_path}{order}/"
    for organism in os.listdir(order_dir):

        
        exact_motif_output = f"{order_dir}{organism}/exact_tata_instances.csv"
        similar_motif_output = f"{order_dir}{organism}/similar_tata_instances.csv"

        #if files already exists, continue, no need to redundantly write data  
        if os.system(f'test -e {exact_motif_output}') == 0 and os.system(f'test -e {similar_motif_output}') == 0: 
            annotated_orgs += 1
        total_orgs += 1

print("TOTAL ORGANISMS: ", total_orgs)
print("ORGANISMS WITH TATA ANNOTATIONS: ", annotated_orgs)
            