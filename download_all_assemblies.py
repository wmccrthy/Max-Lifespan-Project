
from urllib.request import urlopen 
import ssl, re, time, os, csv, zipfile 


# SOME CHEESE TO AVOID SSL CERTIFICATE ERROR THAT WAS OCCURING WHEN TRYING TO WEB-SCRAPE
try:
    _create_unverified_https_context = ssl._create_unverified_context
except AttributeError:
    pass
else:
    ssl._create_default_https_context = _create_unverified_https_context



"""
STEP 1
BUILD ACCESSION DICT FOR EACH ORGANISM (DICT[ORGANISM_SUFFIX] = (CORRESPONDING ACCESSION, species scientific name)
"""
accession_dict = {}

overview_table_local = '/Users/wyattmccarthy/Desktop/MIT Aging Project/Pre-Processing/overview.table.tsv'
overview_table_MIT = '/Mounts/rbg-storage1/users/wmccrthy/overview.table.tsv'

with open(f'{overview_table_MIT}') as file:
    for line in file:
        line = line.strip().split()
        organism_suffix, accession = line[-4], line[-3]
        if organism_suffix == '(bp)': continue #skip header line 
        if organism_suffix == 'Zoo': organism_suffix, accession = line[-6], "DNA Zoo Consortium"
        accession_dict[organism_suffix] = (accession, "_".join([line[0], line[1]]))

"""
STEP 2
BUILD SET OF ALREADY-DOWNLOADED ASSEMBLIES SO WE DON'T REDOWNLOAD 

"""
already_downloaded_assemblies_path = "/Mounts/rbg-storage1/users/nmurphy25/zoonomia"
acquired_assemblies = set()

for organism_directory in os.listdir(f'{already_downloaded_assemblies_path}'):
     acquired_assemblies.add(organism_directory)

#print(acquired_assemblies)


twoBit_base_url = "https://genome.senckenberg.de/download/TOGA/MammalianDNAZooAssemblies/"

twoBit_assembly_path_local = "/Users/wyattmccarthy/Desktop/'MIT Aging Project'/Pre-Processing/2bit_assemblies/"
ncbi_assembly_path_local = "/Users/wyattmccarthy/Desktop/'MIT Aging Project'/Pre-Processing/NCBI_assemblies/"

twoBit_assembly_path_MIT = "/Mounts/rbg-storage1/users/wmccrthy/2bit_assemblies"
ncbi_assembly_path_MIT = "/Mounts/rbg-storage1/users/wmccrthy/NCBI_assemblies"


def ncbi_assembly_download(assembly_accession):
    return f'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{assembly_accession}/download?include_annotation_type=GENOME_FASTA'


#DOWNLOAD ASSEMBLIES FOR ALL TOGA GENOME ACCESSIONS 

for organism_suffix in accession_dict:
    #for each organism 

    accession, species = accession_dict[organism_suffix] #retrieve organism's accession  
    print(accession, species)

    #check this before seeing if in Niall's directory as he has some genomes for species which TOGA used .2bit files for 
    if accession == "DNA Zoo Consortium":
        #if from DNA Zoo Consortium, download corresponding 2bit file from TOGA web directory 

        organism_twoBit_url = twoBit_base_url + f'{organism_suffix}.2bit'
        #wget this^ URL into designated 2bit_assemblies directory 

        #print(f'wget --no-check-certificate -P {twoBit_assembly_path_MIT} {organism_twoBit_url}', '\n')d

        #if file we want to download doesn't already exist 
        if os.system(f'test -e {twoBit_assembly_path_MIT}/{organism_suffix}.2bit') != 0: os.system(f'wget --no-check-certificate -P {twoBit_assembly_path_MIT} {organism_twoBit_url}')
        else: print("Already Downloaded 2bit file for ", species, " ", accession)
        continue 


    #don't redownload data on species for which Niall already downloaded assemblies 
    if species.lower() in acquired_assemblies:
        print("Already Acquired Assembly for: ", species, " ", accession)

        continue 

    #try to catch certain outliers w completely unique accession sources 
    if len(accession) <= 5 or len(accession) >= 20: continue 


    #if not DNA Zoo Consortium, download assembly via NCBI 
    assembly_download_url = ncbi_assembly_download(accession)
    #wget this^ URL into designated NCBI assemblies directory 

    #print(f'wget --header="api-key: 891efcf2208fbb42730d92ebd88f49daff09" --no-check-certificate -P {ncbi_assembly_path_MIT} {assembly_download_url}', '\n')
    #if file not already downloaded / in existence 
    if os.system(f'test -e {ncbi_assembly_path_MIT}/{accession}') != 0: os.system(f'wget --header="api-key: 891efcf2208fbb42730d92ebd88f49daff09" --no-check-certificate -P {ncbi_assembly_path_MIT}/{accession} {assembly_download_url}')
    else: print("Already Downloaded Assembly for ", species, " ", accession)
 

"""
TO UNZIP ALL ASSEMBLIES:
    FOR EACH SUB_DIRECTORY IN NCBI_ASSEMBLIES:
        GET PATH TO ZIPPED FILE 
        EXTRACT CONTENTS OF ZIP FILE TO SUB_DIRECTORY 
"""
for accession_dir in os.listdir(ncbi_assembly_path_MIT):
    #all sub-dirs will initially be at: accession_dir/download?include_annotation_type=GENOME_FASTA’
    print(accession_dir)

    #avoid errors w invalid accessions 
    if "_" not in accession_dir: continue 

    sub_dir = f'{ncbi_assembly_path_MIT}/{accession_dir}'
    #bc of error's first check if ncbi_dataset directory exists within zip_path
    if os.system(f'test -e {sub_dir}/ncbi_dataset') == 0: 
        #   if so: indicates assembly file is already unzipped, so remove the zip contents 
        print("Assembly Already Unzipped at: ", sub_dir, " | Removing Unecessary Zip Content")
        os.system(f'rm -rf {ncbi_assembly_path_MIT}/{accession_dir}/download?include_annotation_type=GENOME_FASTA')
    else:
         #   else: assembly still needs to be unzipped (as indicated by ncbi_dataset directory not existing at this path)
        zip_path = f"{ncbi_assembly_path_MIT}/{accession_dir}/download?include_annotation_type=GENOME_FASTA"
        print("Unzipping: ", zip_path)
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(sub_dir)

    

    
