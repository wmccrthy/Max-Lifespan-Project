"""
FOR GUIDE ON USING ENSEMBL PROGRAMMATICALLY: https://rest.ensembl.org/
    - SPECIFICALLY WE WILL WANT TO USE THE FOLLOWING 'GET' COMMAND: https://rest.ensembl.org/documentation/info/sequence_region
        - TAKES AS INPUT QUERIES OF THE FORMAT: GET sequence/region/:species/:region
        - By that, we just need species (need to test if scientific or informal name) and chromosome/gene with coordinates (figure out what naming style is accepted)

FOR GUIDE ON USING UCSC GENOME BROWSER PROGRAMMATICALLY: https://genome.ucsc.edu/goldenPath/help/api.html
    - ex format of API call: https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chrM;start=123;end=456

FOR GUIDE ON HOW TO GET PROPER GENOME INFORMATION FROM TOGA: 
    - https://github.com/hillerlab/TOGA/discussions/145
    
    intial thinking: 
    - map taxon_id to species_common_name (using overview.table.tsv)
    - get taxon_id given species common name from per species files 
        - webscrape genome id from https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon={taxon_id}
        - 
"""
import requests as rekwests, sys, ssl
try:
    _create_unverified_https_context = ssl._create_unverified_context
except AttributeError:
    pass
else:
    ssl._create_default_https_context = _create_unverified_https_context

# TESTING OF API USAGE FOR GENOME BROWSERS 


"""
# ENSEMBL API call formatting and test 

server = "https://rest.ensembl.org"
# ext = "/sequence/region/human/X:1000000..1000100:1?"
ext =  "/sequence/region/callithrix_jacchus/X:1000000..1000100"

r = rekwests.get(server+ext, headers={ "Content-Type" : "text/plain"})
 
if not r.ok:
  r.raise_for_status()
 
print(r.text)

# question for mass/script-querying Ensemble is: 
#   - what name format is accepted as parameter? (scientific name, informal name, ...?)
#   - how do I get proper chromosome identifier for each species? 
"""

# UCSC GENOME BROWSER API call formatting and test
serv = "https://api.genome.ucsc.edu/list/hubGenomes?hubUrl="
specs = "https://hgdownload.soe.ucsc.edu/hubs/primates/hub.txt"

serv2 = "https://api.genome.ucsc.edu/getData/sequence"
specs2 = "?hubUrl=http://hgdownload.soe.ucsc.edu/hubs/mammals/hub.txt;genome=GCF_002288925.2;chrom=MT;start=1;end=10"
# GCA_002288925.3
# GCF_002288925.2



# GCA_011100555.1

specs3 = "?hubUrl=http://hgdownload.soe.ucsc.edu/hubs/primates/hub.txt;genome=GCA_011100555.1;chrom=CM021918;start=1;end=10"

# r = rekwests.get(serv + specs,  headers={ "Content-Type" : "text/plain"})
# if not r.ok: 
#     r.raise_for_status()
#     sys.exit()
# print(r.text);

# r2 = rekwests.get(serv2 + specs2,  headers={ "Content-Type" : "text/plain"})
# if not r2.ok: 
#     r2.raise_for_status()
# print(r2.text)

# r2 = rekwests.get(serv2 + specs3,  headers={ "Content-Type" : "text/plain"})
# if not r2.ok: 
#     r2.raise_for_status()
# print(r2.text)
global sesh_UCSC
sesh_UCSC = rekwests.Session()
sesh_UCSC.headers.update({ "Content-Type" : "text/plain"})

# sesh_NCBI = rekwests.Session()
# sesh_NCBI.headers.update({"api-key": "891efcf2208fbb42730d92ebd88f49daff09"})

def query_genome_browser_test(isPrimate, accession, chrom): 
    global sesh_UCSC
    api_base = "https://api.genome.ucsc.edu/getData/sequence"
    if isPrimate: api_base += f'?hubUrl=http://hgdownload.soe.ucsc.edu/hubs/primates/hub.txt;genome={accession};chrom={chrom}'
    else: api_base += f'?hubUrl=http://hgdownload.soe.ucsc.edu/hubs/mammals/hub.txt;genome={accession};chrom={chrom}'
    
    #print(api_base)
    # req = rekwests.get(api_base, headers={ "Content-Type" : "text/plain"})
    req = sesh_UCSC.get(api_base)    
    if not req.ok or "<html>" in req.text: return -1 
    else: return 0

def query_genome_browser_hub(isPrimate, accession, chrom, start, end):
    global sesh_UCSC
    api_base = "https://api.genome.ucsc.edu/getData/sequence"
    if isPrimate: api_base += f'?hubUrl=http://hgdownload.soe.ucsc.edu/hubs/primates/hub.txt;genome={accession};chrom={chrom};start={start};end={end}'
    else: api_base += f'?hubUrl=http://hgdownload.soe.ucsc.edu/hubs/mammals/hub.txt;genome={accession};chrom={chrom};start={start};end={end}'
    
    #print(api_base)
    # req = rekwests.get(api_base, headers={ "Content-Type" : "text/plain"})
    try: 
        req = sesh_UCSC.get(api_base)
    except Exception as e:
        print(e)
        sesh_UCSC = rekwests.Session()
        sesh_UCSC.headers.update({ "Content-Type" : "text/plain"})
        return query_genome_browser_hub(isPrimate, accession, chrom, start, end)
    
    if not req.ok or len(req.text) > 10000:
        print("error on genome browser (hub) query for: ", accession, chrom, start, end)
        new_accession = query_ncbi_accession(accession, isPrimate, chrom)
        
        # if new_accession returns None, that is most likely indicative of this genome not being available on genome browser (at least easily available)
        # so we have to return null sequence 
        if new_accession == None: 
            return [accession, None]

        return query_genome_browser_hub(isPrimate, new_accession, chrom, start, end)
        # req.raise_for_status()
        # query for NCBI 
    
    else: 
        #print(req.text, "\n")
        return [accession, req.json()['dna']]
    

# https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{accession}/links
def query_ncbi_accession(accession, isPrimate, chrom):
    print("querying NCBI for alt. accession of ", accession)
    api_base = f'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{accession}/links'
    req = rekwests.get(api_base, headers={"api-key": "891efcf2208fbb42730d92ebd88f49daff09"})
    # req = sesh_NCBI.get(api_base)
    if not req.ok:
        req.raise_for_status()
    else: 
        # print(req.text)
        for link in req.json()["assembly_links"]:
            alt_accession = link["accession"]
            if alt_accession != accession: 
                print("found alt accession: ", alt_accession)

                test_validity = None 
                if isPrimate: test_validity = query_genome_browser_test(True, accession, chrom)
                else: test_validity = query_genome_browser_test(False, accession, chrom)
                #test the alternate accession to see if it is accepted by genome browser 

                if test_validity == -1: return None #if not valid, return None 

                else: return alt_accession #if valid, return alt accession

    
        return None 


#gets valid chromosome name that UCSC will accept as queryable (hopefully)
def query_ncbi_chrom(accession, chrom):
    # print("querying NCBI for chromosomes of ", accession)
    api_base = f'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{accession}/sequence_reports?role_filters=assembled-molecule'
    req = rekwests.get(api_base, headers={"api-key": "891efcf2208fbb42730d92ebd88f49daff09"})
    if not req.ok:
        req.raise_for_status()
    else: 
        for chrom_data in req.json()["reports"]:
            if 'genbank_accession' in chrom_data:
                if chrom in chrom_data['genbank_accession']: 
                    print("found alt chrom: ", chrom_data["genbank_accession"],  " for: ", chrom)

                    return chrom_data["genbank_accession"]
        return None 

def query_ncbi_chrom_list(accession, chroms):
    # print("querying NCBI for chromosomes of ", accession)
    mappings = {}
    api_base = f'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{accession}/sequence_reports'
    req = rekwests.get(api_base, headers={"api-key": "891efcf2208fbb42730d92ebd88f49daff09"})
    if not req.ok:
        print(req.status_code, " error qeurying NCBI for chromosome mappings of: ", accession) 
        return None 
        #handle this by appending .1 to all chroms for now 
       
    else: 
        if len(req.json()) < 2: return None 
        for chrom_data in req.json()["reports"]:
                if 'genbank_accession' in chrom_data:
                    if chrom_data['genbank_accession'][:-2] in chroms:
                        toga_chrom = chrom_data['genbank_accession'][:-2]
                        print("found alt chrom: ", chrom_data["genbank_accession"],  " for: ", toga_chrom)
                        mappings[toga_chrom] = chrom_data["genbank_accession"]
                else: print(chrom_data) #for testing clarity 
        return mappings 
    
#GCA_004027535.1


#weird thing sometimes where invalid query just returns the entire genome browser html page which is not flagged as request error 
#to catch this, just check for len(request.text) as the entire html page is like 40,000 chars 
    
# GCA_004027535.1 PVKX01002725.1
print(query_genome_browser_test(False, "GCA_004027535.1", "PVKX01002725.1"))



#GCA_004027535.1 PVKX01003706 31797 31947