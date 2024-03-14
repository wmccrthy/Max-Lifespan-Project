# THE FOLLOWING CODE ROUGHLY IMPLEMENTS THE EXTRACTION PROCESS FOR A SINGLE ORGANISM 
#   - I WANT TO TRANSLATE THIS SUCH THAT THE PROCESS BELOW IS RAN FOR EACH ORGANISM 
#       FOR MAMMALIAN_ORDER_dir IN HG_38_dir:
#           FOR ORGANISM_dir IN MAMMALIAN_ORDER_dir:
                # WANT 3 FILES: geneAnnotation.bed.gz | codonAlignment.fa.gz | orthologsClassification.tsv.gz | geneAnnotation.gtf.gz 
                # run the below process augmented to include ORGANISM as field of data 


import gzip
import os, sys, genome_browsing, time, csv


# The colors are for intactness:
# - BLUE, "0,0,200": Intact
# - LIGHT_BLUE, "0,200,255": Partially Intact
# - LIGHT_RED, "255,50,50": Lost
# - SALMON, "255,160,120": Uncertain Loss
# - GREY, "130,130,130": Missing
# - BROWN, "159,129,112": Paralogous Projection
intactnessKey = {b'0,0,200':"intact", b'0,200,255':"partially intact", b'255,50,50':'lost', b'255,160,120':'uncertain loss', b'130,130,130':'missing', b'159,129,122':'paralogous'}


cumulativeSet = [] #used to store TOGA orthologous gene/sequence data (has entry for each ortholog)

geneID_to_Index = {} #maps (geneID, organism) to index in cumulative set 
index = 0 #index used in ^ 

orders = {} #dict of orders with species:geneDict pairs where geneDict holds gene:geneAnnotations pairs where geneComponents is an array consisting of [componentType, coordinates]
# fill species:geneDict pairs 
# orders[order] => dictionary mapping species to sub-dictionaries
# orders[order][species] -> dictionary mapping orthologous genes from that species to their components and orthologous types from that species to their frequencies 
# orders[order][species][gene] -> list of annotated components of that gene (form: [component_type, start, end])
# orders[order][species][orthologType] -> # of orthologs from that species of given ortholog type 


ortholog_type_frequency = {b'one2one':0, b'one2many':0, b'many2one':0, b'many2many':0}
#dict for counting frequncies of 

chromosome_mapping = {}
#dict used for mapping organisms TOGA-given chromosome names to UCSC-accepted chromosome names 

# The base path(s) upon which to start parsing TOGA files | change depending on if testing locally or running on servers 
MIT_human_ref_path = "/Mounts/rbg-storage1/users/wmccrthy/human_hg38_reference" 
human_reference_path = "/Users/wyattmccarthy/Desktop/MIT Aging Project/Pre-Processing/human_hg38_reference"


# want format of [organism, gene_id, orthologType, chromosome, start, end, direction (+/-), intactness, sequence]
# organism is gotten from directory, gene_id first from orthologClassification, type from oC, chromosome, start, end, direction, and intactness from geneAnnotation and sequence from codonAlignments 

# FOR each order directory in human_hg38_reference directory 
# Make sure to CHANGE THIS TO MIT_REFERENCE_PATH WHEN TRANSFERRING FILE TO CLUSTERS  


for mammalian_order_dir in os.listdir(human_reference_path):
    if mammalian_order_dir == ".DS_Store": continue #bug fix | os.listdir always includes this file, just skip over

    print(f'{mammalian_order_dir}:') #output the current order being parsed for clarity 

    orders[mammalian_order_dir] = {}#add order to orders dictionary | orders[order] stores a dictionary that will store organism:{} pairs 

    # for each organism/species in order directory 
    for organism_dir in os.listdir(human_reference_path + f'/{mammalian_order_dir}'):
        if organism_dir == ".DS_Store": continue #bug fix | os.listdir always includes this file, just skip over

        orders[mammalian_order_dir][organism_dir] = {} #initialize species dictionary where dict[species][orthologyType] = freq of that orthology type in this species and dict[species][geneID] = [list of all parts (exons, CDS, start codon, ...) of that gene]

        if organism_dir not in chromosome_mapping: chromosome_mapping[organism_dir] = {} #initialize organism's chromosome mapping dictionary that maps TOGA chrom name -> NCBI/UCSC chrom name

        interior_dir_path = human_reference_path + f'/{mammalian_order_dir}/{organism_dir}' #interior_dir_path stores the path to this species (ex. human_hg38_reference/order/species/)
        
        print(organism_dir) #output the current organism being parsed for clarity 

        # iterate through orthologsClassification for current organism 
        with gzip.open(f'{interior_dir_path}/orthologsClassification.tsv.gz') as file:
            for line in file: 
                line = line.strip().split() #format line to get ride of blank spaces and make parseable as list 

                if line[0] == b't_gene':continue #indicates it is the header row (no data)

                geneID, orthologyType = line[3:5] #get columns 3 and 4 (0-indexed) from current line which store geneID and orthologType, respectively 

                if orthologyType == b'one2zero': continue #bad data 

                ortholog_type_frequency[orthologyType] += 1 #update cumulative frequency count for ortholog type 


                if orthologyType in orders[mammalian_order_dir][organism_dir]: orders[mammalian_order_dir][organism_dir][orthologyType] += 1
                else: orders[mammalian_order_dir][organism_dir][orthologyType] = 1 
                #update organism-specific frequency count for ortholog type 


                cumulativeSet.append([organism_dir, geneID, orthologyType, None, None, None, None, None, None]) #add preliminary data to cumulative set
                geneID_to_Index[(geneID, organism_dir)] = index #update index of entry for (geneID, organism) s.t we can access it later and supplement with associated data from geneAnnotation.bed 
                
                #orders[mammalian_order_dir][organism_dir][geneID] = []

                index += 1

        with gzip.open(f'{interior_dir_path}/geneAnnotation.bed.gz') as file:
            for line in file:
                line = line.strip().split() # line now of format: chromsome | start | end | gene_id | irrelevant | direction_to_read 
                chrom, start, end, geneID= line[:4] #pull out chromosome, start and end coordinates 
                direction_to_read = line[5] #pull out orientation of sequence 
                colorVal = line[-4] #pull out colorVal (used to encode intactness)


                if colorVal in intactnessKey: intactness = intactnessKey[colorVal] #get corresponding intactness value from intactness key dict
                else: intactness = None 


                if (geneID, organism_dir) in geneID_to_Index: 
                    cumulativeSet[geneID_to_Index[(geneID, organism_dir)]][3:8] = chrom, start, end, direction_to_read, intactness
                    #update cumulative set entry corresponding to this (geneID and organism)

                if chrom.decode() not in chromosome_mapping[organism_dir]: chromosome_mapping[organism_dir][chrom.decode()] = chrom.decode()
                #update chromosome mappings s.t this organisms TOGA-given chromosome initially maps to itself (used s.t it is not None bc of bugs caused at later point in script)
                
        """
        THIS GETS SPECIFIC PARTS OF EACH GENE (EXON, CDS, START/STOP CODON)
            -for now doesn't seem like we need as these are encompassed by larger sequence 

        with gzip.open(f'{interior_dir_path}/geneAnnotation.gtf.gz') as file:
            for line in file:
                line = line.strip().split()
                region_type, start, end, geneID = line[2], line[3], line[4], line[9][1:-2]
                if geneID in orders[mammalian_order_dir][organism_dir]:
                    orders[mammalian_order_dir][organism_dir][geneID].append([region_type, start, end])
        """ 

        #parse codonAlignments file (get actual sequences)
        with gzip.open(f'{interior_dir_path}/codonAlignments.fa.gz') as file:

            prevIdentifier = None #fasta format uses '>' to denote identifying lines | prevIdentifer will be used to check if found sequence pertains to reference or query genome 
            unAccounted = 0

            for line in file:
                # ord() gets ascii value of passed char; bc gzip reads in bytes we need to check for symbolic chars like such 
                if line[0] == ord('>'):
                    line = line.split(b'|') # formatting line s.t there are no blank spaces | line now of format [geneID, codon, query/reference]
                    for i in range(len(line)): line[i] = line[i].strip() 

                    line[0] = line[0][1:] #remove '>' from geneID of identifier line 
                    prevIdentifier = line #set prevIdentifier s.t we can determine if sequences found thereafter pertain to reference or query 
                else: 
                    # line containing DNA sequence 
                    # we only care abt sequences corresponding to QUERY genome 
                    if prevIdentifier[2] == b'QUERY':

                        geneID = prevIdentifier[0]

                        if (geneID, organism_dir) in geneID_to_Index:
                            cumulativeSet[geneID_to_Index[(geneID, organism_dir)]][-1] = line.strip() #update cumulative set entry with actual DNA sequence 

                        else:unAccounted += 1 # why are some genes from codonAlignments unaccounted for....

    print()

"""
print("Records for Sample Set of 5 Organisms") 
print("Total Data Points: ", len(cumulativeSet))
print("Intact Data Points: ", len([i for i in cumulativeSet if i[-2] == 'intact']))
print("1:1 Data Points: ", len([i for i in cumulativeSet if i[2] == b'one2one']))
print("Intact and 1:1 Data Points: ", len([i for i in cumulativeSet if i[-2] == 'intact' and i[2] == b'one2one']), '\n')
"""

"""
OUTPUT CUMULATIVE SET TO CSV (ensure we do this so we don't have to rerun parts of script): 
"""
with open("cumulativeTOGAset.csv", 'w') as csvfile: 
    writer = csv.writer(csvfile)
    writer.writerow(['organism', 'gene_id', 'orthologType', 'chromosome', 'start', 'end', 'direction (+/-)', 'intactness', 'sequence']) #write header row w orthologous sequence data fields 

    # for each orthologous sequence, ouput its data to csv row 
    for row in cumulativeSet:
        for i in range(len(row)): 
            if type(row[i]) != str and row[i] != None: row[i] = row[i].decode()
        writer.writerow(row)


"""
OUTPUT TOTAL ORTHOLOG FREQUENCY DATA (ensure we do this so we don't have to rerun parts of script):
"""
with open("orthologTypeFrequency.csv", "w") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["category", "1:1", "1:many", "many:1", "many:many"])
    writer.writerow(["total", ortholog_type_frequency[b'one2one'], ortholog_type_frequency[b'one2many'], ortholog_type_frequency[b'many2one'], ortholog_type_frequency[b'many2many']])
    for order in orders:
        one_to_one_cnt, one_to_many_cnt, many_to_one_cnt, many_to_many_cnt = 0,0,0,0
        for organism in orders[order]:
            writer.writerow([f'{organism}', orders[order][organism][b'one2one'], orders[order][organism][b'one2many'], orders[order][organism][b'many2one'], orders[order][organism][b'many2many']])
            one_to_one_cnt += orders[order][organism][b'one2one']
            one_to_many_cnt += orders[order][organism][b'one2many']
            many_to_one_cnt += orders[order][organism][b'many2one']
            many_to_many_cnt += orders[order][organism][b'many2many']
        writer.writerow([f'{order}', one_to_one_cnt, one_to_many_cnt, many_to_one_cnt, many_to_many_cnt])
"""
HAVING BUILT DATA SET OF orders->organism, we can map each organism to the corresponding accession; (STEP 1)

BC SOME TOGA-PROVIDED CHROMOSOME NAMES ARE NOT ACCEPTED/QUERYABLE ON UCSC, WE NEED TO HAVE DICT WHERE DICT[ORGANISM] => DICT OF ORGANISM'S TOGA-GIVEN CHROMOSOMES
THEN, FOR EACH CHROMOSOME OF AN ORGANISM, WE WANT DICT[ORGANISM][TOGA_CHROMOSOME] = CORRESPONDING NCBI/UCSC QUERYABLE CHROMOSOME 
BUILD THIS DICT IN INITIAL PROCESSING (WHEN INITIALLY GOING THRU TOGA FILES)
(STEP 2)

THEN, WE CAN USE ORGANISM ACCESSIONS, NCBI-MAPPED-CHROMOSOME, START/END COORDINATES TO QUERY GENOME BROWSER FOR SEQUENCES OF INTEREST (STEP 3)
"""

"""
STEP 1: 
    ITERATE THROUGH OVERVIEW.TABLE.TSV 
    EXTRACT ACCESSION AND ORGANISM SUFFIX (COLUMNS -3 AND -4)
    MAP ORGANISM_SUFFIX TO THEIR ACCESSION 
"""

accession_dict = {}

overview_table_local = '/Users/wyattmccarthy/Desktop/MIT Aging Project/Pre-Processing/overview.table.tsv'
overview_table_MIT = '/Mounts/rbg-storage1/users/wmccrthy/overview.table.tsv'

with open(f'{overview_table_local}') as file:
    for line in file:
        line = line.strip().split()
        organism_suffix, accession = line[-4], line[-3]
        if organism_suffix == '(bp)': continue #skip header line 
        if organism_suffix == 'Zoo': organism_suffix, accession = line[-6], "DNA Zoo Consortium"
        accession_dict[organism_suffix] = accession


"""
STEP 2:
    FOR EACH ORGANISMS IN CHROM-MAPPING DICT:
        - GET ORGANISM ACCESSION
        FOR EACH CHROM IN CHROM-MAPPING[ORGANISM]:
            - QUERY NCBI FOR ALT. CHROM NAME USING ACCESSION
            - IF FOUND, SET CHROM-MAPPING[ORGANISM][CHROM] -> QEURY CHROM  
"""

"""
AS OF NOW IT APPEARS THAT EVERY SINGLE VALID ALT. CHROMOSOME NAME IS JUST ORIGINAL_CHROM + '.1'
SO FOR NOW WE WILL WAIVE THIS STAGE AS IT IS TIME CONSUMING AND ERROR PRONE 

for organism in chromosome_mapping:
    organism_suffix = organism.split("_")[-1]
    organism_accession = accession_dict[organism_suffix]
    if organism_accession == "DNA Zoo Consortium": continue 
    new_chrom_mappings = genome_browsing.query_ncbi_chrom_list(organism_accession, set(chromosome_mapping[organism].keys()))

    if new_chrom_mappings != None: 
        for chrom in chromosome_mapping[organism]:
            #if toga-given chrom has new mapping, update entry 
            if chrom in new_chrom_mappings: chromosome_mapping[organism][chrom] = new_chrom_mappings[chrom]
    else:
        #handle errors by appending .1 to all chroms for safety purposes (in 5 sample species, all toga-given-chroms map to toga-given-chrom.1 so its a dece bet)
        for chrom in chromosome_mapping[organism]: chromosome_mapping[organism][chrom] = f'{chrom}.1';

    time.sleep(.1)
"""



"""
STEP 3: 
    ITERATE THROUGH MAIN TABULAR DATA SET (ORGANISM, GENEID, CHROM, START, END, ...)
    RETRIEVE ACCESSION[ORGANISM_SUFFIX] AND CHROM, START, END

    QUERY APPROPRIATE GENOME BROWSER HUB (PRIMATE OR MAMMAL) W RETRIEVED ACCESSION, CHROM, AND COORDINATES 
    IF FAILED: 
        - QUERY NCBI API FOR ALT. ACCESSION NAMES (https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{accession}/links)
        - RETRY GENOME BROWSER QUERY
    ELSE: 
        - STORE RETRIEVED SEQUENCE (need to figure out nuances of what we actually want here)
"""
twoBit_base_url = "https://genome.senckenberg.de/download/TOGA/MammalianDNAZooAssemblies/"

#PATHS FOR LOCAL SCRIPTING / MOCK TESTING 
twoBitToFa_path = "/Users/wyattmccarthy/Desktop/twoBitToFa"
twoBit_output_path = "/Users/wyattmccarthy/Desktop/'MIT Aging Project'/Pre-Processing/2bit_Sequences/"
cumulative_set_path = "/Users/wyattmccarthy/Desktop/MIT Aging Project/Pre-Processing/cumulativeTOGAset.csv"

# PATHS FOR SERVER SCRIPTING / ACTUAL RUN
twoBitToFa_MIT_path = "/Mounts/rbg-storage1/users/wmccrthy/twoBitToFa"
twoBit_MIT_output_path = "/Mounts/rbg-storage1/users/wmccrthy/2bit_Sequences/"
cumulative_set_MIT_path = "/Mounts/rbg-storage1/users/wmccrthy/cumulativeTOGAset.csv"

sequences_of_interest = []

not_queryable = set()
coordinate_errors = []

with open(f'{cumulative_set_path}') as file:
  for row in file:
        row = row.split(",")
        organism_id = row[0]
        organism_suffix = organism_id.split("_")[-1]
        if organism_suffix == "organism": continue #skip first row 
        accession = accession_dict[organism_suffix]
        chrom, start, end = row[3], row[4], row[5]

        if accession != "DNA Zoo Consortium":  chrom += '.1' #avoid using alt. chrom names if from DNA Zoo Consortium bc those chromosomes all have proper original names

        #chrom = chromosome_mapping[organism_id][chrom] #map toga-given chrom name to queryable chrom name (this might have been a complete waste, try just mapping to chrom.1)

        if organism_id in not_queryable: #if organism who can't be queried via UCSC's genome browser API, then continue 
            sequences_of_interest.append(row[:3] + [chrom, start, end] + row[6:-1] + [None])
            continue 

        direction = row[-3]
        """
        LOOK 1000 BASE PAIRS +/- from start, depending on direction
        """
        if direction == '+': 
            start, end = str(max(0, int(start)-1000)), start #make sure we don't provide negative coordinates 
        if direction == '-': 
            start, end = end, str(int(end)+1000)

        print("Querying :", organism_id, " Accession: ", accession, " Chromosome: ", chrom, " Coords: ", start, "-", end)


        if accession == 'DNA Zoo Consortium': 
            # CAN USE TwoBitToFa to extract specific sequences from .2bit files; all DNA Zoo accessions have corresponding .2bit files
            
            organism_fa_output_path = twoBit_output_path + f'{organism_suffix}.fa'
            organism_twoBit_url = twoBit_base_url + f'{organism_suffix}.2bit'

            # command = f'{twoBitToFa_path} {organism_twoBit_url}:{chrom}:{start}-{end} {organism_fa_output_path}'
            # print(command)

            os.system(f'{twoBitToFa_path} {organism_twoBit_url}:{chrom}:{start}-{end} {organism_fa_output_path}')

            # next, simply retrieve sequence from {desired_output_path} fasta file (removed apostrophes used in file to allow for system command)
            # these accomodate for spaces in the path names 
            with open(organism_fa_output_path.replace("'", "")) as file:
                extracted_sequence = ""
                for line in file:
                    if line[0] == '>':
                        continue
                    else: 
                        extracted_sequence += line.strip()
                print(organism_id, chrom, start, end, " sequence: ", extracted_sequence, "\n")
            continue 

        # check if in primates order 
        accession_check = None #variable used s.t we can update accession dictionary in instances where the TOGA-given accession is different than what genome browser accepts 
        sequence_of_interest = None 
        if organism_id in orders["Primates"]: 
            accession_check, sequence_of_interest = genome_browsing.query_genome_browser_hub(True, accession, chrom, start, end) #query UCSC primates hub 
        else: 
            accession_check, sequence_of_interest = genome_browsing.query_genome_browser_hub(False, accession, chrom, start, end) #query UCSC mammals hub 
    
        if accession_check != None and accession_check != accession: 
            print("updating accession dict for: ", organism_id, " from ", accession, " to ", accession_check, '\n')
            accession_dict[organism_suffix] = accession_check
        

        sequences_of_interest.append(row[:3] + [chrom, start, end] + row[6:-1] + [sequence_of_interest])

        #see if sequence of interest was not retrieved as this indicates the genome in question is not available via genome browser API and as such, no need to redundantly
        #try and retrieve their data later on 
        if sequence_of_interest == None: 
            #test if it's coordinate issue (queried coordinates out of range) or genome is actually invalid 
            test_validity = None
            if organism_id in orders["Primates"]: test_validity = genome_browsing.query_genome_browser_test(True, accession, chrom)
            else: test_validity = genome_browsing.query_genome_browser_test(False, accession, chrom)

            if test_validity == -1: 
                not_queryable.add(organism_id)
                sequences_of_interest[-1][-1] = "Not Queryable"
                print(organism_id, " not available on genome browser API\n")
            else: 
                #modify sequence of interest data to indicate why there is no data (chromosome coordinate error)
                sequences_of_interest[-1][-1] = "Coordinate Error"
                coordinate_errors.append([organism_id, accession, chrom, start, end])
        

        time.sleep(.001)

print(sequences_of_interest[0])
print(sequences_of_interest[-1])

"""
STEP 4
OUTPUT SEQUENCES OF INTEREST DATA IN SAME FORMAT AS CUMULATIVE SET [ORGANISM, GENEID, ORTHOLOGTYPE, CHROMOSOME, START, END, DIRECTION, INTACTNESS, SEQUENCE]
"""
with open("sequencesOfInterest.csv", "w") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['organism', 'gene_id', 'orthologType', 'chromosome', 'start', 'end', 'direction (+/-)', 'intactness', 'sequence']) #write header row w sequence data fields 

    # for each species, write it's 
    for row in sequences_of_interest:
        writer.writerow(row)

"""
STEP 5 
OUTPUT separate data sets indicating the organisms / queries for which sequence data could not be found; there are 2 known cases here:
    - genome not on UCSC browser so need to download NCBI assembly and extract sequence 
        - output set of genomes/organisms for which this is the case (organism | accession ) format
    - coordinate error (start/end coordinates were out of bounds)
        - this one is really niche and occurs on a chromosonal basis so output set of queries for which this was the case
        - (organism | accession | chromosome | start | end) format 
"""
with open("notQueryable.csv", "w") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["organism", "accession"])
    for organism in not_queryable:
        organism_suffix = organism.split("_")[-1]
        writer.writerow([organism, accession_dict[organism_suffix]])

print(coordinate_errors)
with open("coordinateError.csv", "w") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["organism", "accession", "chromosome", "start", "end"])
    for query in coordinate_errors: 
        writer.writerow(query)



#








