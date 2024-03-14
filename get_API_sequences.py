#EXTENDS UPON extractTOGAdata, using the already created TOGA data set (cumulativeTOGAset.csv)
#avoids chromosome mapping step bc it appears that all TOGA-chroms are valid when '.1' is appended to them (and API-based chrom mapping was expensive)
import time, csv, genome_browsing, re, sys, os
from Bio.Seq import Seq 
from urllib.request import urlopen



# METHOD THAT DISSECTS A TOGA DATA URL, RETURNING A LIST OF THE RELEVANT DIRECTORIES FROM THAT PAGE
# USED TO BUILD/MOCK THE FILE STRUCTURE AS WGET IS NOT DOWNLOADING FROM THE TOP DOWN 
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


"""
HAVING BUILT DATA SET OF orders->organism, we can map each organism to the corresponding accession; (STEP 1)

BC SOME TOGA-PROVIDED CHROMOSOME NAMES ARE NOT ACCEPTED/QUERYABLE ON UCSC, WE NEED TO HAVE DICT WHERE DICT[ORGANISM] => DICT OF ORGANISM'S TOGA-GIVEN CHROMOSOMES
THEN, FOR EACH CHROMOSOME OF AN ORGANISM, WE WANT DICT[ORGANISM][TOGA_CHROMOSOME] = CORRESPONDING NCBI/UCSC QUERYABLE CHROMOSOME 
BUILD THIS DICT IN INITIAL PROCESSING (WHEN INITIALLY GOING THRU TOGA FILES)
(STEP 2; now nullified as it cost too much time and for every invalid TOGA chrom, referred to as chr, UCSC seems to accept chr.1 in query)

THEN, WE CAN USE ORGANISM ACCESSIONS, CHROMOSOME, START/END COORDINATES TO QUERY GENOME BROWSER FOR SEQUENCES OF INTEREST (STEP 3)
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

with open(f'{overview_table_MIT}') as file:
    for line in file: #for each line in file (line's correspond to organism and identifying data)
        line = line.strip().split() #format line into list 
        organism_suffix, accession = line[-4], line[-3] #extract organism identifying suffix and accession 
        if organism_suffix == '(bp)': continue #skip header line 
        if organism_suffix == 'Zoo': organism_suffix, accession = line[-6], "DNA Zoo Consortium" #weird formatting hard-coded fix 
        accession_dict[organism_suffix] = accession #map organism identifying suffix to accession 


"""
STEP 2: (waived)
    FOR EACH ORGANISMS IN CHROM-MAPPING DICT:
        - GET ORGANISM ACCESSION
        FOR EACH CHROM IN CHROM-MAPPING[ORGANISM]:
            - QUERY NCBI FOR ALT. CHROM NAME USING ACCESSION
            - IF FOUND, SET CHROM-MAPPING[ORGANISM][CHROM] -> QEURY CHROM  

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
AUXILLARY STEP HERE: GET ALL ORGANISM IDs from Primates order, s.t we know which UCSC genome hub to query (without having to re-run parts of script that already worked)
"""
primates_set = set(dissect_directory('http://genome.senckenberg.de/download/TOGA/human_hg38_reference/Primates'))


"""
Another auxillary step: create set of organisms with lifespan data so we don't check query for organisms for which there is no lifespan data? 
haven't ended up using this but might add it just for fun
perhaps we will want later? 
"""
lifespan_path = '/Mounts/rbg-storage1/users/wmccrthy/lifespan_data.csv'
lifespan_mappings = {}
with open(f'{lifespan_path}') as file:
    for line in file:
        line = line.strip().split(",")
        species_name = line[0]
        lifespan = line[1]
        if lifespan != "Unknown" and lifespan != "Not Established": lifespan_mappings[species_name] = lifespan


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


"""
method for querying UCSC API for sequence data; takes start and end index parameters, which constrain the method to consider only a portion of the data within 
the cumulative TOGA set 
"""
def retrieve(start_ind, end_ind):
    #add preliminary step of iterating thru temp_API_subset_corrected to create dictionary mapping each species to the set of queries that have already been made for that species 
    #purpose of this is to avoid re-querying for duplicate data on successive script runs, which will save time 
    collected_organisms = {}
    for file in os.listdir("/Mounts/rbg-storage1/users/wmccrthy/dataSubSets/"):
        if file[:4] == "temp": 
            API_data_path = "/Mounts/rbg-storage1/users/wmccrthy/dataSubSets/" + file
            print(file)
            with open(API_data_path) as read_from:
                for line in read_from:
                    line = line.split(",")
                    if line[0] == "organism": continue 
                    if len(line) != 10: continue #skip weird newline data 
                    cur_org = line[0]
                    chrom, start, end, direction = line[4], line[5], line[6], line[-3]
                    if direction == '+': start, end = int(start), int(start) + 999
                    else: start, end = max(0, int(end)-999), int(end)
                    #modify coordinates to match what wud have been queried (since we only query for 1000 bp range)

                    if cur_org not in collected_organisms: 
                        collected_organisms[cur_org] = set((chrom, start, end))
                    else: 
                        collected_organisms[cur_org].add((chrom, start, end))

    print(len(collected_organisms))
    print(sum([len(collected_organisms[org]) for org in collected_organisms]))


    # PATHS FOR SERVER SCRIPTING / ACTUAL RUN
    cumulative_set_MIT_path = "/Mounts/rbg-storage1/users/wmccrthy/cumulativeTOGAset.csv"
    output_path_MIT = f"/Mounts/rbg-storage1/users/wmccrthy/dataSubSets/APIsubSet_{start_ind}_{end_ind}.csv" #output file denoted with start/end index for organization/clarity purposes

    start_ind, end_ind = int(start_ind), int(end_ind)

    not_queryable = set() #set used to keep track of organisms whose assemblies are not available via UCSC API s.t we can skip over their data

    index, data_point = 0,0 #index tracks where within the cumulative data set the loop is, data_point tracks how many points we have written to the output file



    with open(f'{output_path_MIT}', 'w') as write_to_file: #open output file at output path 

        writer = csv.writer(write_to_file)
        writer.writerow(['organism', 'gene_id', 'orthologType', 'accession', 'chromosome', 'start', 'end', 'direction (+/-)', 'intactness', 'sequence']) #write header row w orthologous sequence data fields 
       #write header row 

        with open(f'{cumulative_set_MIT_path}') as file:
            for row in file: #iterate through lines in cumulative TOGA data 
                if index < start_ind:
                    index += 1
                    continue 
                elif index >= end_ind:
                    break 
                #only consider the data points between given 'start_ind' and 'end_ind' 

                row = row.split(",") #split csv line by comma to get line/row in list format 

                organism_id = row[0]
                organism_suffix = organism_id.split("_")[-1]
                #pull out organism identifying information from current line 


                if organism_suffix == "organism" or organism_suffix not in accession_dict: 
                    index += 1
                    continue #skip first row and certain outliers w completely unique accession sources, these will not be formatted properly by the way we split so just skip if coming across these
                    
                accession = accession_dict[organism_suffix] #get current organism's accession from accession_dict 


                if accession == 'DNA Zoo Consortium' or "GIGA" in accession or "http" in accession or len(accession) > 25: 
                    index += 1
                    continue 
                #not getting data for DNA Zoo Consortium based assemblies in this script (also avoids other bugs/unqueryable accessions )


                chrom, og_start, og_end, orig_sequence = row[3], row[4], row[5], row[-1] #pull out chromosome, coordinate(s) and sequence data from current line 


                if "chr" not in chrom and 'v1' not in chrom and "HiC" not in chrom and '.1' not in chrom:  chrom += '.1' 
                #don't reformat chromosomes of standard naming metric (ex. chr1) or those of format incompatible w appended ".1"

                direction = row[-3] #pull out direction/orientation data from current line 
                """
                LOOK 1000 BASE PAIRS +/- from start, depending on direction
                """

                start, end = None, None #start/end coordinates of 1000 base pair range from which we want to extract 
                if direction == '+': 
                    start, end = max(0, int(og_start)-1000), int(og_start)-1 #make sure we don't provide negative coordinates 
                elif direction == '-': 
                    start, end = int(og_end)+1, int(og_end)+1000
                else: continue #no direction indicates bad data so continune 


                if organism_id in collected_organisms: 
                    if (chrom, start, end) in collected_organisms[organism_id] or (chrom, start, end+1) in collected_organisms[organism_id] or (chrom, start-1, end) in collected_organisms[organism_id]: 
                        # print("duplicate data| ", (chrom, start, end), "| don't query again")
                        index += 1
                        continue 
                #dont query for organism data we have already collected
                #includes checks for multiple tuples given that we modified start and end by +1 and -1, respectively, rather than including the start/end coordinate provided by TOGA (to avoid 1 duplicate nucleotide)

                if organism_id in not_queryable: #if organism can't be queried via UCSC's genome browser API, write out relevant data and continue 
                    if direction == '+': end = og_end
                    else: start = og_start
                    curData = row[:3] + [accession, chrom, start, end] + row[6:-1] + ["Not Queryable"]
                    writer.writerow(curData)
                    index += 1
                    data_point += 1
                    continue 

                # print("Querying :", organism_id, " Accession: ", accession, " Chromosome: ", chrom, " Coords: ", start, "-", end)


                # check if in primates order 
                accession_check = None #variable used s.t we can update accession dictionary in instances where the TOGA-given accession is different than what genome browser accepts 
                sequence_of_interest = None 
                if organism_id in primates_set: 
                    accession_check, sequence_of_interest = genome_browsing.query_genome_browser_hub(True, accession, chrom, start, end) #query UCSC primates hub 
                else: 
                    accession_check, sequence_of_interest = genome_browsing.query_genome_browser_hub(False, accession, chrom, start, end) #query UCSC mammals hub 

                if accession_check != None and accession_check != accession: #if found valid, alternate accession for current organism then update accession_dict
                    print("updating accession dict for: ", organism_id, " from ", accession, " to ", accession_check, '\n')
                    accession_dict[organism_suffix] = accession_check
                

                #initialize curData (row to be written to output)
                curData = row[:3] + [accession, chrom, start, end] + row[6:-1] + [None]


                #see if sequence of interest was not retrieved as this indicates two possible cases:
                #   1. genome in question is not queryable via UCSC API browser
                #   2. coordinates of the query were out of bounds (query_end > actual_end)
                if sequence_of_interest == None: 
                    #test if it's coordinate issue (queried coordinates out of range) or genome is actually invalid 
                    #to do this, simply query for general data pertaining to the organism's accession and chrom
                    #if this query returns an error, it indicates the genome is not available on UCSC
                    #otherwise, indicates the genome IS queryable, but the coordinates we queried were NOT 
                    test_validity = None
                    if organism_id in primates_set: test_validity = genome_browsing.query_genome_browser_test(True, accession, chrom)
                    else: test_validity = genome_browsing.query_genome_browser_test(False, accession, chrom)

                    if test_validity == -1: 
                        not_queryable.add(organism_id)
                        curData[-1] = "Not Queryable" 
                        print(organism_id, "not available on genome browser API\n")
                    else: 
                        #modify sequence of interest data to indicate why there is no data (chromosome coordinate error)
                        print("coordinate out of bounds error for: ", accession, chrom, start, end)
                        #these can also be chromosome name errors, given that we waived this mapping part 
                        curData[-1] = "Coordinates Out of Bounds"

                    print("Data Point: ", data_point, " at Index: ", index, "\n")
                    index += 1
                    data_point += 1
                    writer.writerow(curData) 
                    continue 

                #combine sequence provided by TOGA w the 1000 base pair range we've extracted from UCSC
                #when parsing data:
                #   if orientation is "+", look at first 1000 base pairs
                #   if orientation is "-", look at last 1000 base pairs 
                combined_sequence = ""
                if direction == '+': 
                    combined_sequence = sequence_of_interest + orig_sequence
                    end = og_end
                else: 
                    #get reverse complement of retrieved sequence and append to beginning of sequence (new discovery)
                    sequence_of_interest = str(Seq(sequence_of_interest).reverse_complement())
                    # old: combined_sequence = orig_sequence + sequence_of_interest | incorrectly appended sequence 
                    combined_sequence = sequence_of_interest + orig_sequence #new, correct order of retrieved sequence w TOGA-given sequence 
                    start = og_start
                

                #compile all relevant data into list (iterable format), print for debugging / monitoring purposes, and write data to output file 
                curData = row[:3] + [accession, chrom, start, end] + row[6:-1] + [combined_sequence.replace("\n", "").replace('"', '')]
                print("Data Point: ", data_point, " at Index: ", index, "\n")
                writer.writerow(curData)
                
                #increment tracking variables 
                index += 1
                data_point += 1


if __name__ == "__main__":
    args = sys.argv
    print(args)
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])

    



"""
DO NOT NEED THIS ANYMORE (FOR NOW)



STEP 4
OUTPUT SEQUENCES OF INTEREST DATA IN SAME FORMAT AS CUMULATIVE SET [ORGANISM, GENEID, ORTHOLOGTYPE, CHROMOSOME, START, END, DIRECTION, INTACTNESS, SEQUENCE]

with open("sequencesOfInterest.csv", "w") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['organism', 'gene_id', 'orthologType', 'chromosome', 'start', 'end', 'direction (+/-)', 'intactness', 'sequence']) #write header row w sequence data fields 

    # for each species, write it's 
    for row in sequences_of_interest:
        writer.writerow(row)


STEP 5 
OUTPUT separate data sets indicating the organisms / queries for which sequence data could not be found; there are 2 known cases here:
    - genome not on UCSC browser so need to download NCBI assembly and extract sequence 
        - output set of genomes/organisms for which this is the case (organism | accession ) format
    - coordinate error (start/end coordinates were out of bounds)
        - this one is really niche and occurs on a chromosonal basis so output set of queries for which this was the case
        - (organism | accession | chromosome | start | end) format 

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

"""