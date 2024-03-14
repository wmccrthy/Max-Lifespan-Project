import os, sys, csv, Bio
from Bio.Seq import Seq 


"""
SCRIPT FOR RETRIEVING SEQUENCES OF INTEREST FROM 2BIT ASSEMBLIES WHICH TOGA COMPILED FROM "DNA ZOO CONSORTIUM"

retrieving this data is a priority given that TOGA-given chromosome names are guaranteed to align with those in 2bit assembly files, whereas 
TOGA-given chromosome names are NOT guaranteed to align with those in NCBI-derived fasta assembly files 

"""

accession_dict = {}

# overview_table_local = '/Users/wyattmccarthy/Desktop/MIT Aging Project/Pre-Processing/overview.table.tsv'
overview_table_MIT = '/Mounts/rbg-storage1/users/wmccrthy/overview.table.tsv'

with open(f'{overview_table_MIT}') as file:
    for line in file:
        line = line.strip().split()
        organism_suffix, accession = line[-4], line[-3]
        if organism_suffix == '(bp)': continue #skip header line 
        if organism_suffix == 'Zoo': organism_suffix, accession = line[-6], "DNA Zoo Consortium"
        accession_dict[organism_suffix] = accession

"""
STEP 3: 
    ITERATE THROUGH MAIN TABULAR DATA SET (ORGANISM, GENEID, CHROM, START, END, ...)

    RETRIEVE ACCESSION[ORGANISM_SUFFIX] AND CHROM, START, END
        - ONLY RETRIEVING DATA FOR "DNA ZOO CONSORTIUM" ACCESSIONS 
        - if other accession, continue/skip 
"""

def count_TOGA_points():
    cumulative_set_MIT_path = "/Mounts/rbg-storage1/users/wmccrthy/cumulativeTOGAset.csv"
    cnt = 0
    with open(f'{cumulative_set_MIT_path}') as file:
        for row in file:
            cnt += 1
            if cnt % 1000 == 0: print(cnt, "\n")
    print(cnt)

"""
method for pulling sequence data from 2bit assembly files downloaded from TOGA; takes start and end index parameters, which constrain the method to consider only a portion of the data within 
the cumulative TOGA set 
"""
def retrieve(start_ind, end_ind):

    # PATHS FOR SERVER SCRIPTING / ACTUAL RUN
    twoBitToFa_MIT_path = "/Mounts/rbg-storage1/users/wmccrthy/twoBitToFa"
    twoBit_MIT_assembly_path = "/Mounts/rbg-storage1/users/wmccrthy/2bit_assemblies/"
    twoBit_MIT_output_path = "/Mounts/rbg-storage1/users/wmccrthy/2bit_Sequences/"
    cumulative_set_MIT_path = "/Mounts/rbg-storage1/users/wmccrthy/cumulativeTOGAset.csv"
    output_MIT_path = f'/Mounts/rbg-storage1/users/wmccrthy/dataSubSets/2bitSubSet_{start_ind}_{end_ind}.csv'

    index, data_point = 0, 0 #index tracks where within the cumulative data set the loop is, data_point tracks how many points we have written to the output file

    start_ind, end_ind = int(start_ind), int(end_ind)
    
    with open(f'{output_MIT_path}', 'w') as write_to_file: #open file (to write) at output path
        writer = csv.writer(write_to_file)
        writer.writerow(['organism', 'gene_id', 'orthologType', 'accession', 'chromosome', 'start', 'end', 'direction (+/-)', 'intactness', 'sequence']) #write header row w orthologous sequence data fields 
        #write header row 
        
        with open(f'{cumulative_set_MIT_path}') as file: #open cumulative TOGA data 
            for row in file: #for each data point from cumulative TOGA set 

                if index < start_ind:
                    index += 1
                    continue 
                elif index >= end_ind:
                    break 
                # only consider the data points between given 'start_ind' and 'end_ind' 


                row = row.split(",") #split csv line by comma to get line/row in list format 

                organism_id = row[0]
                organism_suffix = organism_id.split("_")[-1]
                #pull out organism identifying information from current line 

                if organism_suffix == "organism" or organism_suffix not in accession_dict: 
                    index += 1
                    continue #skip first row and certain outliers w completely unique accession sources, these will not be formatted properly by the way we split so just skip if coming across these
                

                accession = accession_dict[organism_suffix] #pull organism's accession from accession_dict 

                if accession != "DNA Zoo Consortium": 
                    index += 1
                    continue #skip over entries this script doesn't pertain to 


                chrom, og_start, og_end, orig_sequence = row[3], row[4], row[5], row[-1] #pull out chromosome, coordinate(s) and sequence data from current line

                direction = row[-3]#pull out direction/orientation data from current line 
                """
                LOOK 1000 BASE PAIRS +/- from start, depending on direction
                """
                start, end = None, None #start/end coordinates of 1000 base pair range from which we want to extract 

                if direction == '+': 
                    start, end = str(max(0, int(og_start)-1000)), og_start #make sure we don't provide negative coordinates 
                if direction == '-': 
                    start, end = og_end, str(int(og_end)+1000)


                #print("Data Point: ", index, " Querying :", organism_id, " Accession: ", accession, " Chromosome: ", chrom, " Coords: ", start, "-", end)

                print("Data Point: ", data_point, " at Index: ", index, "\n")

                #increment tracking variables 
                data_point += 1
                index += 1


                if accession == 'DNA Zoo Consortium': 
                    # CAN USE TwoBitToFa to extract specific sequences from .2bit files; all DNA Zoo accessions have corresponding .2bit files
                    
                    #also possible for coordinate range errors on .2bit files, should be more manageable in those instances 
                    organism_fa_output_path = twoBit_MIT_output_path + f'{organism_suffix}.fa'
                    organism_twoBit_path =  twoBit_MIT_assembly_path + f'{organism_suffix}.2bit'

                    #maybe use twoBitInfo to check chromosome length s.t end coordinate is capped at len(chrom) ? 
                    
                    # command = f'{twoBitToFa_path} {organism_twoBit_url}:{chrom}:{start}-{end} {organism_fa_output_path}'

                    #call twoBitToFa command line tool to extract the sequence at given chrom, between given coordinates, and output data to file at given output path
                    os.system(f'{twoBitToFa_MIT_path} {organism_twoBit_path}:{chrom}:{start}-{end} {organism_fa_output_path}')

                    # next, simply retrieve sequence from given output_path fasta file
                    with open(organism_fa_output_path) as file:
                        extracted_sequence = ""
                        for line in file:
                            if line[0] == '>':
                                continue
                            else: 
                                extracted_sequence += line.strip()

                    #combine sequence provided by TOGA w the 1000 base pair range we've extracted from UCSC
                    #when parsing data:
                    #   if orientation is "+", look at first 1000 base pairs
                    #   if orientation is "-", look at last 1000 base pairs                     
                    combined_sequence = None 
                    if direction == '+': 
                        combined_sequence = extracted_sequence + orig_sequence
                        end = og_end
                    else: 
                        #get reverse complement of retrieved sequence and append to beginning of sequence 
                        sequence_of_interest = str(Seq(sequence_of_interest).reverse_complement())
                        combined_sequence = extracted_sequence + orig_sequence
                        start = og_start

                    #compile all relevant data into list (iterable format), print for debugging / monitoring purposes, and write data to output file 
                    cur_data = row[:3] + [accession, chrom, start, end] + row[6:-1] + [combined_sequence.replace("\n", "")]
                    if extracted_sequence == "": 
                        cur_data[-1] = "Coordinate Error"
                    writer.writerow(cur_data)

if __name__ == "__main__":
    args = sys.argv
    print(args)
    # args[0] = current file
    # args[1] = function name
    # args[2:] = function args : (*unpacked)
    globals()[args[1]](*args[2:])

