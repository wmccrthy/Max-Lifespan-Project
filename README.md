(in progress) Overview of Relevant Files

DIRECTORY: *annotation_scripts*

- Directory containing various scripts for the purpose of finding instances (and psuedo-instances) of designated motifs within DNA sequences
- CONTENTS:
    - data_per_motif.py
        - script with various methods for purpose of iterating through promoter region data and outputting file containing orthologs that contain the designated motif (or all motifs)
    - get_all_motif_annotations.py
        - script with various methods that serve the following purpose:
            - given a motif to look for, iterates thru all promoter region data and identifies instances of given motif; assigns a list of motif instances for each organism in which motif is found and outputs file for each organism containing the instances of given motif and their frequency
    - remaining scripts in this directory either informed the above two or are no longer relevant for use

DIRECTORY: *early_training*

- Directory containing preliminary contents for perceiver model that is meant to learn to predict lifespan given DNA sequence
- CONTENTS:
    - get_train_data.py
        - script that randomly samples an inputted quantity of data from cumulativeTOGAset
    - train.py
        - contains model class that incorporates DNABERT-S embedding and lucidrains' PerceiverIO
        - contains methods to assist with tokenizing data and training model

DIRECTORY: *EDA*

- Directory containing some collected data and various scripts for exploratory data analysis
- CONTENTS:
    - gene_sets_data (directory)
        - directory containing csv files with metadata on our gene datasets.
        - MOST RELEVANT CONTENTS:
            - regulatory_bins_cluster_data_normal_z.csv & regulatory_bins_cluster_data_modified_z.csv
                - csv files containing data collected with get_regulatory_embeddings.py (see script for more details on computation)
                - files contain information on mean lifespan, lifespan iqr, z score (mean based)/modified z score (median based), p value, etc for each cluster of each gene dataset
                - statistical significance (p value) in ...modified_z.csv is computed based off modified z score, while in ...normal_z.csv it is computed based off normal z score.
            - gene_rankings.csv
                - csv file containing the following info for each gene dataset:
                      - gene name, average statistical significance of clusters, most statistically significant cluster, least statistically significant cluster, avg #species/lifespans per cluster
    - lifespan_data (directory)
        - directory containing csv files with data on the lifespans for each TOGA species and how those lifespans vary by taxonomic category.
        - MOST RELEVANT CONTENTS:
              - lifespan_data.csv
                  - csv file with lifespan for every TOGA species; species for which we could not find data marked with "Unknown" lifespan
              - orders/families/genera_of_interest.csv 
                  - csv files containing data on taxonomic orders, families, and genera for which there is a species with outlying lifespan.
                  - for details on how these were computed, see query_data.py
              - missing_lifespan_info.csv
                  - csv file containing lifespan information retrieved by Gemini API (google's LLM). Gemini was queried with "What is the max lifespan of {species}? Cite your reference"
                  - data needs to be validated and has not been incorporated into analysis
                  - see missing_lifespan_data.py for details 
    - embeddings_visuals (directory)
        - folder containing embedding visualizations for a sample of our gene datasets. 
    - get_regulatory_embeddings.py
        - collection of methods for embedding (using DNABERT-S) our sequence data then performing PCA and K-means clustering on embeddings for visualization or statistical analysis.
    - gene_type_analysis.py
        - collection of methods for evaluating sets of orthologs that code for regulatory genes
        - main focus is analyzing sequence length distribution and nucleotide composition similarity across orthologs from the same set
    - get_distribution.py
        - collection of methods for evaluating distribution of sequence lengths across cumulative ortholog set and individual ortholog sets (organized by regulatory gene type)
    - trim_data.py
        - collection of methods for cleaning cumulative and individual ortholog sets
    - query_data.py
        - collection of methods for analyzing and getting basic metrics on tabular data we have already curated (pertains only to csv files in EDA directory)
    - THESE ALL PERTAIN TO TOGA-PROVIDED (CODING SEQUENCE) DATA
    - missing_lifespan_data.py
        - collection of methods for identifying what species we don't have lifespan data for and attempting retrieval via various APIs. 

UNORGANIZED SCRIPTS:

- compile_API_data.py
    - script used to compile all API-retrieved (promoter region) data together s.t duplicates were removed and data was generally cleaned
- correct_orientation.py
    - script used to reorganized API-retrieved sequences with ‘-’ orientation as they were initially configured incorrectly
- count_data.py
    - script used to count quantity of API-retrieved data
- data_per_gene.py
    - script used to organize cumulative TOGA ortholog set into individual sets categorized by the regulatory gene that included orthologs code for
- download_all_assemblies.py
    - script used for downloading NCBI assemblies for organisms whose DNA data is not available at any API (that we’ve found)
- extractTOGAdata.py
    - script used for iterating thru files taken directly from TOGA and compiling cumulative set of orthologs w relevant labels
    - set is formatted: [organism, gene_id, orthologType, chromosome, start, end, direction (+/-), intactness, sequence]
- genome_browsing.py
    - script with methods for querying NCBI and UCSC genome browser APIs
- get_2bit_sequences.py
    - script for using TOGA-provided chromosome, start, end coordinates to extract corresponding promoter region from downloaded 2bit assembly files
- get_API_sequences.py
    - script for using TOGA-provided chromosome, start, end coordinates to query for corresponding promoter region from UCSC genome browser
- get_valid_data.py
    - script for cleaning data (removing newlines, weird characters, etc) that disrupt smooth parsing
- getTOGAfileStructure.py
    - script for recreating the TOGA directory structure (human_hg38_reference → orders directory → organism directory → organism files)  within our local storage s.t we could download all TOGA in organized manner
- lifespanScraping.py (there is also a copy of this in EDA directory)
    - script used for web-scraping lifespan data from various online databases
        - primarily AnAge database and then used AnimalDiversityWeb for species who AnAge did not include
        - script also returns source link s.t validity of data can be traced
