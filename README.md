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

DIRECTORY: *model*

- Directory containing preliminary contents for transformer model that is meant to learn to predict lifespan given DNA sequence
- CONTENTS:
    - variety of arbitrarily chosen, small sets of training data
    - get_train_data.py
        - script that randomly samples an inputted quantity of data from cumulativeTOGAset
    - train.py
        - contains a class for tokenizing DNA sequences that supports various token lengths
        - contains a transformer model class implementing a preliminary model
        - contains loop for training on various tokenized configurations of DNA sequences

DIRECTORY: *EDA*

- Directory containing some collected data and various scripts for exploratory data analysis
- CONTENTS:
    - various csv files containing statistical data on our gene datasets. 
    - embeddings_visuals (folder)
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
    - THESE ALL PERTAIN TO TOGA-PROVIDED (CODING SEQUENCE) DATA

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
- lifespanScraping.py
    - script used for web-scraping lifespan data from various online databases
        - primarily AnAge database and then used AnimalDiversityWeb for species who AnAge did not include
        - script also returns source link s.t validity of data can be traced
