
from weblogo import * 
#write sequences to file that can be read by web logo 

path = "/Users/wyattmccarthy/Desktop/MIT Aging Project/Pre-Processing/positive_tata_sequences_Myotis_lucifugus__little_brown_bat__HLmyoLuc1.fa"

file = open(path)
sequences = read_seq_data(file)
logodata = LogoData.from_seqs(sequences)
logooptions = LogoOptions()
logooptions.title = "A Logo Title"
logoformat = LogoFormat(logodata, logooptions)
eps = eps_formatter(logodata, logoformat)
eps_string = eps.decode()
print(eps_string)
# pdf = pdf_formatter(logodata, logoformat)
    
