# FOR WEB-SCRAPING FROM THE FOLLOWING RESOURCES: 
#  - https://animaldiversity.org/
#  - https://genomics.senescence.info/species/index.html 

# FIRST ITERATE THROUGH TOGA DIRECTORY AS FOLLOWS:
#   FOR EACH ORDER:
#       FOR EACH SPECIES OF CUR_ORDER:
#           GRAB SPECIES NAME (JUST SCIENTIFIC NAME, SO WILL NEED TO TRIM STRING)
#           SEARCH FOR SPECIES ON WEBSITE 
#           IF RESULT: 
#                   LOCATE AND RECORD LIFESPAN ON RESULT PAGE 
#                   LOCATE AND RECORD SOURCE/CITATION ON RESULT PAGE 
#           WRITE OUT ENTRY TO TABLE OF FORM: SPECIES | MAX_LIFESPAN | SOURCE 

from urllib.request import urlopen
import csv 

"""
METHOD THAT SCRAPES MAX LIFESPAN OF GIVEN SPECIES (SCIENTIFIC NAME) FROM AnAge Database 
"""
def scrape_max_lifespan(species):
    species_url = f'https://genomics.senescence.info/species/entry.php?species={species}'
    # GIVEN SPECIES SCIENTIFIC NAME, FORMAT URL following websites structure 

    species_page =  urlopen(species_url)
    html_bytes = species_page.read()
    parseable_html = html_bytes.decode('utf-8')
    # GET PARSEABLE HTML FROM SPECIES_URL WITH URLLIB.REQUEST LIBRARY 

    # GET LIFESPAN
    start, end = parseable_html.find("Maximum longevity"), parseable_html.find("years")
    lifespan = str(parseable_html[start:end]).strip().replace("\n", " ").replace("</dt>", "").replace("<dd>", "").strip()
    # TRIM THE HTML S.T IT ONLY CONSISTS OF the maximum lifespan 
    
    if len(lifespan) == 0: return [None, None] #functions to avoid error being thrown when species data is unavailable 
    lifespan = lifespan.split()[2]
    if lifespan == "Not": return  ["Not Established", species_url]
    # further trimming html data s.t it is formatted in a way that makes sense 

    # GET SOURCE | this website lists sources s.t you can retrieve source's sub-url and append to their base url 
    base_url = "https://genomics.senescence.info/species/"
    start, end = parseable_html.find("Source"), parseable_html.find("Sample")
    source_sub_section = parseable_html[start:end]
    start, end = source_sub_section.find('="'), source_sub_section.find('">')
    source = base_url + source_sub_section[start+2:end]
    if "<" in source or ">" in source: source = species_url

    return [float(lifespan), source]

"""
Method that scrapes species family and genus from AnAge database
"""
def scrape_family(species):
    species_url = f'https://genomics.senescence.info/species/entry.php?species={species}'
    # GIVEN SPECIES SCIENTIFIC NAME, FORMAT URL following websites structure 

    species_page =  urlopen(species_url)
    html_bytes = species_page.read()
    parseable_html = html_bytes.decode('utf-8')
    # GET PARSEABLE HTML FROM SPECIES_URL WITH URLLIB.REQUEST LIBRARY 

    # GET family
    start, end = parseable_html.find("Family"), parseable_html.find("Genus")
    family = pull_family(parseable_html[start:end])

    #get genus 
    start, end = parseable_html.find("Genus"), parseable_html.find("<dt>Species</dt>")
    genus = pull_family(parseable_html[start:end])

    return family, genus


"""
METHOD TO TRIM THE HTML returned by web-scrape read S.T IT ONLY CONSISTS OF the species' family (find whats between > and </a>)
"""
def pull_family(html_excerpt):
    #find '">" and '</a>' and pull out whats between
    s, e = None, None
    for i in range(len(html_excerpt)):
        if not s:
            if html_excerpt[i:i+2] == '">': s = i+2
        else: 
            if html_excerpt[i:i+4] == '</a>': 
                e = i
                return html_excerpt[s:e]


scrape_family('Callithrix_jacchus')

"""
scrapes Animal Diversity Web for longevity information given a species 
"""
def scrape_adw(species):
    species_url = f"https://animaldiversity.org/accounts/{species}"
    try: species_page =  urlopen(species_url)
    except: 
        print("invalid URL")
        return [None, None]

    html_bytes = species_page.read()
    parseable_html = html_bytes.decode('utf-8')

    start, end = parseable_html.find("Range lifespan"), parseable_html.find("years")
    # print(parseable_html[start:start+100])
    scraped_lifespan = pull_out_lifespan(parseable_html[start:start+100])
    return [scraped_lifespan, species_url]


"""
GIVEN HTML STRING, pulls out highest present integer (corresponding to lifespan)
Dependent on couple cases:
    - if 'months' in string, convert highest int to years 
    - if male/female in string, ...unsure 
"""
def pull_out_lifespan(html_string):
    months = True if 'month' in html_string else False
    html_string = html_string.split(" ")
    for i in range(len(html_string)):
        html_string[i] = strip_non_numeric(html_string[i])
    html_string = [i for i in html_string if len(i) > 0]
    # print(html_string)
    if len(html_string) > 0: 
        lifespan_val = max([float(i) for i in html_string])
        if months: lifespan_val /= 12 
    else: return None 

    return round(lifespan_val, 1)


def strip_non_numeric(string):
    new = ""
    if new.isalpha(): return ''
    
    for c in string: 
        if c.isdigit() or c == '.': new += c
    return new 
        


"""
lifespan_data = "/Users/wyattmccarthy/Desktop/MIT Aging Project/Everything/EDA/lifespan_data.csv"


Updates missing lifespan data via ADW database info 

total_missing = 0
newly_retrieved = 0 
with open("lifespan_data_updated.csv", "w") as write_to:
    writer = csv.writer(write_to)
    writer.writerow(['Species', 'Lifespan', 'Source'])
    with open(lifespan_data) as read_from:
        for line in read_from:
            line = line.split(",")
            if line[1] == "Unknown" or line[1] == "Not Established":
                lifespan, source = scrape_adw(line[0])
                if lifespan == None: lifespan == "Unknown"
                # print(line[0], lifespan, source)
                line = [line[0], lifespan, source]
                if lifespan: newly_retrieved += 1
                total_missing += 1
            writer.writerow(line)
        # print("Recovered", newly_retrieved, "of", total_missing, "Missing Lifespans")



gets rid of stupid formatting in lifespan csv (newlines, quotation marks, etc)

with open(lifespan_data, "w") as write_to:
    writer = csv.writer(write_to)
    with open("/Users/wyattmccarthy/Desktop/MIT Aging Project/Everything/lifespan_data_updated.csv") as read_from:
        for line in read_from:
            line = line.replace('"', "").replace("\n", "")
            line = line.split(',')
            if len(line) < 3: continue 
            if len(line[1]) == 0: line[1] = 'Unknown'
            writer.writerow(line)
            
"""


#SCRIPT BELOW ITERATES THRU LIFESPAN CSV AND ATTEMPTS TO POPULATE 








# SCRAP CODE (USED TO FIGURE OUT SCRAPE_MAX_LIFESPAN METHOD): 
# test_url = "https://genomics.senescence.info/species/entry.php?species=Callithrix_jacchus"

# test_page =  urlopen(test_url)
# html_bytes = test_page.read()
# parseable_html = html_bytes.decode('utf-8')

# # GET LIFESPAN
# start, end = parseable_html.find("Maximum longevity"), parseable_html.find("years")
# lifespan = str(parseable_html[start:end]).strip().replace("\n", " ").replace("</dt>", "").replace("<dd>", "").strip()
# lifespan = lifespan.split()[2]
# # GET SOURCE | this website lists sources s.t you can retrieve source's sub-url and append to their base url 
# base_url = "https://genomics.senescence.info/species/"
# start, end = parseable_html.find("Source"), parseable_html.find("Sample")
# source_sub_section = parseable_html[start:end]
# start, end = source_sub_section.find('="'), source_sub_section.find('">')
# source = base_url + source_sub_section[start+2:end]
# print(lifespan, ",", source)

# FOR EACH ORDER: 
#   FOR EACH SPECIES: 
#           LET SPECIES_NAME = SPECIES SCIENTIFIC NAME 
#           URL TO SCRAPE: "https://genomics.senescence.info/species/entry.php?species={SPECIES_NAME}
#           RUN ABOVE CODE FOR EXTRACTING MAXIMUM LIFESPAN and SOURCE FROM WEBPAGE 
#           WRITE TO OUTPUT FILE 



