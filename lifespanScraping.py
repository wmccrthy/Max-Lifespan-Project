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

# METHOD THAT SCRAPES MAX LIFESPAN OF GIVEN SPECIES (SCIENTIFIC NAME) FROM AnAge Database 
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



