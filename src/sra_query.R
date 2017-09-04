
## let's try to query SRA instead
library(entrezquery)
q <- '("Homo sapiens"[Organism] OR "Mus musculus"[Organism]) AND rna seq[Strategy]'
sra <- entrezquery::entrez_docsums(q, db = "sra", retmax = 360000)
save(sra, file = sprintf("data/sra_%s.RData", Sys.Date()))

library(dplyr)
library(purrr)
library(xml2)
library(lubridate)

# sra$ExpXml[[7]] %>% extract_expxml()
# xmlstring <- sra$ExpXml[[1]]
#
# node <- body_contents[[1]]
# xml_length(body_contents, only_elements = FALSE)
#
# node <- body_contents[[2]]
# node <- body_contents[[8]]
# node <- body_contents[[9]]
#
# lapply(body_contents, function(x) xml_length(x)>0)
#
# nodename <- xml2::xml_name(node)
# nodename
# node <- xml2::xml_contents(node)
# node
#
# fuck <- function(x) {
#   if(xml2::xml_length(x)>0){
#     nodedata <- get_deepnode(x)
#   } else {
#     nodedata <- get_nodedata(x)
#   }
#   nodedata
# }
#
# nodelist <- lapply(node, fuck)
# nodelist
#
# ## Works if data present
# node[[5]] %>% xml_name()
# deepnode_name <- node[[5]] %>% xml_contents() %>% xml_name()
# deepnode_data <- node[[5]] %>% get_deepnode()
# if(length(deepnode_data)==0) {
#   deepnode_data <- NA
# }
# deepnode_attr_names <- names(deepnode_data)
# names(deepnode_data) <- paste0(deepnode_name, '.', deepnode_attr_name)
#

# Link SRA data to GEO ---------------------------------------------------------

get_bodyconts <- function(xmlstring) {
  xmlstr <- xml2::read_html(xmlstring)
  body <- xml2::xml_child(xmlstr, "body")
  xml2::xml_contents(body)
}


sra <- mutate(sra, body = map(ExpXml, get_bodyconts))
sra$body[[2000]] %>%
  .[[9]] %>%
  xml_text()

sra <- mutate(sra, bioproject = map_chr(body, ~xml_text(.x[[9]])))
sra <- mutate_at(sra, vars(CreateDate, UpdateDate), ymd)
sra_filt <- filter(sra, CreateDate <= "2017/06/19")
bp <- unique(sra_filt$bioproject)

i <- sample(1:length(bp), size = 500)

library(httr)
base <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
url <- sprintf("efetch.fcgi?db=%s&id=%s", "bioproject", paste(bp[i], collapse = ","))
res <- POST(file.path(base, url))
xmldoc <- httr::content(res)
d <- xml_contents(xmldoc)
d <- xml_siblings(d)
uid <- d[[2]] %>%
  xml_attrs()
children <- d[[2]] %>%
  xml_children()
children %>%
  xml_contents() %>%
  xml_contents() %>%
  xml_text()

## Parse xml
d <- XML::xmlParse(xmldoc)
d <- XML::xmlRoot(d)
xml_siblings(d)
