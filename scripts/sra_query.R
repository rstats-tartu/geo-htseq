
## let's try to query SRA instead
library(entrezquery)
q <- '("Homo sapiens"[Organism] OR "Mus musculus"[Organism]) AND rna seq[Strategy]'
sra <- entrezquery::entrez_docsums(q, db = "sra", retmax = 360000)
save(sra, file = sprintf("data/sra_%s.RData", Sys.Date()))
load("data/sra_2017-09-02.RData")

library(dplyr)
library(purrr)
library(xml2)
library(lubridate)

get_bodyconts <- function(xmlstring) {
  xmlstr <- xml2::read_html(xmlstring)
  body <- xml2::xml_child(xmlstr, "body")
  xml2::xml_contents(body)
}

sra <- mutate(sra, body = map(ExpXml, get_bodyconts))

## Make smaller sample
sra_sample <- sample_frac(sra, 0.01)
sra_sample <- mutate(sra_sample, bioproject = map_chr(body, ~xml_text(.x[[9]])))
sra_sample <- mutate_at(sra_sample, vars(CreateDate, UpdateDate), ymd)

## Filter dates
sra_filt <- filter(sra_sample, CreateDate <= "2017/06/19")
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
