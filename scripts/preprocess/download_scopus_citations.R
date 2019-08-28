
# Load libs
source("scripts/_common.R")

## ---- publications ----

if (exists("snakemake")) {
  last_date <- snakemake@params[["last_date"]]
  # "output/publications.csv"
  input_pubs <- snakemake@input[["pubs"]]
  # "output/document_summaries.csv"
  input_document_summaries <- snakemake@input[["document_summaries"]]
} else {
  # "output/publications.csv"
  input_pubs <- "output/publications.csv"
  # "output/document_summaries.csv"
  input_document_summaries <- "output/document_summaries.csv"
}

# Publications
pubs <- read_csv(input_pubs, col_types = "cccccccccccccccccccccccccc")
pubs <- rename(pubs, PubMedIds = Id)

# R/A01_GEO_query.R
ds <- read_csv(input_document_summaries, col_types = "ccccccccccccccDccdccccccc") # GEO HT-seq expr datasets


# Date of the first submission
first_date <- range(ymd(ds$PDAT))[1]

# All HT-seq datasets
# Lump together all non-human and murine taxa
# Convert PDAT to date format
ds_redline <- ds %>% 
  mutate(PDAT = ymd(PDAT),
         model = case_when(
           str_detect(taxon, "Mus musculus|Homo sapiens") ~ "Human and mouse",
           !str_detect(taxon, "Mus musculus|Homo sapiens") ~ "Other taxa"
         )) %>%
  filter(PDAT <= last_date)

write_csv(ds_redline, "output/ds_redline.csv")

pubs <- ds_redline %>% 
  select(Accession, PubMedIds, model) %>% 
  inner_join(pubs) %>% 
  mutate(ISSN = case_when(is.na(ISSN) ~ ESSN,
                          TRUE ~ ISSN))

# Download citations ------------------------------------------------------

dois <- pubs %>% 
  select(Accession, PubMedIds, DOI) %>% 
  distinct() %>% 
  nest(Accession)

scopus <- crul::HttpClient$new(url = 'https://api.elsevier.com', 
                    headers = list(Accept = "application/json",
                                   Encoding = "UTF-8"))

#' Query function
scopus_search <- function(doi, con, sleep = 0.06) { 
  Sys.sleep(sleep)
  resp <- con$get(path = "content/search/scopus", 
          query = list(query = glue::glue("doi({doi})"),
                       apiKey = Sys.getenv("Elsevier_geoseq")))
  message("DOI: ", doi, "; success: ", resp$success())
  return(resp)
}

#' Run query
resps <- dois %>% 
  mutate(resp = map(DOI, scopus_search, con = scopus))

#' Extract citations 
library(jsonlite)
extr_citations <- function(x) {
  parsed <- x$parse(encoding = "UTF-8")
  js <- fromJSON(parsed)
  citedbycount <- js$`search-results`$entry$`citedby-count`
  if (is.null(citedbycount)) return(NA)
  as.numeric(citedbycount[1])
}

message("Download finished, starting parsing.\n")
success <- resps %>% 
  filter(map_lgl(resp, ~ .x$success()))

message("Download ok, let's save damn thing.\n")

#' Select results for saving
scopcitations <- success %>% 
  mutate(citations = map_dbl(resp, extr_citations))

scopcitations <- scopcitations %>% 
  select(PubMedIds, DOI, citations, data) %>% 
  unnest() %>% 
  group_by(PubMedIds, DOI, citations) %>% 
  summarise(Accession = str_c(Accession, collapse = ";"))

#' Save results to file
readr::write_csv(scopcitations, "output/scopus_citedbycount.csv")
