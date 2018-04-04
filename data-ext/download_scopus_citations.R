
# depends on ds_redline object
if (!"ds_redline" %in% ls()) {
  source("R/results_query.R")
}

# depend on p_values object
if (!"p_values" %in% ls()) {
  source("R/results_pvalues.R")
}

## ---- publications ----

pubs <- readRDS("output/publications.rds")
# jif <- read_csv("data/JIF_incites.csv", skip = 1)

# Drop duplicate Id column
pubs <- pubs %>% select(-Id)

pvals_pub <- p_values %>% 
  select(Accession, suppdata_id, pi0)

pubs <- ds_redline %>% 
  select(Accession, PubMedIds, model) %>% 
  inner_join(pubs) %>% 
  mutate(ISSN = case_when(str_length(ISSN) == 0 ~ ESSN,
                          str_length(ISSN) != 0 ~ ISSN))


# Download citations ------------------------------------------------------

dois <- pubs %>% 
  filter(str_detect(model, "Human")) %>% 
  select(Accession, PubMedIds, DOI) %>% 
  distinct() %>% 
  nest(Accession)

scopus <- crul::HttpClient$new(url = 'https://api.elsevier.com', 
                    headers = list(Accept = "application/json",
                                   Encoding = "UTF-8"))

#' Query function
scopus_n_citations <- function(doi, con) { 
  con$get(path = "content/search/scopus", 
          query = list(query = glue::glue("doi({doi})"),
                       apiKey = Sys.getenv("Elsevier_API")))
}

#' Run query
library(jsonlite)
dois <- dois %>% 
  mutate(resp = map(DOI, scopus_n_citations, con = scopus))
dois <- dois %>% 
  filter(map_lgl(resp, ~ .x$status_code == 200)) %>% 
  mutate(parsed = map(resp, ~ .x$parse()),
         json = map(parsed, jsonlite::fromJSON))

#' Extract citations 
extr_citations <- function(x) {
  citedbycount <- x$`search-results`$entry$`citedby-count`
  if (is.null(citedbycount)) return(NA)
  as.numeric(citedbycount[1])
}

#' Select results for saving
scopcitations <- dois %>% 
  mutate(citations = map_dbl(json, extr_citations)) %>% 
  select(PubMedIds, DOI, citations, json)

#' Save results to rds
readr::write_rds(scopcitations, "data/scopus_citedbycount.rds")
