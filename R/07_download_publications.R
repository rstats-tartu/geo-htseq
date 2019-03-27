
# Load libs
source("R/_common.R")
library(entrezquery)

# Check if document summary table is already loaded
ds <- read_rds("output/document_summaries.rds")

# All HT-seq datasets
# Lump together all non-human and murine taxa
# Convert PDAT to date format
fix_pmid <- function(x) {
  str_replace_all(x, "(.{8})", "\\1 ") %>% 
    str_trim() %>% 
    str_split(" ") %>% 
    unlist()
}

# Download publications for all taxa
ds_pmids <- ds %>%
  filter(PDAT <= last_date, 
         str_length(PubMedIds) != 0) %>% 
  select(PubMedIds) %>% 
  mutate(PubMedIds = map(PubMedIds, fix_pmid)) %>% 
  unnest() %>% 
  distinct()

publications <- ds_pmids$PubMedIds %>% 
  entrez_docsums(uid = ., db = "pubmed", wait = 0.33)
publications <- publications %>% 
  unnest()
write_csv(publications, "output/publications.csv")
