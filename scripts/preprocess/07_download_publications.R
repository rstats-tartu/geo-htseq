
# Load libs
source("R/_common.R")
library(entrezquery)

# Check if document summary table is already loaded
message("Loading document summaries")
ds <- read_csv(snakemake@input[[1]])

# All HT-seq datasets
# Lump together all non-human and murine taxa
# Convert PDAT to date format
fix_pmid <- function(x) {
  str_replace_all(x, "(.{8})", "\\1 ") %>% 
    str_trim() %>% 
    str_split(" ") %>% 
    unlist()
}
message("Filtering publications by last date")
ds_pmids <- ds %>%
  filter(PDAT <= "2018-12-31") %>% 
  drop_na(PubMedIds) %>% 
  pull(PubMedIds) %>% 
  str_split(pattern = ";") %>% 
  unlist() %>% 
  unique()

message("Downloading publications for all taxa")
publications <- ds_pmids %>% 
  entrez_docsums(uid = ., db = "pubmed", wait = 0.33)
publications <- publications %>% 
  unnest()
write_csv(publications, snakemake@output[[1]])
