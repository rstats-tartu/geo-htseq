library(tidyverse)
library(lubridate)
library(here)
pvalue_types <- read_csv(here::here("output/pvalue_types.csv"), col_types = "cccdcd")
pvalue_types %>% 
  skimr::skim()

srp_stats <- read_csv(here::here("output/srp_stats_extended.csv"), 
                      col_types = "cccdddddcdddddcdddddcddddd")
srp_stats %>% 
  skimr::skim()

ds_redline <- read_csv(here::here("output/ds_redline.csv"),
                       col_types = "cccccccccccccccccccccccccc")
hist_metadata <- ds_redline %>% 
  select(Accession, PDAT, model, taxon, PubMedIds) %>% 
  mutate_at("PDAT", as_date)
hist_type_with_stats <- srp_stats %>% 
  left_join(hist_metadata) %>% 
  arrange(Accession, suppdata_id) %>% 
  select(Accession, PDAT, suppdata_id, Type, pi0, SRP, fp, rs, ud, everything())

taxons <- read_csv(here::here("output/taxons.csv"), col_types = "cccccc")
jif <- read_csv(here::here("output/journal_IF.csv")) %>% 
  mutate_at("FullJournalName", str_to_lower)
publications <- read_csv(here("output/accessions_with_publications.csv"),
                         col_types = "cDcccccccccccccccccccccccddccc") %>% 
  mutate_at("FullJournalName", str_to_lower)
citations <- read_csv(here::here("output/scopus_citedbycount.csv"), col_types = "ccdc") %>% 
  select(Accession, PubMedIds, everything())

geo_publications <- publications %>% 
  left_join(citations) %>% 
  left_join(jif) %>% 
  select(Accession, PDAT, model, taxon, PubMedIds, PubDate, Source, citations, IF) %>% 
  arrange(Accession)


geo_publications %>% 
  skimr::skim()

ylo_megatab <- hist_type_with_stats %>% 
  left_join(geo_publications)

# write output to csv file
write_csv(ylo_megatab, here("output/megatab.csv"))

# get column types string
ylo_megatab %>% 
  map_chr(typeof) %>% 
  map_chr(str_sub, 1, 1) %>% 
  str_c(collapse = "")

# print out skimr summary
ylo_megatab %>% 
  skimr::skim()
