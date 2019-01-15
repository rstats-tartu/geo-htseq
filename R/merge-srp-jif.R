
library(tidyverse)

# Merge histogram types and pvalue stats data
hist_types <- read_csv("output/hist_types_after_bm_filter.csv")
srp_stats <- read_csv("output/srp_stats.csv")

hist_types_good <- select(hist_types, Accession, suppdata_id = `Supplementary file name`, everything()) %>% 
  filter(type == 1)

srp_stats <- left_join(srp_stats, hist_types_good)


jiffy <- read_csv("data/JIF_incites.csv", skip = 1)
jiffy <- jiffy %>% select(-starts_with("X"))
colnames(jiffy)[c(2, 3, 6)] <- c("FullJournalName", "Source", "IF")
jiffy <- mutate_at(jiffy, c("FullJournalName", "Source"), str_to_lower)

pubs_cit <- read_csv("output/publications_citations.csv")
pubs_cit <- mutate_at(pubs_cit, c("FullJournalName", "Source"), str_to_lower)

pubs_if <- left_join(pubs_cit, jiffy, by = c("FullJournalName", "Source"))
pvalue_impact <- left_join(srp_stats, pubs_if)


write_csv(pvalue_impact, "output/pvalue_impact.csv")


