
library(tidyverse)

# Merge histogram types and pvalue stats data
srp_stats <- read_csv("output/srp_stats.csv")

hist_types_good <- select(hist_types, Accession, suppdata_id = `Supplementary file name`, everything()) %>% 
  filter(Type == "anti-conservative")

srp_stats <- left_join(srp_stats, hist_types_good)


jiffy <- read_csv("data/JIF_incites.csv", skip = 1)
jiffy <- jiffy %>% 
  select(-starts_with("X")) %>% 
  rename(FullJournalName = `Full Journal Title`,
         Source = `JCR Abbreviated Title`,
         IF = `Journal Impact Factor`) %>% 
  filter(!is.na(IF)) %>% 
  mutate_at(c("FullJournalName", "Source"), str_to_lower)

pubs_cit <- read_csv("output/publications_citations.csv", 
                     col_types = "ccccccccccccccccccccccccccccccddddddd")

pubs_cit <- mutate_at(pubs_cit, c("FullJournalName", "Source"), str_to_lower)

pubs_if <- left_join(pubs_cit, jiffy, by = c("FullJournalName", "Source"))
pvalue_impact <- left_join(srp_stats, pubs_if)

write_csv(pvalue_impact, "output/pvalue_impact.csv")
