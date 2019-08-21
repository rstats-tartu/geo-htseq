
if (!any(str_detect(ls(), "pubs$"))) {
  source("scripts/results_publications.R")
}

set.seed(122)
tocheck <- pubs %>%
  filter(str_detect(model, "Human")) %>%
  sample_frac(0.02)

tocheckids <- tocheck %>% 
  group_by(ISSN, Source) %>% 
  nest() %>% 
  arrange(desc(map_int(data, nrow))) %>% 
  unnest() %>% 
  select(ISSN, Source, DOI, Volume, Issue, Pages) %>% 
  distinct()

readr::write_excel_csv(tocheckids, "output/check_pubs_for_pvalues.csv")
