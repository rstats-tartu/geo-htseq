
# depends on ds_redline object
if (!"ds_redline" %in% ls()) {
  source("R/results_query.R")
}

# depend on p_values object
if (!"p_values" %in% ls()) {
  source("R/results_pvalues.R")
}
  
## ---- publications ----

pubs <- read_rds("output/publications.rds")

# Rename Id to PubMedIds
pubs <- pubs %>% rename(PubMedIds = Id)

pvals_pub <- p_values %>% 
  select(Accession, suppdata_id, pi0)

pubs <- ds_redline %>% 
  select(Accession, PubMedIds, model) %>% 
  inner_join(pubs) %>% 
  mutate(ISSN = case_when(str_length(ISSN) == 0 ~ ESSN,
                          str_length(ISSN) != 0 ~ ISSN))

pub_fun <- . %>% 
  select(PubMedIds, Source, ISSN) %>% 
  distinct() %>% 
  group_by(Source, ISSN) %>% 
  summarise(N = n()) %>% 
  ungroup()

pubsum_hs <- pubs %>% 
  filter(str_detect(model, "Human")) %>%
  pub_fun() %>% 
  top_n(20) %>% 
  arrange(desc(N))

pubsum_other <- pubs %>% 
  filter(!str_detect(model, "Human")) %>% 
  pub_fun() %>% 
  top_n(20) %>% 
  arrange(desc(N))

pubsum_pval <- pubs %>% 
  filter(Accession %in% p_values$Accession) %>% 
  pub_fun() %>% 
  top_n(10) %>% 
  arrange(desc(N))

pubplot <- function(data) {
  ggplot(data, aes(reorder(Source, desc(N)), N)) +
  geom_point() +
  labs(y = "N publications") + 
  scale_color_viridis() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.5),
        axis.title.x = element_blank()) + 
  expand_limits(y = 0)
}

pubp_other <- pubsum_other %>% pubplot()

pubp_hs <- pubplot(pubsum_hs) + theme(axis.title.y = element_blank())

pubp_pval <- pubplot(pubsum_pval) + theme(axis.title.y = element_blank())

pg <- lapply(list(pubp_other, pubp_hs, pubp_pval), ggplotGrob)
pg <- add_labels(pg, case = panel_label_case)
pga <- arrangeGrob(grobs = pg, ncol = 3)
grid.draw(pga)

## ---- citations ----

#' Import citation data
scopus <- read_rds("data/scopus_citedbycount.rds")

# Merge publication data with citations and pvalues
publications_citations <- left_join(pubs, scopus) %>% 
  select_if(function(x) !is.list(x))
write_csv(publications_citations, "output/publications_citations.csv")

pubs_citations <- publications_citations %>% 
  select(Accession, PubMedIds, DOI, citations) %>% 
  left_join(pvals_pub)

#' Correlation between pi0 and number of citations
p_cit <- pubs_citations %>% 
  select(-Accession) %>% 
  distinct() %>% 
  filter(pi0 > 0.25) %>% 
  ggplot(aes(pi0, log10(1 + citations))) +
  geom_point() +
  geom_smooth(se = FALSE, method = "lm") +
  labs(x = bquote(Proportion~of~true~nulls~(pi*0)))

mod_citations <- pubs_citations %>% 
  select(-Accession) %>% 
  distinct() %>% 
  filter(pi0 > 0.25) %>% 
  lm(log10(1 + citations) ~ pi0, data = .) %>% 
  broom::tidy()

p_cit_pval <- pubs_citations %>% 
  select(-Accession) %>% 
  distinct() %>% 
  ggplot(aes(log10(1 + citations), fill = is.na(pi0))) +
  geom_histogram(bins = 60, position = "dodge") +
  scale_fill_grey(name = "P values\navailable", labels = c("Yes", "No")) +
  labs(y = "N of articles")

pg <- lapply(list(p_cit_pval, p_cit), ggplotGrob)
pg <- add_labels(pg, case = panel_label_case)
pga <- arrangeGrob(grobs = pg, ncol = length(pg))
grid.draw(pga)

