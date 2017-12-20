
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
  inner_join(pubs)

pubs_sum <- pubs %>% 
  select(PubMedIds, model, Source, FullJournalName) %>% 
  distinct() %>% 
  group_by(model, Source) %>% 
  summarise(N = n())

pubp_hs <- pubs_sum %>% 
  filter(str_detect(model, "Human")) %>% 
  top_n(25) %>%
  ggplot(aes(reorder(Source, desc(N)), N)) +
  geom_point() +
  labs(y = "Number of publications") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

pubp_other <- pubs_sum %>% 
  filter(!str_detect(model, "Human")) %>% 
  top_n(25) %>% 
  ggplot(aes(reorder(Source, desc(N)), N)) +
  geom_point() +
  labs(y = "Number of publications") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.title.x = element_blank())

pubp_pval <- pubs %>% 
  filter(Accession %in% p_values$Accession) %>% 
  select(PubMedIds, model, Source, FullJournalName) %>% 
  distinct() %>% 
  group_by(model, Source) %>% 
  summarise(N = n()) %>% 
  top_n(25) %>% 
  ggplot(aes(reorder(Source, desc(N)), N)) +
  geom_point() +
  labs(y = "Number of publications") +
  scale_y_continuous(breaks = seq(0, 10, 1), limits = c(0, 10)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

pg <- lapply(list(pubp_other, pubp_hs, pubp_pval), ggplotGrob)
pg <- add_labels(pg, case = panel_label_case)
pga <- arrangeGrob(grobs = pg, ncol = 3)
grid.draw(pga)
