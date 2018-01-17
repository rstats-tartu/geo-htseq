
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

jiffy <- read_csv("data/JIF_incites.csv", skip = 1)
jiffy <- jiffy %>% select(-starts_with("X"))
colnames(jiffy)[c(2, 3, 6)] <- c("FullJournalName", "Source", "IF")
jiffy <- jiffy %>% 
  select(c(2:4, 6)) %>% 
  mutate_at(c("FullJournalName", "Source"), str_to_lower)

pub_fun <- . %>% 
  select(PubMedIds, model, Source, ISSN) %>% 
  distinct() %>% 
  group_by(model, Source, ISSN) %>% 
  summarise(N = n()) %>% 
  ungroup()

pubs_sum <- pubs %>% pub_fun()

pubsum_hs <- pubs_sum %>% 
  filter(str_detect(model, "Human")) %>%
  top_n(20)

pubsum_other <- pubs_sum %>% 
  filter(!str_detect(model, "Human")) %>% 
  top_n(20)

pubsum_pval <- pubs %>% 
  filter(Accession %in% p_values$Accession) %>% 
  pub_fun() %>% 
  filter(N >= 2)

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

pubp_hs <- pubsum_hs %>% pubplot() + theme(axis.title.y = element_blank())

pubp_pval <- pubsum_pval %>% pubplot() + theme(axis.title.y = element_blank())

pg <- lapply(list(pubp_other, pubp_hs, pubp_pval), ggplotGrob)
pg <- add_labels(pg, case = panel_label_case)
pga <- arrangeGrob(grobs = pg, ncol = 3)
grid.draw(pga)
