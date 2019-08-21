
#' ## Load libs
source("R/_common.R")

## ---- publications ----

pubs <- read_csv(here("output/publications.csv"), 
                 col_types = "cccccccccccccccccccccccccc")

# Rename Id to PubMedIds
pubs <- pubs %>% rename(PubMedIds = Id)

# Import P value stats
pvals_pub <- read_csv("output/pvalues_pool_pub.csv",
                      col_types = "cccdddddd")

# Merge pubs with document summaries (ds_redline)
ds_redline <- read_csv("output/ds_redline.csv",
                       col_types = "ccccccccccccccccccccccccccccc")

pubs <- ds_redline %>% 
  select(Accession, PubMedIds, model) %>% 
  inner_join(pubs) %>% 
  mutate(ISSN = case_when(
    str_length(ISSN) == 0 ~ ESSN,
    str_length(ISSN) != 0 ~ ISSN)
  )

# Tally journals publishing NGS DE experiments
pub_fun <- . %>% 
  select(PubMedIds, Source, ISSN) %>% 
  distinct() %>% 
  group_by(Source, ISSN) %>% 
  summarise(N = n()) %>% 
  ungroup()

# humans
pubsum_hs <- pubs %>% 
  filter(str_detect(model, "Human")) %>%
  pub_fun() %>% 
  top_n(20) %>% 
  arrange(desc(N))

# other taxa
pubsum_other <- pubs %>% 
  filter(!str_detect(model, "Human")) %>% 
  pub_fun() %>% 
  top_n(20) %>% 
  arrange(desc(N))

# Merge publication info with pvalue stats
pubsum_pval <- pubs %>%
  inner_join(pvals_pub) %>% 
  pub_fun() %>% 
  top_n(10) %>% 
  arrange(desc(N))

# Make a plot
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
scopus <- read_csv(here("output/scopus_citedbycount.csv"), 
                   col_types = "ccdc")

#' Distribution of 
p_citatins <- ggplot(scopus) +
  geom_histogram(aes(citations), bins = 30) +
  scale_x_log10()

# Merge publication data with citations and pvalues
publications_citations <- pubs %>%
  inner_join(pvals_pub) %>% 
  left_join(select(scopus, -Accession)) %>%
  select_if(function(x) !is.list(x))
write_csv(publications_citations, "output/publications_citations.csv")

pubs_citations <- publications_citations %>% 
  select(Accession, PubMedIds, DOI, citations, pi0)

#' Correlation between pi0 and number of citations
pi0_citations <- pubs_citations %>% 
  select(-Accession) %>% 
  distinct() %>% 
  filter(!is.na(pi0), !is.na(citations))

p_cit <- pi0_citations %>% 
  ggplot(aes(pi0, citations)) +
  geom_point() +
  geom_smooth(se = FALSE, method = "lm") +
  labs(x = bquote(Proportion~of~true~nulls~(pi*0))) +
  scale_y_log10()

mod_citations <- pi0_citations %>% 
  lm(log10(citations) ~ pi0, data = .) %>% 
  broom::tidy()

p_cit_pval <- pubs_citations %>% 
  select(-Accession) %>% 
  distinct() %>%
  filter(!is.na(citations)) %>% 
  ggplot(aes(citations, fill = is.na(pi0))) +
  geom_histogram(bins = 20, position = "dodge") +
  scale_fill_grey(name = "Anti-conservative\nP values", labels = c("Yes", "No")) +
  labs(y = "N of articles") +
  scale_x_log10()

pg <- lapply(list(p_cit_pval, p_cit), ggplotGrob)
pg <- add_labels(pg, case = panel_label_case)
pga <- arrangeGrob(grobs = pg, ncol = length(pg), widths = c(2, 1))
grid.draw(pga)

## ---- impactfactor ----

#' Importing impact factor data 
jif <- read_csv("data/JIF_incites.csv", skip = 1)
jif <- jif %>% 
  select(-starts_with("X")) %>% 
  rename(FullJournalName = `Full Journal Title`,
         Source = `JCR Abbreviated Title`,
         IF = `Journal Impact Factor`) %>% 
  filter(!is.na(IF)) %>% 
  mutate_at("Source", str_to_lower)

#' Publivatons and p value stats
pp_stats <- pubs %>%
  inner_join(pvals_pub) %>% 
  mutate_at("Source", str_to_lower)

pp_if <- inner_join(pp_stats, jif, by = "Source")
pp_if %>% 
  ggplot() +
  geom_histogram(aes(IF)) +
  facet_wrap(~ Type, scales = "free_y")

