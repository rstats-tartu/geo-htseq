library(tidyverse)
library(here)
source(here("scripts/histogram_qc.R"))

megatab <- read_csv(here("output/megatab.csv"))
pvalues <- read_rds(here("output/pvalues.rds"))
pvalues_bm <- read_rds(here("output/pvalues_bm.rds")) %>% 
  rename(pvalues_after_filter = pvalues)

safe_QC <- safely(histogram_qc)
breaks = 20

pvalues_qc <- pvalues %>% 
  left_join(pvalues_bm) %>% 
  select(Accession, suppdata_id, starts_with("pvalues")) %>% 
  gather(key, value, starts_with("pvalue")) %>% 
  filter(!map_lgl(value, is.null)) %>% 
  mutate(QC = map(value, safe_QC, breaks = breaks, fdr = 0.05))

pvalues_qc_pass <- pvalues_qc %>% 
  select(Accession, suppdata_id, key, QC) %>%
  mutate(QC = map(QC, "result")) %>% 
  filter(!map_lgl(QC, is.null)) %>% 
  mutate(QC_pass = map_lgl(QC, 1),
         Class = map_chr(QC, 2),
         mids_over_qc = map(QC, 3)) %>% 
  select(Accession, suppdata_id, key, QC_pass, Class, mids_over_qc)

pvalues_qc_pass_plot <- pvalues_qc_pass %>% 
  select(suppdata_id, key, mids_over_qc) %>% 
  filter(map_lgl(mids_over_qc, ~length(.x) >= 1)) %>% 
  unnest(mids_over_qc) %>% 
  ggplot() +
  geom_bar(aes(mids_over_qc)) +
  facet_wrap(~ key)

pvalues_qc_pass_spread <- pvalues_qc_pass %>% 
  select(Accession, suppdata_id, key, Class) %>% 
  mutate_at("key", str_replace, "pvalues", "QCType") %>% 
  spread(key, Class)

megatab_qc <- megatab %>% 
  inner_join(pvalues_qc_pass_spread)

# Update megatab with QC/HC histogram types
megatab_qc %>% 
  write_csv(here("output/megatab.csv"))

megatab_qc %>%
  count(Type, QCType) %>% 
  spread(QCType, n)

pvalues_qc_plots <- pvalues_qc %>% 
  mutate(p = map(value, qc_plot, breaks = breaks)) %>% 
  select(Accession, suppdata_id, key, p)

pvalues_qc_pass %>% 
  select(Accession, suppdata_id, key, Class) %>% 
  left_join(pvalues_qc_plots) %>% 
  mutate(p = map2(p, Class, ~.x + labs(caption = str_c("Assigned type: ", .y))),
         p_saved = pmap(list(suppdata_id, key, p), 
                        ~ggsave(filename = str_c(str_extract(..1, "\\w+"), "_", ..2, ".png"), 
                                plot = ..3, 
                                width = 2, 
                                height = 1.34,
                                path = "output/plots")))
