library(tidyverse)
library(here)

pvalues <- read_rds(here("output/pvalues.rds"))
higher_criticism <- function(x, breaks = 20) {
  x <- na.omit(x)
  counts <- hist(x, breaks = breaks, plot = FALSE)$counts
  b <- 0.05
  qc <- qbinom(1 - b * 0.05, length(x), b)
  counts_over_qc <- counts >= qc
  end_of_peak <- rle(diff(counts_over_qc))$lengths[1] + 1
  !any(counts_over_qc[end_of_peak:length(counts_over_qc)])
}

safe_HC <- safely(higher_criticism)
pvalues_qc <- pvalues %>% 
  select(suppdata_id, features, pvalues) %>% 
  mutate(HC = map(pvalues, safe_HC, breaks = 42))
pvalues_qc_pass <- pvalues_qc %>% 
  mutate(QC_pass = map(HC, "result")) %>% 
  filter(!map_lgl(QC_pass, is.null)) %>% 
  unnest(QC_pass) %>% 
  select(suppdata_id, QC_pass)

megatab <- read_csv(here("output/megatab.csv"))
megatab_qc <- megatab %>% 
  left_join(pvalues_qc_pass) %>% 
  select(suppdata_id, Type, pi0, QC_pass)
megatab_qc %>%
  count(Type, QC_pass)
