library(tidyverse)
library(here)

megatab <- read_csv(here("output/megatab.csv"))
pvalues <- read_rds(here("output/pvalues.rds"))
higher_criticism <- function(x, breaks, fdr) {
  x <- na.omit(x)
  counts <- hist(x, breaks = seq(0, 1, length.out = breaks), plot = FALSE)$counts
  b <- 1 / breaks
  qc <- qbinom(1 - b * fdr, length(x), b)
  counts_over_qc <- counts >= qc
  end_of_peak <- rle(diff(counts_over_qc))$lengths[1] + 1
  !any(counts_over_qc[end_of_peak:length(counts_over_qc)])
}

safe_HC <- safely(higher_criticism)

pvalues_qc <- pvalues %>% 
  select(suppdata_id, features, pvalues) %>% 
  mutate(HC = map(pvalues, safe_HC, breaks = 40, fdr = 0.05))
pvalues_qc_pass <- pvalues_qc %>% 
  mutate(QC_pass = map(HC, "result")) %>% 
  filter(!map_lgl(QC_pass, is.null)) %>% 
  unnest(QC_pass) %>% 
  select(suppdata_id, QC_pass)

megatab_qc <- megatab %>% 
  inner_join(pvalues_qc_pass) %>% 
  select(suppdata_id, Type, pi0, QC_pass)
megatab_qc %>%
  count(Type, QC_pass) %>% 
  spread(QC_pass, n) %>% 
  mutate(`%` = formattable::percent(`TRUE` / (`FALSE` + `TRUE`)))

qc_plot <- function(x, breaks) {
  x <- na.omit(x)
  b <- 1 / breaks
  qc <- qbinom(1 - b * 0.05, length(x), b)
  ggplot(data = NULL) +
    geom_histogram(aes(x), bins = breaks) +
    geom_hline(yintercept = qc, linetype = "dashed") +
    theme(axis.title = element_blank())
}

pvalues_qc_plots <- pvalues_qc %>% 
  mutate(QC_pass = map(HC, "result")) %>% 
  filter(!map_lgl(QC_pass, is.null)) %>% 
  unnest(QC_pass) %>% 
  select(-HC) %>% 
  mutate(p = map(pvalues, qc_plot, breaks = 40))

pvalues_qc_plots %>% 
  mutate(p_saved = map2(p, suppdata_id, ~ggsave(.x, filename = here("output/plots", str_c(str_extract(.y, "\\w+"), ".png")), width = 2, height = 1.34)))
