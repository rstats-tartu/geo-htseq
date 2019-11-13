library(tidyverse)
library(here)

megatab <- read_csv(here("output/megatab.csv"))
pvalues <- read_rds(here("output/pvalues.rds"))

higher_criticism <- function(x, breaks, fdr = 0.05) {
  x <- na.omit(x)
  b <- 1 / breaks
  h <- hist(x, breaks = seq(0, 1, b), plot = FALSE)
  counts <- h$counts
  mids <- h$mids
  qc <- qbinom(1 - b * fdr, length(x), b)
  counts_over_qc <- counts > qc
  if (all(!counts_over_qc)) {
    Class <- "uniform"
    QC_pass <- TRUE
  } else {
    mids_over_qc <- mids[counts_over_qc]
    end_of_peak <- rle(diff(counts_over_qc))$lengths[1] + 2
    QC_pass <- !any(counts_over_qc[end_of_peak:length(counts_over_qc)])
    if (QC_pass) {
      Class <- "anti-conservative"
    } else {
      end_of_peak <- rle(diff(rev(counts_over_qc)))$lengths[1] + 2
      conservative <- !any(rev(counts_over_qc)[end_of_peak:length(counts_over_qc)])
      if (conservative) {
        Class <- "conservative"
      } else {
        Class <- "other"
      }
    }
  }
  
  list(QC_pass, Class, mids_over_qc)
}

x <- runif(10000)
breaks <- 20
higher_criticism(x, breaks = breaks)
qc_plot(x, breaks = breaks)
higher_criticism(p.adjust(x, method = "fdr"), breaks = breaks)
qc_plot(p.adjust(x, method = "fdr"), breaks = breaks)
uniforms <- replicate(10000, higher_criticism(runif(10000), breaks = breaks)[[2]])
table(uniforms)

safe_HC <- safely(higher_criticism)

pvalues_qc <- pvalues %>% 
  select(suppdata_id, features, pvalues) %>% 
  mutate(HC = map(pvalues, safe_HC, breaks = breaks, fdr = 0.05))
pvalues_qc_pass <- pvalues_qc %>% 
  mutate(HC = map(HC, "result")) %>% 
  filter(!map_lgl(HC, is.null)) %>% 
  mutate(QC_pass = map_lgl(HC, 1),
         Class = map_chr(HC, 2),
         mids_over_qc = map(HC, 3)) %>% 
  select(suppdata_id, QC_pass, Class, mids_over_qc)

pvalues_qc_pass_plot <- pvalues_qc_pass %>% 
  select(suppdata_id, mids_over_qc) %>% 
  unnest(cols = c(mids_over_qc)) %>% 
  ggplot() +
  geom_bar(aes(mids_over_qc))

megatab_qc <- megatab %>% 
  inner_join(pvalues_qc_pass) %>% 
  select(suppdata_id, Type, pi0, QC_pass, Class)

megatab_qc %>%
  count(Type, QC_pass) %>% 
  spread(QC_pass, n) %>% 
  mutate(`%` = formattable::percent(`TRUE` / (`FALSE` + `TRUE`)))

megatab_qc %>%
  count(Type, Class) %>% 
  spread(Class, n)

qc_plot <- function(x, breaks) {
  x <- na.omit(x)
  b <- 1 / breaks
  qc <- qbinom(1 - b * 0.05, length(x), b)
  ggplot(data = NULL) +
    geom_histogram(aes(x), binwidth = 1 / breaks, center = 1 / (2 * breaks)) +
    geom_hline(yintercept = qc, linetype = "dashed") +
    theme(axis.title = element_blank())
}

pvalues_qc_plots <- pvalues_qc %>% 
  mutate(QC_pass = map(HC, "result")) %>% 
  filter(!map_lgl(QC_pass, is.null)) %>% 
  unnest(QC_pass) %>% 
  select(-HC) %>% 
  mutate(p = map(pvalues, qc_plot, breaks = breaks))

pvalues_qc_plots %>% 
  mutate(p_saved = map2(p, suppdata_id, ~ggsave(.x, filename = here("output/plots", str_c(str_extract(.y, "\\w+"), ".png")), width = 2, height = 1.34)))
