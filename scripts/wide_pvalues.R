library(tidyverse)
library(sparkline)
library(formattable)

#' ## P value stats
#' Import dataset
imported <- read_csv("output/parsed_suppfiles.csv")
imported <- imported %>% 
  mutate(accession = str_extract(id, "GSE\\d+")) %>% 
  select(accession, everything())

#' Pivot wide
pvalues <- imported %>% 
  filter(!is.na(Type)) %>% 
  mutate(Class = str_replace(Class, "conservative", "cons."),
         Conversion = str_replace(Conversion, "improvement", "impr."))
before <- pvalues %>% filter(Type == "raw")
after <- pvalues %>% filter(Type != "raw") %>% 
  rename_at(vars(Class, pi0, FDR_pval, hist), str_c, "_after") %>% 
  rename(filter_var = Type)
wide <- full_join(before, after) %>% 
  select(accession, id, starts_with("class"), Conversion, starts_with("pi0"), starts_with("hist"), filter_var, Set)


#' Generate html table: beware of very large table!
tab <- wide %>% 
  mutate_at(vars(starts_with("hist")), ~map(.x, ~as.integer(unlist(str_extract_all(.x, "\\d+"))))) %>% 
  mutate_at(vars(starts_with("hist")), ~map(.x, spk_chr, type = "bar")) %>%
  formattable() %>%
  formattable::as.htmlwidget() %>%
  spk_add_deps()

#' Conversion
wide %>% 
  count(Class, Class_after) %>% 
  knitr::kable() %>% 
  kableExtra::kable_styling()

#' ## Experiment metadata
#' Importing experiment metadata (dates and etc) 
document_summaries <- read_csv("output/document_summaries.csv")

#' Importing publications metadata
publications <- read_csv("output/publications.csv") %>% 
  rename(PubMedIds = Id)

#' Single-cell experiments
single_cell <- read_csv("output/single-cell.csv")

#' Citations
citations <- read_csv("output/scopus_citedbycount.csv")

#' Merging publications with citations
pubs <- publications %>% 
  left_join(citations) %>% 
  select(PubMedIds, PubDate, Source, FullJournalName, ISSN, ESSN, citations) %>% 
  mutate(ISSN = case_when(
    is.na(ISSN) ~ ESSN,
    TRUE ~ ISSN
   ))

#' ## Updating P values with some metadata
pvals_wide <- document_summaries %>% 
  select(accession = Accession, PDAT) %>% 
  inner_join(wide) %>% 
  mutate(single_cell = accession %in% single_cell$Accession) %>% 
  mutate(platform = case_when(
    filter_var == "basemean" ~ "deseq",
    filter_var == "aveexpr" ~ "limma",
    filter_var == "logcpm" ~ "edger",
    filter_var == "fpkm" & str_detect(Set, "p_value") ~ "cuffdiff",
    TRUE ~ "other"
  ))

pvals_wide %>% 
  count(platform)
