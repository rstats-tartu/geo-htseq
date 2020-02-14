library(tidyverse)
library(sparkline)
library(formattable)

#' Import dataset
imported <- read_csv("output/parsed_suppfiles.csv")
imported <- imported %>% 
  mutate(accession = str_extract(id, "GSE\\d+")) %>% 
  select(accession, everything())

#' Pivo wide
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
  count(Class, Class_after)

