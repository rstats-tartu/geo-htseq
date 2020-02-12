library(tidyverse)
library(sparkline)
library(formattable)

imported <- read_csv("output/parsed_suppfiles.csv")
imported <- imported %>% 
  mutate(accession = str_extract(id, "GSE\\d+")) %>% 
  select(accession, everything())

p_value_sets <- imported %>% 
  count(Type) %>% 
  mutate(group = case_when(
    Type == "raw" ~ "raw",
    is.na(Type) ~ "no-pvalues",
    TRUE ~ "filtered"
  )) %>% 
  group_by(group) %>% 
  summarise(n = sum(n))

distinct_tables <- imported %>% 
  select(accession, id, Type) %>% 
  filter(is.na(Type) | Type == "raw") %>% 
  distinct()

geos_with_pvalues <- distinct_tables %>% 
  group_by(accession) %>% 
  summarise(pvalues = if_else(any(Type == "raw"), "yes", "no")) %>% 
  mutate(pvalues = case_when(
    is.na(pvalues) ~ "no",
    TRUE ~ pvalues
  )) %>% 
  count(pvalues) 
  
files_with_pvalues <- distinct_tables %>% 
  mutate(file = case_when(
    str_detect(id, "from|sheet") ~ str_replace(id, "^.*(from |sheet-)", ""),
    TRUE ~ id
  )) %>% 
  group_by(accession, file) %>% 
  summarise(pvalues = if_else(any(Type == "raw"), "yes", "no")) %>% 
  ungroup() %>% 
  mutate(pvalues = case_when(
    is.na(pvalues) ~ "no",
    TRUE ~ pvalues
  )) %>% 
  distinct() %>% 
  count(pvalues) 


imported %>% 
  filter(Type == "raw") %>% 
  select(id, Class, Conversion) %>% 
  count(Class, Conversion)



tab_pipe <- . %>% 
  mutate(hist = map(hist, ~as.integer(unlist(str_extract_all(.x, "\\d+")))), 
         hist = map(hist, spk_chr, type = "bar")) %>%
  formattable() %>%
  formattable::as.htmlwidget() %>%
  spk_add_deps()
