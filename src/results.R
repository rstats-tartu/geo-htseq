
library(tidyverse)
library(lubridate)
library(stringr)
library(magrittr)
library(formattable)
library(Biobase)
library(grid)
library(gridExtra)

# Load helper functions
source("lib/helpers.R")

# Set panel label case
panel_label_case <- "upper"


# RNA-seq-dynamics --------------------------------------------------------

## @knitr rna-seq-dynamics

library(tidyverse)
library(lubridate)
library(stringr)
library(magrittr)
library(formattable)

## src/A01b_GEO_RNA-seq-dynamics.R
load("data/GEO_RNA-seq-dynamics.RData") # all GEO HT-seq expr datasets
ds_all <- ds

## src/A01_GEO_query.R
load("data/ds.RData") # mouse and human GEO HT-seq expr datasets

ds_merge <- bind_rows(ds_all, ds, .id = "id")

pdat <- ds_merge %>% 
  select(id, PDAT, PubMedIds) %>% 
  mutate(pub = nchar(PubMedIds)!=0) %>% 
  group_by(id, PDAT) %>% 
  summarise(N = n(),
            pub = sum(pub)) %>% 
  mutate_at(vars(N, pub), cumsum)

pdat %<>% gather(key, value, -PDAT, -id) 

pdat %<>% 
  ungroup() %>% 
  mutate(id = if_else(id==1, "All taxa", "Human and mouse"))

geop <- pdat %>% 
  ggplot(aes(ymd(PDAT), value, linetype=key)) + 
  geom_line() +
  facet_wrap(~id) +
  ylab("Number of GEO series") +
  xlab("Publication date") +
  scale_linetype_discrete(labels=c("All series","Series with\npublications")) +
  theme(legend.position = c(0.75, 0.77),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank())

## formattable::percent()
perc_mmhs <- percent(round(nrow(ds)/nrow(ds_all), 1), digits = 0)
ppub <- group_by(pdat, id, key) %>% 
  summarise_at("value", max) %>% 
  group_by(id) %>% 
  summarise(ppub = value[key=="pub"]/value[key=="N"]) %>% 
  ungroup %>% 
  summarise(ppub = mean(ppub)) %>% 
  .$ppub %>% 
  round(digits=1) %>% 
  percent(digits=0)


# Missingsuppfiles --------------------------------------------------------
# In this section we download and analyse supplementary file names

## @knitr missingsuppfiles
# src/A02_download_suppfiles.R
load("data/suppfilenames_2017-06-19.RData")

failed_suppfiles <- suppfilenames %>% 
  filter(map_lgl(SuppFileNames, ~inherits(.x,"try-error"))) %>% 
  arrange(PDAT) %>% 
  group_by(PDAT) %>%
  summarise(N = n()) %>% 
  mutate(N = cumsum(N))

failed_suppfiles <- suppfilenames %>% 
  filter(map_lgl(SuppFileNames, ~inherits(.x,"try-error"))) %>% 
  select(PDAT, PubMedIds) %>% 
  mutate(pub = nchar(PubMedIds)!=0) %>% 
  arrange(PDAT) %>% 
  group_by(PDAT) %>% 
  summarise(N = n(),
            pub = sum(pub)) %>% 
  mutate_at(vars(N, pub), cumsum) %>% 
  gather(key, value, -PDAT) 

failed_suppfiles$id <- "Human or mouse series\nlacking supplementary files"

fsupp <- ggplot(failed_suppfiles, aes(ymd(PDAT), value, linetype = key)) +
  geom_line() +
  facet_wrap(~id) +
  xlab("Publication date") +
  ylab("Number of GEO series") +
  scale_linetype_discrete(labels=c("All series", "Series with\npublications")) +
  theme(legend.position = c(0.5, 0.77),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank())

no_suppfile <- suppfilenames %>% 
  mutate(has_suppfile = !map_lgl(SuppFileNames, ~inherits(.x,"try-error"))) %>% 
  select(has_suppfile) %>% 
  table

perc_wsuppfile <- percent(round(1-(no_suppfile[1]/sum(no_suppfile)), 2), digits = 0)


# Queryfig ----------------------------------------------------------------

## @knitr queryfig

pg <- lapply(list(geop, fsupp), ggplotGrob)
pg <- add_labels(pg, case = panel_label_case)
pga <- arrangeGrob(grobs = pg, ncol = 2, widths = c(2, 1))
grid.draw(pga)


# Commonfilenames ---------------------------------------------------------

## @knitr commonfilenames

# Single most common filename: filelist.txt
most_common_filename <- suppfilenames %>% 
  filter(!map_lgl(SuppFileNames, ~inherits(.x,"try-error"))) %>% 
  unnest(SuppFileNames) %>% 
  group_by(SuppFileNames) %>% 
  summarise(N=n())

# Supplemental file names with more than N=10 occurences
cf <- suppfilenames %>% 
  filter(!map_lgl(SuppFileNames, ~inherits(.x,"try-error"))) %>% 
  unnest(SuppFileNames) %>%
  mutate(common_filenames = str_replace(SuppFileNames, "GSE[0-9]*_", ""),
         common_filenames = str_replace(common_filenames, "\\.gz$", ""), # remove also .gz
         common_filenames = tolower(common_filenames))

cfn <- group_by(cf, common_filenames) %>% 
  summarise(N=n()) %>% 
  arrange(desc(N)) %>% 
  filter(N > 10)

cfp <- ggplot(cfn, aes(common_filenames, N)) +
  geom_point() +
  scale_x_discrete(limits = rev(cfn$common_filenames)) +
  scale_y_log10() +
  coord_flip() + 
  xlab("Common stubs of SuppFileNames\n(>10 occurences) ") +
  ylab("Number of files")

# plot commonfilenames ggplot
cfp

# Raw data as supplemental file plot (decide wether to show)
cfraw <- filter(cf, str_detect(common_filenames, "raw.tar")) %>%
  arrange(PDAT) %>%
  group_by(PDAT) %>%
  summarise(N = n()) %>%
  mutate(N = cumsum(N)) %>%
  ggplot(aes(ymd(PDAT), N, group = 1)) +
  geom_line()

# Percent GEO ids submit raw data as supplemental file
supp_raw_perc <- percent(n_distinct(filter(cf, str_detect(common_filenames, "raw.tar"))$Accession)/n_distinct(suppfilenames$Accession), 0)

# Filter downloaded supplementary file names ------------------------------
# In this section, we filter supplementary file names for patterns: 
# we are looking anly for tabular data. 
suppfiles_of_interest <- suppfilenames %>% 
  unnest(SuppFileNames) %>%
  filter(!str_detect(tolower(SuppFileNames), out_string1),
         !str_detect(tolower(SuppFileNames), paste0(out_string2, "(\\.gz|\\.bz2)?$", collapse = "|"))) %>% 
  select(Accession, SuppFileNames, FTPLink, PDAT) %>% 
  mutate(filext = str_extract(tolower(SuppFileNames), "\\.[:alpha:]+([:punct:][bgz2]+)?$")) 

suppf_oi_perc <- percent(1-(n_distinct(suppfiles_of_interest$Accession)/n_distinct(suppfilenames$Accession)), 0)



