
library(tidyverse)
library(magrittr)
library(RCurl)

# Load helper functions
source("lib/get_GEO_funs.R")

# Query
query <- 'expression profiling by high throughput sequencing[DataSet Type]'

# Get Ids
Ids <- get_GEO_Ids(query = query, retmax = 13000)

# Get query summaries -----------------------------------------------------

sumcont <- get_GEO_DocSums(Ids)

ds <- sumcont %>% lapply(extract_gds_docsums) %>% bind_rows()
save(ds, file = "data/GEO_RNA-seq-dynamics.RData")
geo <- ds$PDAT %>% table %>% cumsum() %>% data_frame(date=names(.), studies=.)

library(lubridate)
ggplot(geo, aes(ymd(date), studies, group=1)) +
  geom_line() +
  ylab("Number of GEO series") +
  xlab("Date")

pdat <- ds %>% 
  select(PDAT, PubMedIds) %>% 
  mutate(pub=nchar(PubMedIds)!=0) %>% 
  group_by(PDAT) %>% 
  summarise(N = n(),
            pub = sum(pub)) %>% 
  mutate_at(vars(N, pub), cumsum)
  
pdat %<>% gather(key,value, -PDAT) 

pdat %>% 
  ggplot(aes(ymd(PDAT), value, linetype=key)) + 
  geom_line() +
  ylab("Number of GEO series") +
  xlab("Date") +
  scale_linetype_discrete(labels=c("GEO series","GEO series with\npublications")) +
  theme(legend.position = c(.35, 0.8),
        legend.background = element_blank(),
        legend.title = element_blank())

