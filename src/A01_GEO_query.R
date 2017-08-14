
# Run GEO query for HT sequencing data --------------------------------------------------------
library(tidyverse)
library(magrittr)

# @knitr query-string
query <- 'expression profiling by high throughput sequencing[DataSet Type] AND (Homo sapiens[Organism] OR Mus musculus[Organism])'

# @knitr run-query
# Load helper functions
source("lib/get_GEO_funs.R")

# get Ids
Ids <- get_ids(query, retmax = 13000)

# Get query summaries -----------------------------------------------------
sumcont <- get_docsums(Ids)

# XML stuff ------------------------------------------------------

# extract docsums data much much faster
# library(XML)

# Gateway error in previous step may fuck up this step!!!, if this occurs, please just rerun.
ds <- sumcont %>% lapply(extract_docsums) %>% bind_rows()
# ds is short of document summaries
save(ds, file = "data/ds.RData")
