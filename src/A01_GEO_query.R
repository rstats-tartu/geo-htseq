
# Run GEO query for HT sequencing data --------------------------------------------------------
library(dplyr)

# @knitr query-string
query <- 'expression profiling by high throughput sequencing[DataSet Type]'

# @knitr run-query
# Load helper functions
source("lib/get_GEO_funs.R")

# get Ids
Ids <- get_ids(query, retmax = 15000)

# Get query summaries -----------------------------------------------------
sumcont <- get_docsums(Ids)

# Gateway error in previous step may fuck up this step!!!, if this occurs, please just rerun.
ds <- lapply(sumcont, extract_docsums) %>% bind_rows()

# ds is short of document summaries
save(ds, file = "data/ds.RData")
