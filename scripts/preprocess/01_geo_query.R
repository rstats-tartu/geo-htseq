
# Run GEO query for HT sequencing data --------------------------------------------------------
## entrezquery can be installed from github
if(!require(entrezquery)) devtools::install_github("tpall/entrezquery")
library(entrezquery)
library(readr)

## ---- query-string ----
query <- 'expression profiling by high throughput sequencing[DataSet Type]'
  
## ----- run-query -----
# ds is short of document summaries
httr::set_config(httr::config(http_version = 0, timeout = 6000))
ds <- entrez_docsums(query = query, db = "gds", retmax = 30000)
  
## ----- save downloaded filenames ----- 
# create output dir if not present
if (!dir.exists("output")) dir.create("output")

write_rds(ds, path = snakemake@output[[1]])
