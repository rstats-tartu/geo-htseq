
# Run GEO query for HT sequencing data --------------------------------------------------------
## entrezquery can be installed from github
if(!require(entrezquery)) devtools::install_github("tpall/entrezquery")
library(entrezquery)
library(readr)

geo_query <- function(out_path) {
  
  ## ---- query-string ----
  query <- 'expression profiling by high throughput sequencing[DataSet Type]'
  
  ## ----- run-query -----
  # ds is short of document summaries
  ds <- entrez_docsums(query = query, db = "gds", retmax = 23000)
  
  ## ----- save downloaded filenames ----- 
  # create output dir if not present
  if (!dir.exists("output")) {
    dir.create("output")
  }
  
  write_rds(ds, path = out_path)
}

geo_query(snakemake@output[[1]])
