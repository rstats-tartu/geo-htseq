
# Run GEO query for HT sequencing data --------------------------------------------------------
## Library can be installed from github
# devtools::install_github("tpall/entrezquery")
library(entrezquery)

## ---- query-string ----
query <- 'expression profiling by high throughput sequencing[DataSet Type]'

## ----- run-query -----
# ds is short of document summaries
ds <- entrez_docsums(query = query, db = "gds", retmax = 16000)

## ----- save downloaded filenames ----- 
# create output dir if not present
if (!dir.exists("output")) {
  dir.create("output")
}

saveRDS(ds, file = "output/document_summaries.rds")

