
# Run GEO query for HT sequencing data --------------------------------------------------------
## Library can be installed from github
# devtools::install_github("tpall/entrezquery")
library(entrezquery)

## ---- query-string ----
query <- 'expression profiling by high throughput sequencing[DataSet Type]'

## ----- run-query -----
ds <- entrez_docsums(query = query, db = "gds", retmax = 15000)

# ds is short of document summaries
save(ds, file = "data/ds.RData")

