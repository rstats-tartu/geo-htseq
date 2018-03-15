library(dplyr)
glue::glue('
#!/bin/bash\n
output/suppdata.rds: R/_common.R output/document_summaries.rds output/gsem.rds lib/munge_geo.R lib/checkFullRank.R lib/text_funs.R\n
\tRscript R/A06_munge_suppfiles.R\n
output/gsem.rds: R/_common.R\n
\tRscript R/A05_munge_series_matrixfiles.R\n
{readr::read_lines("output/suppfilenames_filtered.txt")}: output/suppfilenames_filtered.rds\n
\tRscript R/A04_download_suppfiles.R\n
output/supfilenames_filtered.rds output/supfilenames_filtered.txt: R/_common.R output/suppfilenames.rds R/out_strings.R\n
\tRscript R/A03_filter_suppfilenames.R\n
output/suppfilenames.rds: output/document_summaries.rds R/A02_download_suppfilenames.R R/_common.R\n
\tRscript R/A02_download_suppfilenames.R\n
output/document_summaries.rds: R/A01_GEO_query.R\n
\tRscript R/A01_GEO_query.R'
           ) %>% 
  readr::write_lines(path = 'Makefile')
