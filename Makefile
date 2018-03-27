#!/bin/bash

suppfiles=$(cat output/suppfilenames_filtered.txt)
suppfiles_rds=$(ls output/suppl | sed 's/$/.rds/g')

_main.html: index.Rmd 01_introduction.Rmd 02_methods.Rmd 03_results.Rmd \
04_discussion.Rmd 05_references.Rmd output/document_summaries.rds \
output/gsem.rds output/suppdata.rds output/suppfilenames.rds \
output/publications.rds output/document_summaries.rds
	Rscript -e "rmarkdown::render_site(input = 'index.Rmd', output_format = 'bookdown::html_document2')"

output/suppdata.rds: R/_common.R output/document_summaries.rds output/gsem.rds lib/munge_geo.R lib/checkFullRank.R lib/text_funs.R
	Rscript R/06_munge_suppfiles.R

output/gsem.rds: R/_common.R
	Rscript R/05_munge_series_matrixfiles.R

output/suppl_rds/{echo $suppfiles_rds}: echo $suppfiles
	Rscript R/munge_geo_pvalue_to_rds.R

echo $suppfiles: output/suppfilenames_filtered.rds
	Rscript R/04_download_suppfiles.R

output/supfilenames_filtered.rds output/supfilenames_filtered.txt: R/_common.R output/suppfilenames.rds R/out_strings.R
	Rscript R/03_filter_suppfilenames.R

output/suppfilenames.rds: output/document_summaries.rds R/02_download_suppfilenames.R R/_common.R
	Rscript R/02_download_suppfilenames.R

output/publications.rds: R/_common.R output/document_summaries.rds
	Rscript R/07_download_publications.R

output/document_summaries.rds: R/01_GEO_query.R
	Rscript R/01_GEO_query.R
