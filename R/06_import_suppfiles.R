
# Load libs, local_suppfile_folder, and helper functions
source("R/_common.R")
source("R/munge_geo.R")

library(digest)
library(glue)
library(GEOquery)

import_suppdata <- function(supptab, out_path) {
  supptab <- read_rds(supptab)
  bad <- c("GSE121228_gene_expression_anno.txt.gz", "GSE113074_Raw_combined.annotated_counts.tsv.gz", "GSE113074_Corrected_combined.annotated_counts.tsv.gz", "GSE121737_early_and_medium_bud.repGene.txt.gz")
  supptab <- supptab %>% filter(!(suppfiles %in% bad))
  supptab %>%
    mutate(result = map(suppfiles, ~ try(munge_geo_pvalue(file.path(local_suppfile_folder, .x))))) %>% 
    write_rds(path = out_path)
}

import_suppdata(snakemake@input[[1]], snakemake@output[[1]])

