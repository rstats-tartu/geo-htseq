#!/usr/bin/env Rscript
xls <- commandArgs(trailingOnly = TRUE)

# Load libs and settings --------------------------------------------------
source("R/_common.R")
library(GEOquery)

# GEO query results and document summaries --------------------------------
ds <- read_rds("output/document_summaries.rds")

# Load series matrix data frames ----------------------------
gsem <- read_rds("output/gsem.rds")

# Extract series matrixes
gsem <- gsem %>% 
  mutate(series_matrix = map(gse, "result"))

# Identify and filter out Accessions with missing gsematrices
gsem_missing_or_faulty <- gsem %>% 
  filter(map_lgl(series_matrix, ~class(.x) != "ExpressionSet"))

gsem <- gsem %>% 
  filter(map_lgl(series_matrix, ~class(.x) == "ExpressionSet"))

## Read in local supplemental tables ---------------------------------------

# File name extensions in downloaded supplementary files
supptabs <- geofile_df(local_suppfile_folder, "suppfiles")

# Unzip gz xls(x)? files, keep originals
xls_gz <- filter(supptabs, str_detect(suppfiles, "xls(x)?.gz$"))$suppfiles
xls_gz <- file.path(local_suppfile_folder, xls_gz)

# List files again to exclude compressed xls files
supptabs <- geofile_df(local_suppfile_folder, "suppfiles") %>% 
  filter(!str_detect(suppfiles, "xls(x)?.gz$"))

# Duplicated files
supptabs_duplicated <- supptabs %>% 
  mutate(files = str_replace(suppfiles, "\\.(gz|bz2)$","")) %>% 
  filter(duplicated(files))

# Merge gsem ExpressionSets to supptabs
supptabs <- select(gsem, Accession, series_matrix) %>% 
  nest(series_matrix, .key = "matrixfiles") %>% 
  left_join(supptabs, .)

# Files that crash R session
bad <- c("GSE93374_Merged_all_020816_DGE.txt.gz", 
         "GSE88931_RNA-seq_MergedReadCounts.tsv.gz",
         "GSE55385_transcripts_GSE.tsv.gz",
         "GSE74549_ChIP_1kbWindows_correctedReadCount.txt.gz",
         "GSE77213_Nguyen_GEO_TN03_tallies.xls",
         "GSE77213_Nguyen_GEO_TN05_tallies_total.xls",
         "GSE60012_100bpTiles_RRBS_Mouse.txt.gz",
         "GSE60415_heatmap-upload.xls",
         "GSE67516_RNA_seq_rep1_diffExp_analysis.xls",
         "GSE53298_processed_data.xlsx",
         "GSE53350_SAGE_KomatsuFurukawa_Processed.xls",
         "GSE60483_AQ_exp.tab.gz",
         "GSE53260_Supplementary_S2.xls",
         "GSE53260_Supplementary_S1.xls")

# Remove 'bad' files
supptabs <- supptabs %>% filter(!(suppfiles %in% bad))

source("R/munge_geo.R")
source("R/checkFullRank.R")
source("R/text_funs.R")

if (xls) {
  st <- supptabs %>%
    filter(str_detect(suppfiles, "xls(x)?")) %>% 
    mutate(result = map(suppfiles, ~ try(munge_geo_pvalue(file.path(local_suppfile_folder, .x)))))
  write_rds(st, path = "output/suppdata_xls.rds")
} else {
  st <- supptabs %>%
    filter(!str_detect(suppfiles, "xls(x)?")) %>% 
    mutate(result = map(suppfiles, ~ try(munge_geo_pvalue(file.path(local_suppfile_folder, .x)))))
  write_rds(st, path = "output/suppdata_allothers.rds")
}
