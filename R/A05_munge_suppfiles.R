## Load libs ---------------------------------------------------------------

library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(Biobase)
source("lib/helpers.R")
source("lib/munge_geo.R")
source("lib/checkFullRank.R")
source("lib/text_funs.R")

## GEO query results and document summaries --------------------------------
## source("src/A01_GEO_query.R")
if(!"ds" %in% ls()){
  load("data/ds.RData")
}

## Load series matrix data frames ----------------------------
update_geoseriesmatrix_files <- FALSE
if(update_geoseriesmatrix_files){
  source("src/A04_munge_series_matrixfiles.R")
} else {
  load("data/gsem.RData")
}

## Identify and filter out Accessions with missing gsematrices
gsem_missing_or_faulty <- filter(gsem, map_lgl(gsematrix, ~class(.x) != "ExpressionSet")) 
gsem <- filter(gsem, map_lgl(gsematrix, ~class(.x) == "ExpressionSet"))

## Read in local supplemental tables ---------------------------------------

## File name extensions in downloaded supplementary files
local_suppfile_folder <- "/Volumes/Media/srp_example/data/counts" # "data/counts"
supptabs <- geofile_df(local_suppfile_folder, "suppfiles")

## Unzip gz xls(x)? files, keep originals
need_to_unzip_xls_files <- FALSE
if(need_to_unzip_xls_files){
  xls_gz <- filter(supptabs, str_detect(suppfiles, "xls(x)?.gz$"))$suppfiles
  xls_gz <- file.path(local_suppfile_folder, xls_gz)
  system(sprintf("gzip -dk %s", paste(xls_gz, collapse = " ")))
}

## List files again to exclude compressed files
supptabs <- geofile_df(local_suppfile_folder, "suppfiles") %>% 
  filter(!str_detect(suppfiles, "xls(x)?.gz$"))

## Duplicated files
supptabs_duplicated <- supptabs %>% 
  mutate(files = str_replace(suppfiles, "\\.(gz|bz2)$","")) %>% 
  filter(duplicated(files))

## One more time. Remove non tabular filetypes from downloaded files
non_tabular_formats <- c(".r",".rdata",".xml",".bed12",".json",".rda",".sam",
                         ".fasta",".html",".tar",".hdf5",".psl")

supptabs <- filter(supptabs, !(get_filext(suppfiles) %in% non_tabular_formats))

## Files that may crash R session
bad <- c("GSE93374_Merged_all_020816_DGE.txt.gz", 
         "GSE88931_RNA-seq_MergedReadCounts.tsv.gz",
         "GSE55385_transcripts_GSE.tsv.gz",
         "GSE74549_ChIP_1kbWindows_correctedReadCount.txt.gz",
         "GSE77213_Nguyen_GEO_TN03_tallies.xls",
         "GSE77213_Nguyen_GEO_TN05_tallies_total.xls",
         "GSE60012_100bpTiles_RRBS_Mouse.txt.gz",
         "GSE60415_heatmap-upload.xls",
         "GSE67516_RNA_seq_rep1_diffExp_analysis.xls")

## Remove 'bad' files
supptabs <- filter(supptabs, !(suppfiles %in% bad))

## Merge gsem ExpressionSets to supptabs
supptabs <- select(gsem, Accession, gsematrix) %>% 
  nest(gsematrix, .key = "matrixfiles") %>% 
  left_join(supptabs, .)

import_supptabs <- TRUE
if(import_supptabs){
  start <- Sys.time()
  st <- mutate(supptabs, result = map(suppfiles, ~ try(munge_geo2(.x, dir = local_suppfile_folder))))
  end <- Sys.time()
  start-end
  
  # Time difference of -1.556035 hours
  save(st, file = "data/st.RData")
  # load("data/st.RData")
}
