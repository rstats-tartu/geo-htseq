
# Load libs and settings --------------------------------------------------
source("R/_common.R")

# GEO query results and document summaries --------------------------------
ds <- read_rds("output/document_summaries.rds")

# Load series matrix data frames ----------------------------
gsem <- read_rds("output/gsem.rds")

# Identify and filter out Accessions with missing gsematrices
gsem_missing_or_faulty <- gsem %>% 
  filter(map_lgl(series_matrix, ~class(.x) != "ExpressionSet"))
gsem <- gsem %>% 
  filter(map_lgl(series_matrix, ~class(.x) == "ExpressionSet"))

## Read in local supplemental tables ---------------------------------------

# File name extensions in downloaded supplementary files
supptabs <- geofile_df(local_suppfile_folder, "suppfiles")

# Unzip gz xls(x)? files, keep originals
need_to_unzip_xls_files <- FALSE
if (need_to_unzip_xls_files) {
  xls_gz <- filter(supptabs, str_detect(suppfiles, "xls(x)?.gz$"))$suppfiles
  xls_gz <- file.path(local_suppfile_folder, xls_gz)
  system(sprintf("gzip -dk %s", paste(xls_gz, collapse = " ")))
}

# List files again to exclude compressed xls files
supptabs <- geofile_df(local_suppfile_folder, "suppfiles") %>% 
  filter(!str_detect(suppfiles, "xls(x)?.gz$"))

# Duplicated files
supptabs_duplicated <- supptabs %>% 
  mutate(files = str_replace(suppfiles, "\\.(gz|bz2)$","")) %>% 
  filter(duplicated(files))

# Remove non tabular filetypes from downloaded files
out_string1 <- c("filelist","annotation","readme","error","raw.tar","csfasta",
                 "bam","sam","bed","[:punct:]hic","hdf5","bismark","map",
                 "barcode","peaks")
out_string2 <- c("tar","gtf","(big)?bed(\\.txt|12|graph|pk)?","bw","wig",
                 "hic","gct(x)?","tdf","gff(3)?","pdf","png","zip","sif",
                 "narrowpeak","fa", "r$", "rda(ta)?$")

supptabs <- supptabs %>% 
  filter(!str_detect(tolower(suppfiles), str_c(out_string1, collapse = "|")),
         !str_detect(tolower(suppfiles), str_c(out_string2, "(\\.gz|\\.bz2)?$", 
                                                   collapse = "|")))

# Merge gsem ExpressionSets to supptabs
supptabs <- select(gsem, Accession, series_matrix) %>% 
  nest(series_matrix, .key = "matrixfiles") %>% 
  left_join(supptabs, .)

# Use only files 
supptabs <- supptabs %>% 
  filter(Accession %in% filter(ds, PDAT <= last_date)$Accession)

# Files that crash R session
bad <- c("GSE93374_Merged_all_020816_DGE.txt.gz", 
         "GSE88931_RNA-seq_MergedReadCounts.tsv.gz",
         "GSE55385_transcripts_GSE.tsv.gz",
         "GSE74549_ChIP_1kbWindows_correctedReadCount.txt.gz",
         "GSE77213_Nguyen_GEO_TN03_tallies.xls",
         "GSE77213_Nguyen_GEO_TN05_tallies_total.xls",
         "GSE60012_100bpTiles_RRBS_Mouse.txt.gz",
         "GSE60415_heatmap-upload.xls",
         "GSE67516_RNA_seq_rep1_diffExp_analysis.xls")

# Remove 'bad' files
supptabs <- supptabs %>% filter(!(suppfiles %in% bad))

source("lib/munge_geo.R")
source("lib/checkFullRank.R")
source("lib/text_funs.R")
library(Biobase)

import_supptabs <- TRUE
if(import_supptabs) {
  start <- Sys.time()
  st <- supptabs %>% 
    mutate(result = map(suppfiles, ~ try(munge_geo_pvalue(file.path(local_suppfile_folder, .x)))))
  end <- Sys.time()
  time_diff <- start-end
  
  system(paste("echo ' Run took", abs(time_diff), units(time_diff), "' >> log.txt"))
  saveRDS(st, file = "output/suppdata.rds")
  
}
