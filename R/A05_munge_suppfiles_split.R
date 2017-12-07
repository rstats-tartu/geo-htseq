
# Load libs ---------------------------------------------------------------

source("R/_common.R")

# GEO query results and document summaries --------------------------------
# source("src/A01_GEO_query.R")
if(!"ds" %in% ls()){
  load("data/ds.RData")
}

# Load series matrix data frames ----------------------------
update_geoseriesmatrix_files <- FALSE
if(update_geoseriesmatrix_files){
  source("src/A04_munge_series_matrixfiles.R")
} else {
  load("data/gsem.RData")
}

# Identify and filter out Accessions with missing gsematrices
gsem_missing_or_faulty <- gsem %>% 
  filter(map_lgl(gsematrix, ~class(.x) != "ExpressionSet"))
gsem <- gsem %>% 
  filter(map_lgl(gsematrix, ~class(.x) == "ExpressionSet"))

## Read in local supplemental tables ---------------------------------------

# File name extensions in downloaded supplementary files
local_suppfile_folder <- "/Volumes/Media/srp_example/data/counts" # "data/counts"
supptabs <- geofile_df(local_suppfile_folder, "suppfiles")

# Unzip gz xls(x)? files, keep originals
need_to_unzip_xls_files <- FALSE
if(need_to_unzip_xls_files){
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

# Files that may crash R session
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

# Merge gsem ExpressionSets to supptabs
supptabs <- select(gsem, Accession, gsematrix) %>% 
  nest(gsematrix, .key = "matrixfiles") %>% 
  left_join(supptabs, .)

source("lib/munge_geo.R")
source("lib/checkFullRank.R")
source("lib/text_funs.R")
library(Biobase)

import_supptabs <- FALSE
if(import_supptabs){
  start <- Sys.time()
  st <- mutate(supptabs, result = map2(matrixfiles, suppfiles, ~ try(munge_geo(.x, .y, dir = local_suppfile_folder))))
  end <- Sys.time()
  start-end
  
  # Time difference of -1.556035 hours
  save(st, file = "data/st.RData")
  # load("data/st.RData")
}

# # P values ----------------------------------------------------------------
# tabs <- filter(res, !map_lgl(result, is.data.frame))
# tabs <- tabs[rep(1:nrow(tabs), mutate(tabs, n = map_int(result, length))$n), 
#              c("Accession", "countfiles", "PDAT", "excelfiles")] %>%
#   bind_cols(data_frame(result = unlist(tabs$result, recursive = FALSE))) %>% 
#   filter(!map_lgl(result, is.null))
# 
# pvalues <- filter(tabs, map_lgl(result, is.tibble))
# pvalues %>% 
#   filter(map_lgl(countfiles, str_detect, pattern = "xls"))
# 
# pvalues %>% 
#   mutate(features = map_int(result, nrow)) %>% 
#   ggplot(aes(features)) +
#   geom_histogram(bins=30) +
#   scale_x_log10() +
#   geom_vline(xintercept = 8000, linetype = 3) +
#   xlab("Number of P values")
# 
# # Fit models -----------------------------------------------------
# # Filter ExpressionSets
# esets <- tabs %>% filter(map_chr(result, ~class(.x)[1])=="ExpressionSet")
# 
# # Identify vars for model fitting
# esets %<>% mutate(groups = map(result, ~get_model_eset(pData(.x))))
# 
# library(edgeR)
# # GSE91063_Torrent_run_1_RPKM_cutoff.xlsx has only 3720 reads
# esets %<>% 
#   filter(map_lgl(result, ~nrow(exprs(.x))>nrowthreshold)) %>% 
#   mutate(dge = pmap(list(result, groups, countfiles), function(x,y,z) {message(z); DGEList(counts=exprs(x), samples = y, genes = fData(x))})) %>% 
#   filter(!map_lgl(groups, is.null))
# 
# # Make design formulas and model.matrix
# make_formula <- . %>% colnames %>% paste(collapse="+") %>% paste("~",.) %>% formula()
# esets %<>% mutate(design = map(groups, make_formula)) 
# esets %<>% mutate(design = map2(design, groups, ~model.matrix(object=.x, data=.y))) 
# 
# # Filter genes with low counts: 3 or more samples with cpm of 10/L. L, # millions of counts in the smallest library
# esets %<>% mutate(keep = map(dge, ~rowSums(cpm(.x)>10/min(.x$samples$lib.size)/1e6) >= 3),
#                   dge_filt = map2(dge, keep, ~.x[.y,]))
# 
# # Recalculate library sizes and calculate normalisation factors
# esets %<>% mutate(dge_filt = map(dge_filt, ~{.x$samples$lib.size <- colSums(.x$counts); .x}))
# esets %<>% mutate(dge_filt = map2(countfiles, dge_filt, ~{message(.x); calcNormFactors(.y)}))
# esets %<>% mutate(dge_filt = pmap(list(dge_filt, design, countfiles), 
#                                   function(x,y,z) {message(z); try(estimateDisp(x,y))}))
# 
# # Fit models
# esets %<>% 
#   filter(!map_lgl(dge_filt, ~inherits(.x, "try-error"))) %>%
#   mutate(fit = pmap(list(dge_filt, design, countfiles), function(x,y,z) {message(z); try(glmFit(x,y))})) %>% 
#   filter(!map_lgl(fit, ~inherits(.x, "try-error"))) %>%
#   mutate(results = map(fit, ~glmLRT(.x, coef=2:ncol(.x$design))),
#          toptags = map(results, topTags, n=Inf))
# save(esets, file="data/edgeR_fits_esets.RData")
# # load("data/edgeR_fits_esets.RData")
# # Extract edgeR p values
# 
# edger_pvalues <- esets %>% 
#   mutate(tt = map(toptags, ~.x@.Data[[1]]),
#          pvalue = map(tt, "PValue")) %>% 
#   select(Accession, countfiles, PDAT, pvalue)
# 
# # Filter submitter provided p values
# # pv <- pvalues %>% filter(countfiles=="GSE78078_Processed_file_cd38_negative_kd_oe.xlsx") %$% result[[1]]
# basemean_filter <- function(pv, basemeanthresh=10) {
#   
#   colnames(pv) <- make.names(colnames(pv), unique = T)
#   if(!all(is.na(pv$basemean))){
#     pv <- filter(pv, basemean>basemeanthresh)
#   }
#   pv %>% 
#     select(-basemean) %>% 
#     filter(complete.cases(.)) %>% 
#     as.list()
# }
# 
# pvalues %>% filter(map_lgl(result, ~ncol(.x)>2))
# 
# filter_pvalues <- pvalues %>% 
#   mutate(pvalue = map2(result, countfiles, ~{message(.y); basemean_filter(.x)}),
#          n=map_int(pvalue, ~length(.x)))
# 
# submitter_pvalues <- bind_cols(filter_pvalues[rep(1:nrow(filter_pvalues), filter_pvalues$n),1:3],
#           filter_pvalues$pvalue %>% unlist(recursive=F) %>% data_frame(pvalue=.))
# 
# # Results -----------------------------------------------------------------
# # Number of samples and features
# save(dims, file = "~/Dropbox/Presentation/data/dims.RData")
# save(dims, file = "data/dims.RData")
# save(edger_pvalues, submitter_pvalues, tabs, file = "data/pvalues.RData")

