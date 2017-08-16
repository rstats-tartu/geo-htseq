
# Load libs ---------------------------------------------------------------

library(tidyverse)
library(stringr)
library(GEOquery)
library(magrittr)
library(ggplot2)
nrowthreshold <- 8000

# GEO query results and document summaries --------------------------------
# source("src/A01_GEO_query.R")
load("data/ds.RData")

# Produces gsem -- series matrix data frames
# local_suppfile_folder <- "/Volumes/Media/srp_example/data/counts" # "data/counts"
update_geoseriesmatrix_files <- FALSE
if(update_geoseriesmatrix_files){
  source("src/local_matrixfiles.R")
} else {
  load("data/gsem.RData")
  load("data/supptabs.RData")
}

# Read in local supplemental tables ---------------------------------------
source("lib/munge_geo.R")
source("lib/read_tabs.R") # delim import fun
source('lib/text_funs.R')
source("lib/checkFullRank.R")


# File name extensions
file_ext <- supptabs %>% 
  mutate(filext = str_replace(tolower(countfiles), "\\.(gz|bz2)$", "") %>% str_extract("\\.[_A-z0-9]*$")) %>% 
  group_by(filext) %>% 
  summarise(N = n()) %>%
  ungroup() %>% 
  mutate(perc = (N/sum(N))*100) %>% 
  arrange(desc(N))

file_ext_plot <- ggplot(file_ext, aes(filext, perc)) +
  geom_point() +
  scale_x_discrete(limits = file_ext$filext) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("Percent of files") +
  xlab("Supplementary file extention")

# Duplicated files
supptabs_duplicated <- supptabs %>% 
  mutate(files = str_replace(countfiles, "\\.(gz|bz2)$","")) %>% 
  filter(duplicated(files))

# Files that may crash R session
bad <- c("GSE93374_Merged_all_020816_DGE.txt.gz", 
         "GSE88931_RNA-seq_MergedReadCounts.tsv.gz",
         "GSE55385_transcripts_GSE.tsv.gz",
         "GSE74549_ChIP_1kbWindows_correctedReadCount.txt.gz",
         "GSE77213_Nguyen_GEO_TN03_tallies.xls",
         "GSE77213_Nguyen_GEO_TN05_tallies_total.xls",
         "GSE60012_100bpTiles_RRBS_Mouse.txt.gz")

local_suppfile_folder <- "/Volumes/Media/srp_example/data/counts/"
start <- Sys.time()
st <- supptabs %>% 
  left_join(select(ds, Accession, PDAT)) %>% 
  filter(!countfiles %in% bad,
         !str_detect(tolower(countfiles),"\\.r([dat]*)?(\\.gz)?$")) %>% 
  mutate(result = map2(data, countfiles, ~ try(munge_geo(.x, .y, path = local_suppfile_folder))))
end <- Sys.time()
start-end
# Time difference of -27.8913 mins
system("say the job is completed, boss")
save(st, file = "data/st.RData")
# load("data/st.RData")

# Add sheet names to excel tables
res <- mutate(st, sheetnames = map(result, ~if(!is.data.frame(.x)) {str_c("_", names(.x))} else {""}),
       excelfiles = map2(countfiles, sheetnames, str_c, sep = "")) %>% 
  select(-sheetnames)

dims <- res %>% 
  filter(map_lgl(result, is.data.frame)) %>%
  select(-data) %>% 
  unnest

dims <- left_join(dims, gsem) %>% 
  group_by(excelfiles) %>%
  mutate(idcols = columns - samples) %>% 
  filter(idcols>=0, idcols == min(idcols)) %>% 
  ungroup %>% 
  select(-matrixfiles) %>% 
  distinct

dims %>%
  select(-excelfiles, -idcols) %>% 
  filter(features > nrowthreshold) %>%
  summarise_at(vars(features, samples), funs("mean","median"))

# P values ----------------------------------------------------------------
tabs <- filter(res, !map_lgl(result, is.data.frame))
tabs <- tabs[rep(1:nrow(tabs), mutate(tabs, n = map_int(result, length))$n), 
             c("Accession", "countfiles", "PDAT", "excelfiles")] %>%
  bind_cols(data_frame(result = unlist(tabs$result, recursive = FALSE))) %>% 
  filter(!map_lgl(result, is.null))

pvalues <- filter(tabs, map_lgl(result, is.tibble))
pvalues %>% 
  filter(map_lgl(countfiles, str_detect, pattern = "xls"))

pvalues %>% 
  mutate(features = map_int(result, nrow)) %>% 
  ggplot(aes(features)) +
  geom_histogram(bins=30) +
  scale_x_log10() +
  geom_vline(xintercept = 8000, linetype = 3) +
  xlab("Number of P values")

# Fit models -----------------------------------------------------
# Filter ExpressionSets
esets <- tabs %>% filter(map_chr(result, ~class(.x)[1])=="ExpressionSet")

# Identify vars for model fitting
esets %<>% mutate(groups = map(result, ~get_model_eset(pData(.x))))

library(edgeR)
# GSE91063_Torrent_run_1_RPKM_cutoff.xlsx has only 3720 reads
esets %<>% 
  filter(map_lgl(result, ~nrow(exprs(.x))>nrowthreshold)) %>% 
  mutate(dge = pmap(list(result, groups, countfiles), function(x,y,z) {message(z); DGEList(counts=exprs(x), samples = y, genes = fData(x))})) %>% 
  filter(!map_lgl(groups, is.null))

# Make design formulas and model.matrix
make_formula <- . %>% colnames %>% paste(collapse="+") %>% paste("~",.) %>% formula()
esets %<>% mutate(design = map(groups, make_formula)) 
esets %<>% mutate(design = map2(design, groups, ~model.matrix(object=.x, data=.y))) 

# Filter genes with low counts: 3 or more samples with cpm of 10/L. L, # millions of counts in the smallest library
esets %<>% mutate(keep = map(dge, ~rowSums(cpm(.x)>10/min(.x$samples$lib.size)/1e6) >= 3),
                  dge_filt = map2(dge, keep, ~.x[.y,]))

# Recalculate library sizes and calculate normalisation factors
esets %<>% mutate(dge_filt = map(dge_filt, ~{.x$samples$lib.size <- colSums(.x$counts); .x}))
esets %<>% mutate(dge_filt = map2(countfiles, dge_filt, ~{message(.x); calcNormFactors(.y)}))
esets %<>% mutate(dge_filt = pmap(list(dge_filt, design, countfiles), 
                                  function(x,y,z) {message(z); try(estimateDisp(x,y))}))

# Fit models
esets %<>% 
  filter(!map_lgl(dge_filt, ~inherits(.x, "try-error"))) %>%
  mutate(fit = pmap(list(dge_filt, design, countfiles), function(x,y,z) {message(z); try(glmFit(x,y))})) %>% 
  filter(!map_lgl(fit, ~inherits(.x, "try-error"))) %>%
  mutate(results = map(fit, ~glmLRT(.x, coef=2:ncol(.x$design))),
         toptags = map(results, topTags, n=Inf))
save(esets, file="data/edgeR_fits_esets.RData")
# load("data/edgeR_fits_esets.RData")
# Extract edgeR p values

edger_pvalues <- esets %>% 
  mutate(tt = map(toptags, ~.x@.Data[[1]]),
         pvalue = map(tt, "PValue")) %>% 
  select(Accession, countfiles, PDAT, pvalue)

# Filter submitter provided p values
# pv <- pvalues %>% filter(countfiles=="GSE78078_Processed_file_cd38_negative_kd_oe.xlsx") %$% result[[1]]
basemean_filter <- function(pv, basemeanthresh=10) {
  
  colnames(pv) <- make.names(colnames(pv), unique = T)
  if(!all(is.na(pv$basemean))){
    pv <- filter(pv, basemean>basemeanthresh)
  }
  pv %>% 
    select(-basemean) %>% 
    filter(complete.cases(.)) %>% 
    as.list()
}

pvalues %>% filter(map_lgl(result, ~ncol(.x)>2))

filter_pvalues <- pvalues %>% 
  mutate(pvalue = map2(result, countfiles, ~{message(.y); basemean_filter(.x)}),
         n=map_int(pvalue, ~length(.x)))

submitter_pvalues <- bind_cols(filter_pvalues[rep(1:nrow(filter_pvalues), filter_pvalues$n),1:3],
          filter_pvalues$pvalue %>% unlist(recursive=F) %>% data_frame(pvalue=.))

# Results -----------------------------------------------------------------
# Number of samples and features
save(dims, file = "~/Dropbox/Presentation/data/dims.RData")
save(dims, file = "data/dims.RData")
save(edger_pvalues, submitter_pvalues, tabs, file = "data/pvalues.RData")

