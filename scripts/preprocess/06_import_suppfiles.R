
# Load libs, local_suppfile_folder, and helper functions
source("scripts/_common.R")
source("scripts/munge_geo.R")

library(digest)
library(glue)
library(GEOquery)

import_suppdata <- function(supptab, outpath, bad) {
  supptabs <- read_rds(supptab)
  # Remove 'bad' files
  supptabs <- filter(supptabs, !(suppfiles %in% bad))
  supptabs %>%
    mutate(result = map(suppfiles, ~ try(munge_geo_pvalue(file.path(local_suppfile_folder, .x))))) %>% 
    write_rds(path = outpath)
}

import_suppdata(supptab = snakemake@input[[1]], 
                outpath = snakemake@output[[1]],
                bad = snakemake@params[["bad"]])

