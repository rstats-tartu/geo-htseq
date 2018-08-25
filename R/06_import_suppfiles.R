
# Load libs, settings and functions
source("R/_common.R")
source("R/munge_geo.R")
source("R/checkFullRank.R")
source("R/text_funs.R")
pacman::p_load(digest, glue)
pacman::p_load_gh("seandavi/GEOquery")

import_suppdata <- function(supptab, out_path) {
  supptab <- read_rds(supptab)
  supptab %>%
    mutate(result = map(suppfiles, ~ try(munge_geo_pvalue(file.path(local_suppfile_folder, .x))))) %>% 
    write_rds(path = out_path)
}

import_suppdata(snakemake@input[[1]], snakemake@output[[1]])

