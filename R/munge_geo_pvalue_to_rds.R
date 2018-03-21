
library(tidyverse)
library(Biobase)
source("lib/munge_geo.R")

suppl <- list.files("output/suppl", full.names = TRUE)

safe_munge_geo_pvalue <- safely(munge_geo_pvalue)
munge_geo_pvalue_wrap <- function(x) {
  filename <- str_c(basename(x), ".rds")
  path <- file.path("output/suppl_rds", filename)
  if (file.exists(path)) stop("File exists!")
  x <- safe_munge_geo_pvalue(x)
  write_rds(x, path)
}

map(suppl, safely(munge_geo_pvalue_wrap))
