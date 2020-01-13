
library(tidyverse)
library(entrezquery)

download_suppfiles <- function(suppfilenames_filtered) {
  
  suppfilenames_filtered <- read_rds(suppfilenames_filtered)
  
  # Filter out local files
  local_matrix_files <- list.files("output/matrix")
  local_suppl_files <- list.files("output/suppl")
  
  suppfilenames_tobe_downloaded <- filter(suppfilenames_filtered, !(files %in% c(local_suppl_files, local_matrix_files)))
  
  # Safe wrap download function
  safe_download <- safely(~ download_gsefile(.x, dest = "output"))
  
  # Start download
  mutate(suppfilenames_tobe_downloaded, d = map(files, safe_download))
}

download_suppfiles(snakemake@input[[1]])
