
pacman::p_load(tidyverse)
pacman::p_load_gh(entrezquery)

download_suppfiles <- function(data_path) {
  
  sfn <- read_rds(data_path)
  
  # Filter out local files
  local_matrix_files <- list.files("output/matrix")
  local_suppl_files <- list.files("output/suppl")
  
  sfn_missing <- sfn %>% 
    filter(!(files %in% c(local_suppl_files, local_matrix_files)))
  
  # Safe wrap download function
  safe_download <- safely(~ download_gsefile(.x, dest = "output"))
  
  # Start download
  sfn_missing %>% 
    ungroup() %>% 
    mutate(d = map(files, safe_download))
}

download_suppfiles(snakemake@input[[1]])
