
pacman::p_load(tidyverse)

merge_suppdata <- function(rds_files, out_path) {
  
  # Check for missing files
  imported <- parse_number(str_extract(unlist(rds_files), "[:digit:]+"))
  probs <- setdiff(1:100, imported)
  
  if(length(probs) > 0) {#"Warning! Missing files"
  
    logfiles <- list.files("output/slurm_log/", full.names = TRUE)
    logfiles <- logfiles[parse_number(str_extract(logfiles, "_[:digit:]+.out")) %in% probs] %>% 
      map(read_lines)
  
    return_bad_boy <- function(x) {
      last_file <- max(which(map_lgl(x, str_detect, "^output")))
      error <- (last_file + 1):(last_file + 2)
      data_frame(last_file = x[last_file], error = x[error])
    }
  
    errors <- data_frame(logfiles) %>% 
      mutate(error = map(logfiles, return_bad_boy)) %>% 
      unnest(error) %>% 
      filter(map_lgl(error, ~str_length(.x) != 0)) %>% 
      mutate(error = case_when(
        str_detect(error, "caught") ~ "segfault",
        str_detect(error, "slurmstepd") ~ "memory"
      )) %>% 
      distinct() 
  
    filter(errors, error == "segfault") %>% 
      mutate(last_file = str_remove(last_file, "output/suppl/")) %>% 
      pull(last_file)
    filter(errors, error == "memory")
  }
  
  if (length(probs) == 0) {
    suppdata_split <- map(rds_files, read_rds)
    suppdata <- bind_rows(suppdata_split)
    write_rds(suppdata, out_path)
  }
}

merge_suppdata(snakemake@input, snakemake@output[[1]])

