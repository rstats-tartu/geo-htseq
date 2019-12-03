library(tidyverse)
library(here)

run_dir <- "output/read_runs/"
run_files <- list.files(run_dir)
file_path <- str_c(run_dir, run_files)
runs <- tibble(file_path) %>% 
  mutate(file_size = file.size(file_path)) %>% 
  filter(file_size > 0)


summarise_reads <- function(path){
  data <- read_tsv(path, col_types = "ccciccccccddDccc")
  data %>% 
    select(-ends_with("accession"), -experiment_alias, -sample_title) %>% 
    select(Accession = study_alias, everything()) %>% 
    group_by(Accession, tax_id, scientific_name, instrument_model, library_layout, library_strategy, library_source, library_selection, PDAT = first_public) %>% 
    summarise_at(vars(read_count, base_count), list(mean = mean, median = median, sd = sd, min = min, max = max, n = length)) %>% 
    select(-base_count_n) %>% 
    rename(n = read_count_n)
}

safely_summarise_reads <- safely(summarise_reads)
runs_sum <- runs %>% 
  mutate(sum = map(file_path, safely_summarise_reads))

transposed <- runs_sum %>% 
  pull(sum) %>% 
  transpose()
read_runs <- transposed$result %>% 
  bind_rows()
read_runs %>% 
  write_csv(here("output/read_runs.csv"))
