
# Load libs
source("R/_common.R")
# Load out strings
source("R/out_strings.R")
pacman::p_load_gh("tpall/entrezquery")

filter_suppfilenames <- function(data_path, rds_path, txt_path) {
  
  # Load supplementary file names
  suppfilenames <- read_rds(data_path)
  
  # Munge file list
  sfn <- suppfilenames %>% 
    group_by(Accession) %>%
    mutate(result = map(dirlist, "result"),
           files = map(result, "file"),
           type = map(result, "type")) %>% 
    select(Accession, files, type) %>% 
    filter(map_chr(files, class) != "NULL") %>% 
    unnest()
  
 
  
  # Filter files of interest
  sfn <- sfn %>%
    filter(!str_detect(tolower(files), str_c(out_string1, collapse = "|")),
           !str_detect(tolower(files), str_c(out_string2, "(\\.gz|\\.bz2)?$", 
                                             collapse = "|")))
  
  write_rds(sfn, rds_path)
  
  ## Save text file with file names
  sfn %>% 
    mutate(files = file.path(type, files)) %>% 
    pull(files) %>% 
    unique() %>% 
    str_c(collapse = " ") %>% 
    write_lines(path = txt_path)
}

filter_suppfilenames(snakemake@input[[1]], snakemake@output[[1]], snakemake@output[[2]])

