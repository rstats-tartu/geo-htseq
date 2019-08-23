
source("scripts/_common.R")

# Load out strings
source("scripts/out_strings.R")

filter_suppfilenames <- function(suppfilenames, outpath) {
  
  # Load supplementary file names
  suppfilenames <- read_rds(suppfilenames)
  
  # Munge file list
  suppfilenames_unnested <- suppfilenames %>% 
    group_by(Accession) %>%
    mutate(result = map(dirlist, "result"),
           files = map(result, "file"),
           type = map(result, "type")) %>% 
    select(Accession, files, type) %>% 
    filter(map_chr(files, class) != "NULL") %>% 
    unnest()
  
  # Filter files of interest
  suppfilenames_filtered <- suppfilenames_unnested %>%
    filter(!str_detect(tolower(files), str_c(out_string1, collapse = "|")),
           !str_detect(tolower(files), str_c(out_string2, "(\\.gz|\\.bz2)?$", 
                                             collapse = "|"))) %>% 
    ungroup()
  
  # more than six hundred series matrix files are downloaded from supplementary files folder, relabel them as matrix files
  suppfilenames_filtered  <- mutate(suppfilenames_filtered, 
                                    type = case_when(
                                      str_detect(files, "series_matrix.txt.gz$") ~ "matrix",
                                      TRUE ~ "suppl"
                                    ))
  
  write_rds(suppfilenames_filtered, outpath)
}

filter_suppfilenames(snakemake@input[[1]], 
                     snakemake@output[[1]])

