
# Load libs
source("R/_common.R")
# Load out strings
source("R/out_strings.R")

filter_suppfilenames <- function(suppfilenames, suppfilenames_filtered_out) {
  
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
                                             collapse = "|")))
  
  write_rds(suppfilenames_filtered, suppfilenames_filtered_out)
}

filter_suppfilenames(snakemake@input[["suppfilenames"]], 
                     snakemake@output[["suppfilenames_filtered_out"]])

