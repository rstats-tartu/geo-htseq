
source("R/_common.R")
library(GEOquery)

munge_matrixfiles <- function(out_path) {
  # Munge matrixfiles -------------------------------------------------------
  # local_matrixfile_folder and local_suppfile_folder come from _common.R
  matrixfiles <- local_matrixfile_folder %>% 
    dir(pattern = "GSE[[:digit:]]*", full.names = TRUE) %>% 
    data_frame(series_matrix_file = .) %>% 
    mutate(Accession = str_extract(toupper(series_matrix_file), "GSE[[:digit:]]*")) %>% 
    select(Accession, everything())
  
  suppfiles <- local_suppfile_folder %>% 
    dir(pattern = "GSE[[:digit:]]*", full.names = TRUE) %>% 
    data_frame(files = .) %>% 
    mutate(Accession = str_extract(toupper(files), "GSE[[:digit:]]*")) %>% 
    select(Accession, everything())
  
  # Unique series matrix files
  matrixfiles <- left_join(suppfiles, matrixfiles) %>% 
    select(-files) %>%
    distinct()
  
  safely_getGEO <- safely(
    function(path) {
      message(path)
      getGEO(filename = path, getGPL = FALSE)
    }
  )
  
  weird_stuff <- matrixfiles %>% 
    filter(!str_detect(series_matrix_file, "GSE[0-9]*(-GPL[0-9]*)?_series_matrix.txt.gz"))
  
  gsem <- matrixfiles %>%
    filter(str_detect(series_matrix_file, "GSE[0-9]*(-GPL[0-9]*)?_series_matrix.txt.gz")) %>% 
    mutate(gse = map(series_matrix_file, safely_getGEO, path = local_matrixfile_folder))
  
  write_rds(gsem, path = out_path)
}

munge_matrixfiles(snakemake@output[[1]])
