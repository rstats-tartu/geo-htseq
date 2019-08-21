
source("scripts/_common.R")
library(GEOquery)

munge_matrixfiles <- function(out_path) {
  # Munge matrixfiles -------------------------------------------------------
  # local_matrixfile_folder and local_suppfile_folder come from _common.R
  acc <- "GSE[[:digit:]]+"
  matrixfiles <- local_matrixfile_folder %>% 
    dir(pattern = acc, full.names = TRUE) %>% 
    data_frame(series_matrix_file = .) %>% 
    mutate(Accession = str_extract(toupper(series_matrix_file), acc)) %>% 
    select(Accession, everything())
  
  suppfiles <- local_suppfile_folder %>% 
    dir(pattern = acc, full.names = TRUE) %>% 
    data_frame(files = .) %>% 
    mutate(Accession = str_extract(toupper(files), acc)) %>% 
    select(Accession, everything())
  
  # Unique series matrix files
  matrixfiles <- left_join(suppfiles, matrixfiles) %>% 
    select(-files) %>%
    distinct()
  
  safely_getGEO <- safely(
    function(path) {
      message(path)
      stopifnot(file.size(path) > 0)
      getGEO(filename = path, getGPL = FALSE)
    }
  )
  gsm <- str_c(acc, "(-GPL[[:digit:]]+)?_series_matrix.txt.gz")
  weird_stuff <- matrixfiles %>% 
    filter(!str_detect(basename(series_matrix_file), gsm))
  
  gsem <- matrixfiles %>%
    filter(str_detect(basename(series_matrix_file), gsm)) %>% 
    mutate(gse = map(series_matrix_file, safely_getGEO))

  write_csv(weird_stuff, path = file.path(dirname(out_path), "not_real_series_matrix.csv"))
  write_rds(gsem, path = out_path)
}

munge_matrixfiles(snakemake@output[[1]])
