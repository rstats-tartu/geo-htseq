
library(tidyverse)
library(GEOquery)

# Munge matrixfiles -------------------------------------------------------
local_suppfile_folder <- "/Volumes/Media/srp_example/data/counts" # "data/counts"
local_matrixfile_folder <- "/Volumes/Media/srp_example/data/matrix" # "data/matrix"

matrixfiles <- dir(local_matrixfile_folder) %>% 
  data_frame(series_matrix_file = .) %>% 
  mutate(Accession = str_extract(toupper(series_matrix_file), "GSE[[:digit:]]*")) %>% 
  select(Accession, everything())

suppfiles <- dir(local_suppfile_folder) %>% 
  data_frame(files = .) %>% 
  mutate(Accession = str_extract(toupper(files), "GSE[[:digit:]]*")) %>% 
  select(Accession, everything())

# Unique series matrix files
matrixfiles <- left_join(suppfiles, matrixfiles) %>% 
  select(-files) %>%
  distinct()

my_getGEO <- function(x, path = "data/matrix/") {
  message(x)
  path <- file.path(path, x)
  try(getGEO(filename = path, getGPL = FALSE))
}

gsem <- matrixfiles %>% 
  mutate(series_matrix = map(series_matrix_file, my_getGEO, path = local_matrixfile_folder))
saveRDS(gsem, file = "data/gsem.rds")
