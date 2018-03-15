
source("R/_common.R")
library(GEOquery)

# Munge matrixfiles -------------------------------------------------------
# local_matrixfile_folder and local_suppfile_folder come from _common.R
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

my_getGEO <- function(x, path) {
  message(x)
  path <- file.path(path, x)
  try(getGEO(filename = path, getGPL = FALSE))
}

gsem <- matrixfiles %>%
  mutate(gse = map(series_matrix_file, my_getGEO, path = local_matrixfile_folder))
write_rds(gsem, path = "output/gsem.rds")
