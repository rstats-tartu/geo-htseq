
library(dplyr)
library(purrr)
library(stringr)
library(GEOquery)

# Munge matrixfiles -------------------------------------------------------
local_suppfile_folder <- "/Volumes/Media/srp_example/data/counts" # "data/counts"
local_matrixfile_folder <- "/Volumes/Media/srp_example/data/matrix" # "data/matrix"

matrixfiles <- list.files(local_matrixfile_folder)
matrixfiles <- data_frame(matrixfiles) %>% 
  mutate(Accession = str_extract(toupper(matrixfiles), "GSE[[:digit:]]*"))

suppfiles <- list.files(local_suppfile_folder)
suppfiles <- data_frame(suppfiles) %>% 
  mutate(Accession = str_extract(toupper(suppfiles), "GSE[[:digit:]]*"))

my_getGEO <- function(x, path = "data/matrix/") {
  message(x)
  try(getGEO(filename = file.path(path, x), getGPL = FALSE))
}

# Unique series matrix files
matrixfiles <- left_join(suppfiles, matrixfiles) %>% 
  select(-suppfiles) %>%
  distinct()

gsem <- mutate(matrixfiles, gsematrix = map(matrixfiles, my_getGEO, path = local_matrixfile_folder))
save(gsem, file = "data/gsem.RData")
