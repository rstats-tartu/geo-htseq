
library(tidyverse)
library(stringr)
library(magrittr)
library(GEOquery)

# Munge matrixfiles -------------------------------------------------------
local_suppfile_folder <- "/Volumes/Media/srp_example/data/counts" # "data/counts"
local_matrixfile_folder <- "/Volumes/Media/srp_example/data/matrix" # "data/matrix"

matrixfiles <- list.files(local_matrixfile_folder)
matrixfiles <- data_frame(matrixfiles) %>% 
  mutate(Accession = str_extract(toupper(matrixfiles), "GSE[[:digit:]]*"))

countfiles <- list.files(local_suppfile_folder)
countfiles <- data_frame(countfiles) %>% 
  mutate(Accession = str_extract(toupper(countfiles), "GSE[[:digit:]]*"))

my_getGEO <- function(x, path = "data/matrix/") {
  message(x)
  try(getGEO(filename = file.path(path, x), getGPL = FALSE))
}

# Unique series matrix files
matrixfiles <- left_join(countfiles, matrixfiles) %>% 
  select(-countfiles) %>%
  distinct()

gsem <- mutate(matrixfiles, gsematrix = map(matrixfiles, my_getGEO, path = local_matrixfile_folder))
gsem <- filter(gsem, map_lgl(gsematrix, function(x) class(x)=="ExpressionSet")) 

gsem <- mutate(gsem, samples = map_int(gsematrix, ~nrow(pData(.x))))
save(gsem, file = "data/gsem.RData")
