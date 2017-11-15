
library(dplyr)
library(lubridate)
library(GEOquery)

newquery <- FALSE
if(newquery) {
  source("src/A01_GEO_query.R")
} else {
  load("data/ds.RData")
}

# GEOquery ----------------------------------------------------------------

# Folder where downloaded supplementary files are stored
local_suppfile_folder <- "/Volumes/Media/srp_example/data/counts" # "data/counts"
local_matrixfile_folder <- "/Volumes/Media/srp_example/data/matrix" # "data/counts"

# Downloaded series matrix files
matrixfiles_acc <- data_frame(matrixfile = list.files(local_matrixfile_folder)) %>% 
  mutate(Accession = str_extract(matrixfile, "GSE[[:digit:]]*"))

# Downloaded supplementary files
suppfiles_acc <- data_frame(SuppFileNames = list.files(local_suppfile_folder)) %>% 
  mutate(Accession = str_extract(SuppFileNames, "GSE[[:digit:]]*")) %>% 
  nest(SuppFileNames)

# Missing series matrix files
gsematrix <- ds %>% 
  mutate_at("PDAT", ymd) %>% 
  filter(PDAT <= ymd("2017-06-19"),
         str_detect(taxon, "Mus musculus|Homo sapiens")) %>%
  select(Accession, PDAT) %>% 
  distinct %>% 
  anti_join(matrixfiles_acc) 

# Download supplemental file GEOs missing series matrix file
gsematrix <- mutate(gsematrix, 
                    gsem = map(Accession, ~ try(getGEO(GEO = .x, 
                                                       destdir = local_matrixfile_folder, 
                                                       getGPL = FALSE))))
