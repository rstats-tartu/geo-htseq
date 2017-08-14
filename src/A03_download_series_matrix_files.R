
library(tidyverse)
library(lubridate)
library(stringr)
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

# Download missing series matrix files
gsematrix <- ds %>% 
  mutate_at("PDAT", ymd) %>% 
  filter(PDAT<=ymd("2017-06-19")) %>%
  select(Accession, PDAT) %>% 
  distinct %>% 
  anti_join(matrixfiles_acc) 

# supplemental file GEOs missing series matrix file
gsematrix <- mutate(gsematrix, gsem = map(Accession, ~ try(getGEO(GEO = .x, 
                                                                  destdir = local_matrixfile_folder, 
                                                                  getGPL = FALSE))))

# Missing series matrix files are behind password
missing_series_matrix <- gsematrix %>% 
  filter(map_lgl(gsem, ~inherits(.x, "try-error")))

missing_series_matrix %>% 
  filter(!is.na(PDAT)) %>% 
  count(PDAT) %>% 
  mutate(n = cumsum(n),
         PDAT = ymd(PDAT)) %>% 
  ggplot(aes(PDAT, n, group = 1)) + 
  geom_line() +
  ggtitle("GEO series with missing series matrix file") +
  ylab("Number of GEO series")

