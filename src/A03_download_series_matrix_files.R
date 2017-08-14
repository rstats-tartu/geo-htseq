
library(tidyverse)
library(stringr)
library(GEOquery)

source("src/A01_GEO_query.R")

# GEOquery ----------------------------------------------------------------

# Folder where downloaded supplementary files are stored
local_suppfile_folder <- "/Volumes/Media/srp_example/data/counts" # "data/counts"

# Downloaded series matrix files
matrixfiles_acc <- data_frame(matrixfile = list.files("data/matrix/")) %>% 
  mutate(Accession = str_extract(matrixfile, "GSE[[:digit:]]*"))

# Downloaded supplementary files
suppfiles_acc <- data_frame(SuppFileNames = list.files(local_suppfile_folder)) %>% 
  mutate(Accession = str_extract(SuppFileNames, "GSE[[:digit:]]*")) %>% 
  nest(SuppFileNames)

# Download missing series matrix files
gsematrix <- ds %>% 
  filter(PDAT<="2017/06/22") %>%
  select(Accession, PDAT) %>% 
  distinct %>% 
  # left_join(suppfiles_acc,.) %>% # keep only GEOs with supplemental files
  anti_join(matrixfiles_acc) %>% # supplemental file GEOs missing series matrix file
  mutate(gsem = map(Accession, ~ try(getGEO(GEO = .x, destdir = "data/matrix/", getGPL = FALSE))))

# Missing series matrix files are behind password
missing_series_matrix <- gsematrix %>% 
  filter(map_lgl(gsem, ~inherits(.x, "try-error")))

library(lubridate)
missing_series_matrix %>% 
  filter(!is.na(PDAT)) %>% 
  count(PDAT) %>% 
  mutate(n = cumsum(n),
         PDAT = ymd(PDAT)) %>% 
  ggplot(aes(PDAT, n, group = 1)) + 
  geom_line() +
  ggtitle("GEO series with missing series matrix file") +
  ylab("Number of GEO series")

