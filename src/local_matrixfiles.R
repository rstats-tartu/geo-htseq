
library(tidyverse)
library(stringr)
library(magrittr)
library(GEOquery)

# Munge matrixfiles -------------------------------------------------------
local_suppfile_folder <- "/Volumes/Media/srp_example/data/counts" # "data/counts"
local_matrixfile_folder <- "/Volumes/Media/srp_example/data/matrix" # "data/matrix"
matrixfiles <- list.files(local_matrixfile_folder)
matrixfiles <- data_frame(matrixfiles) %>% 
  mutate(Accession = str_extract(matrixfiles, "GSE[[:digit:]]*"))

countfiles <- list.files(local_suppfile_folder)
countfiles <- data_frame(countfiles) %>% 
  mutate(Accession = str_extract(countfiles, "GSE[[:digit:]]*"))

# matrixfiles <- matrixfiles[str_extract(matrixfiles, "GSE[[:digit:]]*") %in% suppfiles_acc]

my_getGEO <- function(x, path = "data/matrix/") {
  message(x)
  try(getGEO(filename = file.path(path, x), getGPL = FALSE))
}

# Unique series matrix files
matrixfiles <- inner_join(countfiles, matrixfiles) %>% 
  select(-countfiles) %>%
  distinct()

gsem <- mutate(matrixfiles, gsematrix = map(matrixfiles, my_getGEO, path = local_matrixfile_folder))
gsem <- filter(gsem, map_lgl(gsematrix, function(x) class(x)=="ExpressionSet")) 
rm(matrixfiles)

# Expand platform data
platforms <- function(x){ 
  if(all(c("instrument_model", "library_strategy") %in% colnames(x))){
    grouped_df <- group_by(x, instrument_model, library_strategy) 
  } else {
    return(summarise(x, samples = n()))
  }
  summarise(grouped_df, samples = n())
}

gpl_gsem <- mutate(gsem, samples = map(gsematrix, ~platforms(pData(.x))))


# Download GPL files for platform info ------------------------------------

library(GEOquery)

# Get platform info out of GPL
local_softfile_folder <- "/Volumes/Media/srp_example/data/soft"
gpl <- data_frame(gpl = list.files(local_softfile_folder, full.names = T)) %>% 
  mutate(platform = map_chr(gpl, ~{message(.x); getGEO(filename = .x)@header$title}))
gpl <- separate(gpl, platform, c("instrument", "species"), sep = "\\(")
w <- warnings()

# Update missing species info manually
gpl[c(2, 98),]$species <- c("Homo sapiens", "Mus musculus")
gpl[c(69, 93),]$species
gpl <- mutate(gpl, gpl = map_chr(gpl, str_extract, pattern="GPL[0-9]*"),
                species = map_chr(species, ~str_replace(.x, "\\)",""))) %>% 
  mutate_at(vars(instrument, species), function(x) map_chr(x, str_trim))


# gpl %>% filter(instrument=="Illumina HiSeq 4000")
gsem <- mutate(gsem, gpl = map_chr(gsematrix, annotation)) %>% 
  left_join(gpl)

# Let's throw out platforms for other species than mouse and human
gsem <- filter(gsem, !is.na(instrument))

gsem <- mutate(gsem, samples = map_int(gsematrix, ~nrow(pData(.x))))
save(gsem, file = "data/gsem.RData")

# Read in supp files ------------------------------------------------------
source("lib/read_tabs.R") # delim import fun
source('lib/text_funs.R')
source("lib/checkFullRank.R")

supptabs <- inner_join(countfiles, gsem) %>% select(-matrixfiles)
supptabs <- select(supptabs, Accession, countfiles, gsematrix) %>% nest(gsematrix)
save(supptabs, file = "data/supptabs.RData")
