
source("src/A01_GEO_query.R")
source("lib/getDirListing.R")
library(tidyverse)
library(magrittr)
library(stringr)
library(GEOquery)

# Retrieve filenames 
start <- Sys.time()
suppfilenames <- ds %>%
  mutate(SuppFileNames = map(FTPLink, ~try(getDirListing(file.path(.x, "suppl/")))))
end <- Sys.time()
save(suppfilenames, file = "data/suppfilenames_2017-06-19.RData")
start-end 
# Time difference of -7 hours


# Files with name extentions that we are NOT going to download and --------

## @knitr out_strings
out_string1 <- c("filelist|annotation|readme|error|raw.tar|csfasta|bam|sam|bed|[:punct:]hic|hdf5|bismark")
out_string2 <- c("tar","gtf","(big)?bed(\\.txt|12|graph|pk)?","bw","wig","hic","gct(x)?","tdf",
         "gff(3)?","pdf","png","zip","sif","narrowpeak","fa")
out_strings <- unique(c(unlist(str_split(out_string1, "\\|")), out_string2))

# Potentially interesting supplemental files ------------------------------

## @knitr interesting-supplemental-files
suppfiles_of_interest <- suppfilenames %>% 
  unnest(SuppFileNames) %>%
  filter(!str_detect(tolower(SuppFileNames), out_string1),
         !str_detect(tolower(SuppFileNames), paste0(out_string2, "(\\.gz|\\.bz2)?$", collapse = "|"))) %>% 
  select(Accession, SuppFileNames, FTPLink, PDAT) %>% 
  mutate(filext = str_extract(tolower(SuppFileNames), "\\.[:alpha:]+([:punct:][bgz2]+)?$")) 

suppf_oi_perc <- percent(1-(n_distinct(suppfiles_of_interest$Accession)/n_distinct(suppfilenames$Accession)))

## @knitr suppfile-download-code-n
# Folder where downloaded supplementary files are stored
local_suppfile_folder <- "/Volumes/Media/srp_example/data" # "data/counts"

# Download if suppfile is missing, some big files are moved into separate folder
local_files <- list.files(local_suppfile_folder, recursive = T) %>% 
  str_replace(".*\\/","")

data_frame(local_files) %>% 
  filter(str_detect(tolower(local_files), "peak")) %$% 
  local_files

missing_local <- suppfiles_of_interest %>% 
  filter(!str_replace(tolower(SuppFileNames),"\\.(gz|bz2)$","") %in% 
           str_replace(tolower(local_files), "\\.(gz|bz2)$",""))

# Download supplementary files --------------------------------------------
missing_local %>% 
  mutate(map2(FTPLink, SuppFileNames, ~ geo_supp_dwnl(.x, .y, local_suppfile_folder)))
