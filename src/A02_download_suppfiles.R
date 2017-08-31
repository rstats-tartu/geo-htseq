
freshgeoquery <- FALSE
if(freshgeoquery && !file.exists("data/ds.RData")){
  source("src/A01_GEO_query.R")
} else{
  load("data/ds.RData")
}

source("lib/getDirListing.R")
library(tidyverse)
library(lubridate)
library(magrittr)
library(formattable)
library(stringr)
library(GEOquery)

# Retrieve filenames 
fullsuppfilenamesquery <- FALSE
if(fullsuppfilenamesquery){
  # getting directory listings for complete dataset may take long time: ~7 hours
  suppfilenames <- ds %>%
    filter(ymd(PDAT) <= last_date, 
           str_detect(taxon, "Mus musculus|Homo sapiens")) %>% 
    mutate(SuppFileNames = map(FTPLink, ~try(getDirListing(file.path(.x, "suppl/")))))
  save(suppfilenames, file = "data/suppfilenames_2017-06-19.RData")
} else {
  load("data/suppfilenames_2017-06-19.RData")
}

# Accessions from query missing from suppfilenames data.frame
missing_suppfilenames <- filter(ds, ymd(PDAT) <= last_date, 
                           str_detect(taxon, "Mus musculus|Homo sapiens"),
                           !(Accession %in% suppfilenames$Accession))
missing_suppfilenames <- mutate(missing_suppfilenames, SuppFileNames = map(FTPLink, ~try(getDirListing(file.path(.x, "suppl/")))))
suppfilenames <- bind_rows(suppfilenames, missing_suppfilenames)

# Files with name extentions that we are NOT going to download and --------

## @knitr out_strings
out_string1 <- c("filelist|annotation|readme|error|raw.tar|csfasta|bam|sam|bed|[:punct:]hic|hdf5|bismark|map|barcode|peaks")
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

## @knitr suppfile-download-code-n
# Folder where downloaded supplementary files are stored
local_suppfile_folder <- "/Volumes/Media/srp_example/data" # "data/counts"

# Download if suppfile is missing, some big files are moved into separate folder, hence recursive = T
local_files <- list.files(local_suppfile_folder, recursive = T) %>% 
  str_replace(".*\\/","")

missing_local <- suppfiles_of_interest %>% 
  filter(!str_replace(tolower(SuppFileNames),"\\.(gz|bz2)$","") %in% 
           str_replace(tolower(local_files), "\\.(gz|bz2)$",""))

# Download supplementary files --------------------------------------------
mutate(missing_local, map2(FTPLink, SuppFileNames, ~ geo_supp_dwnl(.x, .y, local_suppfile_folder)))

## Some readme files still passed the filter, let's remove those ----------
slippedin_readmes <- list.files(local_suppfile_folder, full.names = T)[str_detect(tolower(list.files(local_suppfile_folder)), "readme")]
system(sprintf("rm %s", paste(slippedin_readmes, collapse = " ")))
