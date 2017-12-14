
library(tidyverse)
library(httr)
library(lubridate)
source("lib/getDirListing.R")

if(!"ds" %in% ls()){
  ds <- readRDS("data/document_summaries.rds")
}

ds_filtered <- filter(ds, 
                      ymd(PDAT) <= "2017-06-19", 
                      str_detect(taxon, "Mus musculus|Homo sapiens"))

start <- Sys.time()
suppfilenames <- mutate(ds_filtered, url = file.path(FTPLink, "suppl/")) %>% 
  sample_n(10) %>% 
  mutate(r = map(url, ~{message(.x); try(httr::GET(.x))}))
end <- Sys.time()
end-start
cat(sprintf("Downloading filenames took %s hours\n", end-start), file = "suppfilename.log", append = T)

saveRDS(suppfilenames, file = "data/suppfilenames_2017-11-27.rds")

