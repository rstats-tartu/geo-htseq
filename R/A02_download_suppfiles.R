
library(dplyr)
library(purrr)
library(tidyr)
library(httr)
library(lubridate)
library(stringr)
source("lib/getDirListing.R")

if(!"ds" %in% ls()){
  load("data/ds.RData")
}

ds_filtered <- filter(ds, 
                      ymd(PDAT) <= "2017-06-19", 
                      str_detect(taxon, "Mus musculus|Homo sapiens"))

start <- Sys.time()
suppfilenames <- mutate(ds_filtered, url = file.path(FTPLink, "suppl/")) %>% 
  sample_n(100) %>% 
  mutate(r = map(url, ~{message(.x); try(httr::GET(.x))}))
end <- Sys.time()
end-start
cat(sprintf("Downloading filenames took %s hours\n", end-start), file = "suppfilename.log", append = T)

save(suppfilenames, file = "data/suppfilenames_2017-08-25.RData")
suppfilenames  <- filter(suppfilenames, !map_lgl(r, inherits, "try-error")) %>% 
  mutate(suppfiles = map(r, get_dirlist))
