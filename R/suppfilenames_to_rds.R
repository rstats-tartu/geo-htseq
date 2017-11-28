
source("lib/getDirListing.R")
load("data/suppfilenames_2017-11-27.RData")

suppfilenames <- suppfilenames %>% 
  mutate(suppfiles = map2(Accession, r, ~ {message(.x); try(get_dirlist(.y))}))
saveRDS(suppfilenames, "data/suppfilenames_2017-11-27.rds")
