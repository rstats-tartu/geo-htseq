
# Load libs
source("R/_common.R")
library(entrezquery)

suppfilenames <- readRDS("output/suppfilenames.rds")

safe_download <- safely(~ download_gsefiles(.x, dest = "output"))
out <- suppfilenames %>% 
  group_by(Accession) %>%
  mutate(result = map(dirlist, "result")),
         files = map(result, "file"))
out %>% 
  select(Accession, files) %>% 
  unnest() %>% 
  filter(!stringr::str_detect(files, "RAW.tar|filelist.txt")) %>% 
  group_by(Accession) %>%
  do(dwnl = safe_download(.$files))
