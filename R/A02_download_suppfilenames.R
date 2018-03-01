
# Load libs
source("R/_common.R")
library(entrezquery)

if (!"ds" %in% ls()) {
  ds <- readRDS("output/document_summaries.rds")
}

ds_filtered <- ds %>% 
  filter(PDAT <= last_date)

start <- Sys.time()
get_safely <- possibly(get_dirlist, data_frame)
ds_filtered <- ds_filtered %>%
  sample_n(10) %>% 
  mutate(dirlist = map(Accession, ~ {message(.x); get_safely(.x)}))
end <- Sys.time()
cat(sprintf("Downloading filenames took %s hours\n", end-start), file = "suppfilename.log", append = T)

safe_download <- possibly(~ download_gsefiles(.x, dest = "output"), "fail")
ds_filtered <- ds_filtered %>% 
  mutate(files = map(dirlist, "file"),
         download = map(files, ~ {message(paste(.x, collapse = "\n")); safe_download(.x)}))

saveRDS(suppfilenames, file = "output/suppfilenames.rds")

