# Load libs
source("R/_common.R")
library(entrezquery)

# Check if document summary table is already loaded
if (!"ds" %in% ls()) {
  ds <- readRDS("output/document_summaries.rds")
}

# Filter by data
ds_filtered <- ds %>% 
  filter(PDAT <= last_date)

# Start download
get_safely <- possibly(get_dirlist, data_frame)
start <- Sys.time()
suppfilenames <- ds_filtered %>%
  mutate(dirlist = map(Accession, ~ {message(.x); get_safely(.x)}))
end <- Sys.time()
cat(sprintf("Downloading filenames took %s hours\n", end-start), file = "suppfilename.log", append = T)

# Save downloaded dataset
saveRDS(suppfilenames, file = "output/suppfilenames.rds")
