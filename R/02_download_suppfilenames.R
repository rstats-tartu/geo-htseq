
# Load libs
source("R/_common.R")
library(entrezquery)

# Check if document summary table is already loaded
ds <- readRDS("output/document_summaries.rds")

# Filter by data
ds_filtered <- ds %>% 
  filter(PDAT <= last_date)

# Start download supplementary and matrix files
get_safely <- safely(~ get_dirlist(.x, dir = c("suppl", "matrix")))

start <- Sys.time()
suppfilenames <- ds_filtered %>%
  mutate(dirlist = map(Accession, ~ {message(.x); get_safely(.x)}))
end <- Sys.time()

# Save downloaded dataset
saveRDS(suppfilenames, file = "output/suppfilenames.rds")
