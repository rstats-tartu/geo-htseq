
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
get_safely <- safely(get_dirlist)

start <- Sys.time()
suppfilenames <- ds_filtered %>%
  mutate(dirlist = map(Accession, ~ {message(.x); get_safely(.x)}))
end <- Sys.time()

# Save downloaded dataset
saveRDS(suppfilenames, file = "output/suppfilenames.rds")

# Send remainder
msg <- sprintf("Hi!\nDownload took %s hours and ended at %s.\nBest regards,\nYour Computer.", 
               end - start, Sys.time())
cmd <- sprintf("echo '%s' | mail -s 'Downloading filenames finished!' tapa741@gmail.com", msg)
system(cmd)

