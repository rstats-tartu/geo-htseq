
# Load libs
source("R/_common.R")
library(entrezquery)

download_suppfilenames <- function(data_path, out_path, last_date) {
  
  # Check if document summary table is already loaded
  ds <- read_rds(data_path)
  
  # Filter by data
  ds_filtered <- filter(ds, PDAT <= ymd(last_date))
  
  # Start download supplementary and matrix files
  get_safely <- safely(~ get_dirlist(.x, dir = c("suppl", "matrix")))
  
  suppfilenames <- ds_filtered %>%
    mutate(dirlist = map(Accession, ~ {message(.x); get_safely(.x)}))
  
  # Save downloaded dataset
  write_rds(suppfilenames, path = out_path)
}

download_suppfilenames(snakemake@input[[1]], snakemake@output[[1]], last_date = snakemake@params[["last_date"]])

