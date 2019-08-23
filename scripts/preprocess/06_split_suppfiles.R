
# Load libs, local_suppfile_folder, and helper function geofile_df
source("scripts/_common.R")

munge_suppfiles <- function(suppfilenames_filtered, gsem, bad) {
  
  # supplementary files of interest expected to be locally available
  suppfilenames_filtered <- read_rds(suppfilenames_filtered)
  
  # Load series matrix data frames -----
  gsem <- read_rds(gsem)
  
  # Extract series matrixes
  gsem <- mutate(gsem, series_matrix = map(gse, "result"))
  
  # Filter Accessions with gsematrices present
  gsem <- filter(gsem, map_lgl(series_matrix, ~class(.x) == "ExpressionSet"))
  
  ## Read in local supplemental tables -----
  
  # File names of downloaded supplementary files
  supptabs <- geofile_df(local_suppfile_folder, "suppfiles")
  
  # Keep only files present in filtered suppfilenames
  supptabs <- filter(suppfilenames_filtered, type == "suppl") %>% 
    rename(suppfiles = files) %>% 
    inner_join(supptabs)
  
  # Merge gsem ExpressionSets to supptabs
  supptabs <- dplyr::select(gsem, Accession, series_matrix) %>% 
    nest(series_matrix, .key = "matrixfiles") %>% 
    left_join(supptabs, .)
  
  # Remove 'bad' files
  supptabs <- filter(supptabs, !(suppfiles %in% bad))
  
  mutate(supptabs, splits = rep_len(1:100, nrow(supptabs))) %>% 
    group_by(splits) %>% 
    nest() %>% 
    mutate(data = map2(data, splits, ~write_rds(.x, path = glue::glue("output/tmp/supptabs_{.y}.rds"))))
}

munge_suppfiles(suppfilenames_filtered = snakemake@input[["suppfilenames_filtered"]],
                gsem = snakemake@input[["gsem"]], bad = snakemake@params[["bad"]])
