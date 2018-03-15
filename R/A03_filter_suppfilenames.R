
# Load libs
message("Loading libs\n")
source("R/_common.R")
library(entrezquery)

# Load supplementary file names
message("Loading supplementary file names\n")
suppfilenames <- readRDS("output/suppfilenames.rds")

# Munge file list
message("Munging file list\n")
sfn <- suppfilenames %>% 
  group_by(Accession) %>%
  mutate(result = map(dirlist, "result"),
         files = map(result, "file"),
         type = map(result, "type")) %>% 
  select(Accession, files) %>% 
  filter(map_chr(files, class) != "NULL") %>% 
  unnest()

# Load out strings
message("Loading out strings\n")
source("R/out_strings.R")

# Filter files of interest
message("Filtering files of interest\n")
sfn <- sfn %>%
  filter(!str_detect(tolower(files), str_c(out_string1, collapse = "|")),
         !str_detect(tolower(files), str_c(out_string2, "(\\.gz|\\.bz2)?$", 
                                                   collapse = "|"))) 

write_rds(sfn, "output/suppfilenames_filtered.rds")

## save text file with file names
sfn %>% 
  mutate(files = case_when(
    str_detect(files, "soft.gz$") ~ file.path("soft", files),
    TRUE ~ file.path("suppl", files)
  )) %>% 
  pull(files) %>% 
  unique() %>% 
  str_c(collapse = " ") %>% 
  write_lines(path = "output/suppfilenames_filtered.txt")
