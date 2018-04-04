
# Load libs
source("R/_common.R")
library(entrezquery)

# Load supplementary file names
suppfilenames <- read_rds("output/suppfilenames.rds")

# Munge file list
sfn <- suppfilenames %>% 
  group_by(Accession) %>%
  mutate(result = map(dirlist, "result"),
         files = map(result, "file"),
         type = map(result, "type")) %>% 
  select(Accession, files, type) %>% 
  filter(map_chr(files, class) != "NULL") %>% 
  unnest()

# Load out strings
source("R/out_strings.R")

# Filter files of interest
sfn <- sfn %>%
  filter(!str_detect(tolower(files), str_c(out_string1, collapse = "|")),
         !str_detect(tolower(files), str_c(out_string2, "(\\.gz|\\.bz2)?$", 
                                                   collapse = "|")))

write_rds(sfn, "output/suppfilenames_filtered.rds")

## Save text file with file names
sfn %>% 
  mutate(files = file.path(type, files)) %>% 
  pull(files) %>% 
  unique() %>% 
  str_c(collapse = " ") %>% 
  write_lines(path = "output/suppfilenames_filtered.txt")
