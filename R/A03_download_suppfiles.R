
# Load libs
source("R/_common.R")
library(entrezquery)

# Load supplementary file names
suppfilenames <- readRDS("output/suppfilenames.rds")

# Safe wrap download function 
safe_download <- safely(~ download_gsefiles(.x, dest = "output"))

# Munge file list
sfn <- suppfilenames %>% 
  group_by(Accession) %>%
  mutate(result = map(dirlist, "result"),
         files = map(result, "file"),
         type = map(result, "type")) %>% 
  select(Accession, files) %>% 
  filter(map_chr(files, class) != "NULL") %>% 
  unnest()

# Load 
source("R/out_strings.R")

start <- Sys.time()
sfn_filtered <- sfn %>%
  filter(!str_detect(tolower(files), str_c(out_string1, collapse = "|")),
         !str_detect(tolower(files), str_c(out_string2, "(\\.gz|\\.bz2)?$", 
                                                   collapse = "|"))) 

sfn_filtered %>%
  group_by(Accession) %>%
  do(safe_download(.$files))
end <- Sys.time()

# Send remainder
msg <- sprintf("Hi Taavi!\nSupplementary files download took %s hours and ended 
               at %s.\nBest regards,\nYour Computer.", 
               end - start, Sys.time())
cmd <- sprintf("echo '%s' | mail -s 'Downloading files finished!' tapa741@gmail.com", msg)
system(cmd)
