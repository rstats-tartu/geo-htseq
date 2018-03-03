
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

# Safe wrap download function
safe_download <- safely(~ download_gsefile(.x, dest = "output"))

# Start download
message("Sownloading files\n")
start <- Sys.time()
sfn %>% 
  ungroup() %>% 
  mutate(d = map(files, safe_download))
end <- Sys.time()

# Job finished, send email
message("Job finished, sending email\n")
msg <- sprintf("Subject:Downloading files finished!\n\nHi Taavi!\n\nSupplementary files download took %s and ended at %s.\n\nBest regards,\nYour Computer.", 
               format(difftime(end, start)), Sys.time())
cmd <- sprintf("echo -e '%s' | sendmail tapa741@gmail.com", msg)
system(cmd)
message("End")
