
library(tidyverse)
sfn <- read_rds("output/suppfilenames_filtered.rds")

# Safe wrap download function
safe_download <- safely(~ download_gsefile(.x, dest = "output"))

# Start download
message("Downloading files\n")
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