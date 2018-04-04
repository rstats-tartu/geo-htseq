
library(tidyverse)

rds_files <- list.files("output/tmp", full.names = TRUE)
suppdata_split <- map(rds_files, read_rds)
suppdata <- bind_rows(suppdata_split)
write_rds(suppdata, "output/suppdata.rds")
