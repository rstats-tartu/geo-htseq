###
#Counts how many pvalues where found
###

library(dplyr)

source("R/_common.R")

get_stat <- function(path,title) {
  ###
  #Gets data from single file and writes statistics to file
  ###
  
  print(path)
  ## ---- loadst -----
  st <- readRDS(path)
  
  # imported files
  imported_geos <- select(st, Accession) %>%  n_distinct()

  st_unnested <- unnest(st, result)

  st_unnested <- unnest(st_unnested, sheets)

  # Add sheet names to xls files
  st_unnested <- st_unnested %>% 
    mutate(suppdata_id = case_when(
      str_length(sheets) > 0 ~ str_c(suppfiles, "-sheet-", sheets),
      str_length(sheets) == 0 ~ suppfiles
    ))

  st_pvalues <- filter(st_unnested,!map_lgl(pvalues,is.null),
                       !map_lgl(pvalues, inherits, "try-error"))

  pvalue_geos <- select(st_pvalues, Accession) %>%  n_distinct()

  write(paste(title,imported_geos,pvalue_geos,nrow(st_unnested),nrow(st_pvalues),sep="\t"),
        "output.txt",append=TRUE)
  
}

    
write(paste("Source","Imported_GEOs","Pvalue_GEOs","Imported_datasets","Pvalue_datasets",sep="\t"),
      "output.txt") 

#get_stat(file.path("output","tmp_tabs_data","suppdata_1.rds"),"Ori")  
get_stat(file.path("output","Taavis_import_versio","suppdata.rds"),"Original")
get_stat(file.path("output","Taavis_import_versio","suppdata_regex.rds"),"Regex") 
get_stat(file.path("output","Taavis_import_versio","suppdata_regex_xls.rds"),"Regex_xls") 

