 
# Funs --------------------------------------------------------------------

get.delim_fs <- function(path) {
  delim <- try( reader::get.delim(path, n = 50, large = 50), silent = T)
  if(inherits(delim, "try-error") | length(delim) == 0){
    return("\t")
  }
  return(delim)
}


my_get_delim <- function(path){
  tab <- reader::n.readLines(path, 50)
  
  get_max_names <- function(x) names(x[x==max(x)])
  
  delims <- str_extract_all(tab, ",|;|[[:space:]]") 
  allempty <- delims %>% sapply(is_empty) %>% all
  if(allempty){
    return("")
    }
  
  delims %<>% 
    sapply(table, simplify = FALSE) %>%
    sapply(get_max_names) %>% unlist %>% 
    table %>% get_max_names 
  
  if(length(delims) > 1){
     ws <- str_detect(delims, " ")
     return(delims[!ws])
  }
  return(delims)
}

# read delimited tables ---------------------------------------------------

read_tabs <- function(path) {
  
  if(stringr::str_detect(path, "gct$")){
    tab <- CePa::read.gct(path)
    return(tab)
  }
  
  if(stringr::str_detect(path, "csv$")){
    
    tab <- read.csv(path, stringsAsFactors = FALSE)
    
    if(ncol(tab)>1){
      return(tab)
    }
  }
  
  delimiter <- my_get_delim(path)
  
  tab <- try(read.delim(path, sep = delimiter, stringsAsFactors = FALSE))

  if(all(is.na(as.numeric(tab[1,])))){
    tab <- try(read.delim(path, sep = delimiter, skip = 1, stringsAsFactors = FALSE))
    if(inherits(tab, "try-error")){
      return(tibble())
    }
    return(tab)
  }
  
  if(ncol(tab)==1) {
    tab <- read.delim(path, skip = 0, stringsAsFactors = FALSE)
    
    if(ncol(tab)==1){
      tab <- read.delim(path, skip = 1, stringsAsFactors = FALSE)
    }
    return(tab)
  }
  
  return(tab)
}


# read excel sheets -------------------------------------------------------

read_excel_fs <- function(path) {
  
  sheets <- try(readxl::excel_sheets(path), silent = TRUE)
  
  rex <- function(path, sheet){
    tab <- try(readxl::read_excel(path, sheet), silent = TRUE)
    if(inherits(tab, "try-error")){
      return(tibble())
    }
    return(tab)
  }
  
  if(!inherits(sheets, "try-error")){
    tablist <- lapply(sheets, function(x) rex(path, sheet = x))
    names(tablist) <- sheets
  } else {
    return(tibble())
  }
  
  if(length(tablist) == 1){
    tab <- tablist[[1]]
    return(tab)
  }
  
  return(tablist)
}

# Unnest tables from xlsfiles ---------------------------------------------
unnest_xls_sheets <- function(tib){
  suppt <- mutate(tib, islist = map_lgl(tables, is.vector))
  xs <- filter(suppt, islist) %>% mutate(n = map_int(tables, length)) 
  sfn <- transmute(xs, countfiles = map2(countfiles, n, rep)) %>% unnest
  tb <- xs$tables %>% unlist(recursive = F) 
  sfn$sheet <- names(tb) %>% rm_punct_tolower
  tb %<>% data_frame(tables = .)
  sfn %<>% mutate(countfiles = paste(countfiles, sheet, sep="_")) %>% select(countfiles) 
  xlsfiles <- bind_cols(sfn, tb)
  suppt %>% filter(!islist) %>% 
    bind_rows(xlsfiles) %>% 
    select(-islist)
}

# read all together -------------------------------------------------------

readGEOtabs <- function(path){
  
  message(path)
  system(paste("echo '", path, "' >> log.txt"))
  
  if(stringr::str_detect(path, "xlsx?")){
    tab <- read_excel_fs(path)
    return(tab)
  }
  
  if(stringr::str_detect(path, "gct$")){
    tab <- CePa::read.gct(path)
    return(tab)
  }
  
  if(stringr::str_detect(path, "\\.gz$")){
    path <- paste0("zcat < ", path)
  }
  
  tab <- try(data.table::fread(path, data.table = F))
  
  if(inherits(tab, "try-error")){
    message <- tab[1]
    system(paste("echo '", path, "\n", message,"' >> log.txt"))
    return(tibble())
  }
   
  if(tibble::has_rownames(tab)){
      tab <- tibble::rownames_to_column(tab)
    }
  
  tab <- try(tibble::as_tibble(tab))
  
  if(inherits(tab, "try-error")){
    message <- cat("as_tibble: ", tab[1])
    system(paste("echo '", path, "\n", message,"' >> log.txt"))
    return(tibble())
  }
  
  return(tab)
}

