


# read excel sheets -------------------------------------------------------

read_excelfs <- function(path) {
  
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


# read_geotabs ------------------------------------------------------------

read_geotabs <- function(path){
  
  message(path)
  system(paste("echo '", path, "' >> log.txt"))
  
  if(stringr::str_detect(path, "xlsx?")){
    tab <- read_excelfs(path)
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

# get_pvalues_basemean ----------------------------------------------------
#' Extract p value and when available basemean columns from geo data frames.
#' @description Function tries to identify and subscript raw p value and base mean columns from data frame. When multiple basemean columns are detected, function calculates rowmeans.
#' @param x data frame imported from Entrez GEO supplementary file.
#' @return Data frame with pvalue columns and basemean column.
#' @examples \notrun {
#' library(data.table)
#' ## Dataset misses base mean column
#' x <- fread(sprintf("zcat < %s", "inst/extdata/GSE66793_cmp1-3-geo-gene.tsv.gz"))
#' pv <- get_pvalues_basemean(as.data.frame(x))
#' }
#' @import dplyr
#' @import stringr
get_pvalues_basemean <- function(x){
  
  ## Remove NA and duplicated columns
  colns <- colnames(x)
  colns <- colns[!is.na(colns)] 
  colns <- colns[!duplicated(colns)]
  
  ## If matrix convert to data.frame
  if(is.matrix(x)){
    x <- as.data.frame(x)
  }
  
  ## Fix colnames 
  if(str_detect(tail(colns, 1), "V[0-9]")){
    colnames(x) <- c(tail(colns, 1), setdiff(colns, tail(colns, 1)))
  }
  
  x <-  x[colns]
  pval_col <- str_detect(colns, "^p(-)?val") & !str_detect(colns, "[Aa]dj|FDR|Corrected")
  
  ## If P value column present, check if P values are between 0 and 1 
  if(any(pval_col)){
    
    pvals_ok <- vapply(x[pval_col], function(x) {all(na.omit(x)>=0 & na.omit(x)<=1)}, logical(1))
    
    if(!all(pvals_ok)){
      
      stop("Detected P values not in 0 to 1 range or otherwise malformed!")
      
    }
    
  } else {
    
    return(NULL) 
    
  }
  
  basemean_col <- str_detect(colns, "base[Mm]ean")
  # restab <- restab[mapply(any, pval_col, basemean_col)]
  
  ## Summarise multiple basemean columns
  if(any(basemean_col)) {
    x <- mutate(x, bmean = rowMeans(select(x, matches("basemean"))))
  }
  
  ## Select only basemean and pvalue olumns
  x <- select(x, matches("bmean|p(-)?val"))
  
  ## Rename pvalue column
  colnames(x)[str_detect(colnames(x),"p(-)?val")] <- "pvalue"
  
  ## Rename basemean column
  if(any(str_detect(colnames(x), "bmean"))){
    colnames(x)[str_detect(colnames(x),"bmean")] <- "basemean"
  } else {
    x$basemean <- NA
  }
  
  return(x)
}

# munge_geo function ------------------------------------------------------

geosupplement <- function(pvals = list(), dims = list(), eset = ExpressionSet()) {
  
  ## Assign outputs to list
  geo <- list(pvalues = pvals, dims = dims, eset = eset)
  
  ## Set the name for the class
  class(geo) <- append(class(geo), "geosupplement")
  return(geo)
}

#' If geo supplementary table has only normalised read counts, return dims.
#' When table has pvalues return pvalues along with basemean if possible.
#' When table contains raw read counts, return eset.
#' @param countfile GSE supplemental file name
#' @param eset nested esets
#' @param path path to GSE supplemental file, defaults to .
#' @import purrr
#' @import dplyr
#' @import GEOquery
munge_geo <- function(eset, countfile, dir = ".") {
  
  ## Assemble path to supplementary file
  path <- file.path(dir, countfile)
  
  ## Import supplemental file, 
  supptab <- read_geotabs(path)
  
  ## Convert to list if single table
  if(is.data.frame(supptab)|is.matrix(supptab)){
    supptab <- list(supptab) 
  }
  
  ## Check for pvalues
  pvalues <- map(supptab, ~try(get_pvalues_basemean(.x)))
  
  ## Exclude tables which failed to import
  pvalues <- pvalues[!map_lgl(pvalues, ~inherits(.x, "try-error"))]
  
  if(any(!map_lgl(pvalues, is.null))){
    
    supptab <- supptab[map_lgl(pvalues, is.null)]
    
    pvalues <- pvalues[!map_lgl(pvalues, is.null)]
    
    if(length(supptab)==0){
      
      return(pvalues)
    }
  }
  
  ## Munge series matrix to eset
  esets <- eset %>% as.list %>% unlist
  
  ## Get sample names from eset
  titles <- map(esets, ~as.character(pData(.x)[,"title"]))
  
  ## Get counts from table
  counts <- map(titles, ~get_counts(supptab[[1]], .x))
  
  ## Test for integers
  gplmatch <- map_lgl(counts, ~testInteger(.x))
  
  ## If not integers and no pvalues, return only dim
  if(sum(gplmatch)==0){

    dims <- unlist(map(supptab, dim))
    
    dims <- data_frame(features = dims[1], columns = dims[2])
    
    if(all(map_lgl(pvalues, is.null))){
      
      return(dims)
    
      }
    
    return(c(pvalues, dims))
  }
  
  ## Return raw read counts
  exprs <- counts[gplmatch][[1]]
  
  ## Add feature names
  rownames(exprs) <- make.names(rownames(exprs), unique = T)
  
  ## Create eset
  title <- pData(esets[gplmatch][[1]])[, 1, drop = FALSE]
  
  title$title  <- rm_punct_tolower(title$title)
  
  exprs <- exprs[, title$title]
  
  colnames(exprs) <- rownames(title)
  
  phenoData <- new("AnnotatedDataFrame", data = pData(esets[gplmatch][[1]]))
  
  neweset <- ExpressionSet(assayData = exprs,
                           phenoData = phenoData,
                           experimentData = experimentData(esets[gplmatch][[1]]),
                           annotation = annotation(esets[gplmatch][[1]]))
  
  return(c(pvalues, neweset))
}

munge_geo2 <- function(countfile, dir = ".") {
  
  ## Assemble path to supplementary file
  path <- file.path(dir, countfile)
  
  ## Import supplemental file, 
  supptab <- read_geotabs(path)
  
  if(inherits(supptab, "try-error")){
    stop("Table import error!")
  }
  
  ## Convert to list if single table
  if(is.data.frame(supptab)|is.matrix(supptab)){
    supptab <- list(supptab) 
  }
  
  ## Check for pvalues
  pvalues <- map(supptab, ~try(get_pvalues_basemean(.x)))
  heads <- map(supptab, head)
  features <- map_int(supptab, nrow)
  columns <- map_int(supptab, ncol)
  supptab_names <- names(supptab)
  
  if(length(supptab_names)!=0){
    sheets <- as.list(supptab_names)
  } else {
    sheets <- list("")
  }
  
  data_frame(sheets, heads, pvalues, features, columns)
}
