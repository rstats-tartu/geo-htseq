
#' @param path path to excel file
#' @param sheet name of the excel worksheet
#' 
rex <- function(path, sheet){
  
  tab <- try(readxl::read_excel(path, sheet), silent = TRUE)
  
  if (inherits(tab, "try-error")) {
    return(tibble())
  }
  return(tab)
}

# read excel sheets -------------------------------------------------------
#' Import all sheets from excel file.
#' @param path Path to excel file.
#' @return If file contains only one sheet returns data_frame. When more than one sheet is present, returns list of data_frames.
#' @import readxl
#' @import tibble
read_excelfs <- function(path) {
  
  sheets <- try(readxl::excel_sheets(path), silent = TRUE)
  
  if (inherits(sheets, "try-error")) {
    return(tibble()) 
  }
  
  tabs <- purrr::map(sheets, function(x) rex(path, sheet = x))
  names(tabs) <- sheets
  
  if (length(tabs) == 1) {
    tabs <- tabs[[1]]
  }
  return(tabs)
}


# Find duplicated columns -------------------------------------------------

find_duplicated_columns <- function(x){
  hashs <- vapply(x, function(x) digest::digest(x), character(1))
  duplicated(hashs)
}

# read_geotabs ------------------------------------------------------------
#' Import Entrez GEO supplemetary tables
#' @param path Path to supplementary file.
#' @return Returns data frame or list of data frames.
#' @import stringr
#' @import data.table
#' @import CePa
#' @import tibble
read_geotabs <- function(path){
  
  message(path)
  system(paste("echo '", path, "' >> log.txt"))
  
  if (stringr::str_detect(path, "xlsx?")) {
    tab <- read_excelfs(path)
    return(tab)
  }
  
  if (stringr::str_detect(path, "gct$")) {
    tab <- CePa::read.gct(path)
    return(tab)
  }
  
  if (stringr::str_detect(path, "\\.gz$")) {
    path <- paste0("zcat < ", path)
  }
  
  tab <- try(data.table::fread(path, data.table = FALSE))
  
  if (inherits(tab, "try-error")) {
    message <- tab[1]
    system(paste("echo '", path, "\n", message,"' >> log.txt"))
    return(tibble())
  }
  
  # Check for duplicated colnames
  coln <- colnames(tab)
  dups <- any(duplicated(coln))
  
  if (dups) {
    colnames(tab) <- make.unique(coln, sep = "_")
  }
  
  if (tibble::has_rownames(tab)) {
    tab <- tibble::rownames_to_column(tab)
  }
  
  tab <- try(tibble::as_tibble(tab))
  
  if (inherits(tab, "try-error")) {
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
  
  # Remove columns with all NA
  x <- x[colSums(!is.na(x)) > 0]
  
  # If matrix convert to data.frame
  if (is.matrix(x)) {
    x <- as.data.frame(x)
  }
  
  # Remove NA and duplicated column names
  colns <- colnames(x)
  colns <- colns[!is.na(colns)] 
  colns <- colns[!duplicated(colns)]
  
  if (length(colns) == 0) {
    stop("No colnames!")
  }
  
  # Fix colnames 
  if (stringr::str_detect(tail(colns, 1), "V[0-9]")) {
    colnames(x) <- c(tail(colns, 1), setdiff(colns, tail(colns, 1)))
  }
  
  x <- x[colns]
  
  # Now convert column names to lower case
  colnames(x) <- stringr::str_to_lower(colnames(x))
  
  pval_col <- stringr::str_detect(colns, "^p-?val(ue)?s?$") & !stringr::str_detect(colns, "adj|fdr|corrected")
  
  # Return NULL if P values not present
  if (!any(pval_col)) {
    return(NULL)
  }
  
  # Check if P values are between 0 and 1
  pvals_ok <- vapply(x[pval_col], function(x) {x <- range(x, na.rm = TRUE); all(x >= 0 & x <= 1)}, logical(1))
  
  if (!all(pvals_ok)) {
    stop("Detected P values not in 0 to 1 range or otherwise malformed!")
  }
  
  # Summarise multiple basemean columns
  if (any(stringr::str_detect(colns, "base[Mm]ean"))) {
    x <- dplyr::mutate(x, bmean = rowMeans(dplyr::select(x, dplyr::matches("base[Mm]ean"))))
  }
  
  # Select only basemean and pvalue olumns
  x <- dplyr::select(x, matches("bmean|^p(-)?val"))
  
  # Rename pvalue column
  colnames(x)[stringr::str_detect(colnames(x), "p(-)?val")] <- "pvalue"
  
  # Rename basemean column
  if (any(stringr::str_detect(colnames(x), "bmean"))) {
    colnames(x)[stringr::str_detect(colnames(x), "bmean")] <- "basemean"
  } else {
    x$basemean <- NA
  }
  return(x)
}

# munge_geo function ------------------------------------------------------

#' If geo supplementary table has only normalised read counts, return dims.
#' When table has pvalues return pvalues along with basemean if possible.
#' When table contains raw read counts, return eset.
#' @param suppfile GSE supplemental file name
#' @param eset nested esets
#' @param path path to GSE supplemental file, defaults to .
#' @import purrr
#' @import dplyr
#' @import Biobase
#' @import BiocGenerics
#' @import GEOquery
munge_geo <- function(eset, suppfile, dir = ".") {
  
  ## Assemble path to supplementary file
  path <- file.path(dir, suppfile)
  
  ## Import supplemental file, 
  supptab <- read_geotabs(path)
  
  ## Convert to list if single table
  if (is.data.frame(supptab) | is.matrix(supptab)) {
    supptab <- list(supptab) 
  }
  
  ## Check for pvalues
  pvalues <- purrr::map(supptab, ~try(get_pvalues_basemean(.x)))
  
  ## Exclude tables which failed to import
  pvalues <- pvalues[!purrr::map_lgl(pvalues, ~inherits(.x, "try-error"))]
  
  if (any(!purrr::map_lgl(pvalues, is.null))) {
    supptab <- supptab[purrr::map_lgl(pvalues, is.null)]
    pvalues <- pvalues[!purrr::map_lgl(pvalues, is.null)]
    if (length(supptab) == 0) {
      return(pvalues)
    }
  }
  
  ## Munge series matrix to eset
  esets <- eset %>% as.list %>% unlist
  
  ## Get sample names from eset
  titles <- purrr::map(esets, ~as.character(Biobase::pData(.x)[,"title"]))
  
  ## Get counts from table
  counts <- purrr::map(titles, ~get_counts(supptab[[1]], .x))
  
  ## Test for integers
  gplmatch <- purrr::map_lgl(counts, ~testInteger(.x))
  
  ## If not integers and no pvalues, return only dim
  if (sum(gplmatch) == 0) {
    
    dims <- unlist(purrr::map(supptab, dim))
    
    dims <- dplyr::data_frame(features = dims[1], columns = dims[2])
    
    if (all(purrr::map_lgl(pvalues, is.null))) {
      return(dims)
    }
    
    return(c(pvalues, dims))
  }
  
  ## Return raw read counts
  exprs <- BiocGenerics::counts[gplmatch][[1]]
  
  ## Add feature names
  rownames(exprs) <- make.names(rownames(exprs), unique = TRUE)
  
  ## Create eset
  title <- Biobase::pData(esets[gplmatch][[1]])[, 1, drop = FALSE]
  
  title$title  <- rm_punct_tolower(title$title)
  
  exprs <- Biobase::exprs[, title$title]
  
  colnames(exprs) <- rownames(title)
  
  phenoData <- new("AnnotatedDataFrame", data = Biobase::pData(esets[gplmatch][[1]]))
  
  neweset <- Biobase::ExpressionSet(assayData = exprs,
                                    phenoData = phenoData,
                                    experimentData = Biobase::experimentData(esets[gplmatch][[1]]),
                                    annotation = Biobase::annotation(esets[gplmatch][[1]]))
  
  return(c(pvalues, neweset))
}

#' Import pvalue column from Entrez GEO supplemental file.
#' @param suppfile Supplemental file name.
#' @param dir Directory of supplemental files. Defaults to current working directory.
#' @param n_head Number of rows to return from table head. defaults to n = 100.
#' @return A data_frame with columns: sheet, heads, pvalues, features and columns.
#' @import purrr
#' @import dplyr
munge_geo_pvalue <- function(suppfile, n_head = 100) {
  
  # Import supplemental file, 
  supptab <- read_geotabs(suppfile)
  
  # Check if import was successful or exit
  if (inherits(supptab, "try-error")) stop("Table import error!")
  
  # Convert to list if single table
  if (is.data.frame(supptab) || is.matrix(supptab)) {
    supptab <- list(supptab) 
  }
  
  # Check for pvalues
  pvalues <- purrr::map(supptab, ~try(get_pvalues_basemean(.x)))
  heads <- purrr::map(supptab, head, n = n_head)
  features <- purrr::map_int(supptab, nrow)
  columns <- purrr::map_int(supptab, ncol)
  
  # Get excel sheet names when available, otherwise return empty string
  supptab_names <- names(supptab)
  
  if (length(supptab_names) != 0) {
    
    sheets <- as.list(supptab_names)
    
  } else {
    
    sheets <- list("")
  }
  
  ## Collect data into data_frame
  dplyr::data_frame(sheets, heads, pvalues, features, columns)
}


# munge_geo2 --------------------------------------------------------------

munge_geo2 <- function(countfile, dir = ".") {
  
  # Assemble path to supplementary file
  path <- file.path(dir, countfile)
  
  # Import supplemental file
  supptab <- read_geotabs(path)
  
  if (inherits(supptab, "try-error")) {
    stop("Table import error!")
  }
  
  if (is.data.frame(supptab) | is.matrix(supptab)) {
    supptab <- list(supptab) 
  }
  
  # Check for pvalues
  pvalues <- map(supptab, ~try(get_pvalues_basemean(.x)))
  heads <- map(supptab, head)
  features <- map_int(supptab, nrow)
  columns <- map_int(supptab, ncol)
  supptab_names <- names(supptab)
  
  if (length(supptab_names) != 0) {
    sheets <- as.list(supptab_names)
  } else {
    
    sheets <- list("")
  }
  
  data_frame(sheets, heads, pvalues, features, columns)
  
}
