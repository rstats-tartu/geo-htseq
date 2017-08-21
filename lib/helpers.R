

# myversion function ------------------------------------------------------

myversion <- function(pkg){
  paste0("v", package.version(pkg))
}

# Function to format numbers ------
myround <- function(x, digits = 1){
  if(digits < 1)
    stop("This is intended for the case digits >= 1.")
  
  if(length(digits) > 1) {
    digits <- digits[1]
    warning("Using only digits[1]")
  }
  
  tmp <- sprintf(paste("%.", digits, "f", sep=""), x)
  
  # deal with "-0.00" case
  zero <- paste0("0.", paste(rep("0", digits), collapse=""))
  tmp[tmp == paste0("-", zero)] <- zero
  
  tmp
}

# Figure panel labels ----
add_labels <- function(groblist, case = c("upper", "lower"), ...){
  
  n_panels <- length(groblist)
  
  if(case=="upper"){
    label <- LETTERS[1:n_panels]
  } 
  
  if (case=="lower"){
    label <- letters[1:n_panels]
  }
  
  groblabels <- lapply(label, textGrob, x = unit(0,"npc"), hjust = 0, ...)
  mapply(function(x, y) arrangeGrob(x, y, ncol=1, heights = c(1, 12)), groblabels, groblist)
}

# Ref labels --------------------------------------------------------------
ref_labels <- function(label, case = c("upper", "lower")){
  if(case=="upper"){
    return(toupper(label))
  } 
  
  if (case=="lower"){
    return(tolower(label))
  } 
}


# geofile dataframe -------------------------------------------------------
#' Takes geo file name with accession name string and extracts Accession nr.
#' @param dir directory to look into for file names.
#' @param var a chracter string with column name for file names.
#' @return data_frame with two columns: 'Accession' and var.
#' @import dplyr
#' @example \notrun{
#' df <- geofile_df("data", "gsematrixfile")
#' }
#' 
geofile_df <- function(dir, var){
  
  ## List files and extract filenames
  files <- list.files(dir)
  acc <- stringr::str_extract(stringr::str_to_upper(files), "GSE[[:digit:]]*")

  ## Add accession nr and file names from dir into data_frame
  df <- data_frame(acc, files)
  setNames(df, c("Accession", var))
}

# get file extension ------------------------------------------------------

get_filext <- function(x) {
  x <- sub("\\.(gz|bz2)$", "", tolower(x))
  m <- regexpr("\\.[_A-z0-9]*$", x)
  regmatches(x, m)
}


