
#' @param countfile GSE supplemental file name
#' @param eset nested esets
#' @param path path to GSE supplemental file, defaults to .
munge_geo <- function(eset, countfile, path = ".") {
  # message(countfile)
  # Import supplemental file
  fullpath <- file.path(path, countfile)
  supptab <- readGEOtabs(fullpath)
  
  # Convert to list if single table
  if(is.data.frame(supptab)|is.matrix(supptab)){
    supptab <- list(supptab) 
  }
  
  # Check for pvalues
  pvalues <- map(supptab, ~try(get_pvalues_basemean(.x)))
  pvalues <- pvalues[!map_lgl(pvalues, ~inherits(.x, "try-error"))]
  
  if(any(!map_lgl(pvalues, is.null))){
    message("Found P-values!")
    supptab <- supptab[map_lgl(pvalues, is.null)]
    pvalues <- pvalues[!map_lgl(pvalues, is.null)]
    if(length(supptab)==0){
      return(pvalues)
    }
  }
  
  # Munge eset
  esets <- eset %>% as.list %>% unlist
  titles <- map(esets, ~as.character(pData(.x)[,"title"]))
  counts <- map(titles, ~get_counts(supptab[[1]], .x))
  
  # Test for integers
  gplmatch <- map_lgl(counts, ~testInteger(.x))
  
  if(sum(gplmatch)==0){
    dims <- unlist(map(supptab, dim))
    dims <- data_frame(features=dims[1], columns=dims[2])
    message("No raw reads!")
    
    if(all(map_lgl(pvalues, is.null))){
      return(dims)
    }
    
    return(c(pvalues, dims))
  }
  
  # message("Returning raw reads!")
  exprs <- counts[gplmatch][[1]]
  rownames(exprs) <- make.names(rownames(exprs), unique = T)
  
  title <- pData(esets[gplmatch][[1]])[,1, drop = FALSE]
  title$title  <- rm_punct_tolower(title$title)
  exprs <- exprs[, title$title]
  colnames(exprs) <- rownames(title)
  phenoData <- new("AnnotatedDataFrame", data = pData(esets[gplmatch][[1]]))
  neweset <- ExpressionSet(assayData = exprs,
                           phenoData = phenoData,
                           experimentData = experimentData(esets[gplmatch][[1]]),
                           annotation = annotation(esets[gplmatch][[1]]))
  
  c(pvalues, neweset)
}
