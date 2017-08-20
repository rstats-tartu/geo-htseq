
library(dplyr)
library(purrr)
library(magrittr)
library(stringr)
library(tidyr)

# Add srp to histograms ---------------------------------------------------

hist_srp <- function(ggp, SRP, pi0) {
  y <- max(layer_data(ggp)$ymax)
  size <- theme_get()$text[["size"]]/5
  family <- theme_get()$text[["family"]]
  ggp + annotate("text", 
                 label = paste("SRP ==", round(SRP, 2)),
                 x = .5,
                 y = y*.95,
                 size = size,
                 family = family,
                 parse = T) +
    annotate("text", 
             label = paste("pi*0 ==", round(pi0, 2)),
             x = .5,
             y = y*.85,
             size = size,
             family = family,
             parse = T)
}

# Fail-safe wrapper to calculate srp --------------------------------------

my_srp <- function(x, id=NULL, ...) {
  message(id)
  srp <- try(srp(x, ...), silent=T)
  if(inherits(srp, "try-error")){
    message(srp[1])
    return(srp[1])
  }
  srp
}


# Plot p value histogram --------------------------------------------------

phist <- function(df, id = NULL){
  message(id)
  lab <- str_split_fixed(id, "_", 2) %>% 
    str_split(pattern = "(?<=.xlsx)") %>% 
    vapply(paste, collapse = "\n", character(1)) %>% 
    paste(collapse="\n")
  p <- ggplot(df, aes(value)) +
    geom_histogram(breaks = seq(0, 1, 0.01)) +
    xlab("p values")
  yax <- max(layer_data(p)$ymax)*.9
  # p + annotate("text", x = 0.5, y = yax, label = lab, size = 1.5) 
  p + ggtitle(lab)
}

# Save grid plots ---------------------------------------------------------

save_gridplots <- function(gtab, filepath="Rplot", grdevice="pdf", ...){
  
  plotfun <- function(gtab, grdevice, path, ...){
    
    do.call(grdevice, list(file=path, ...))
    grid.draw(gtab)
    dev.off()
  }
  
  if(!is.gtable(gtab)){
    for(i in seq_along(gtab)){
      path <- file.path(paste(c(paste0(filepath, i), grdevice), collapse = "."))
      plotfun(gtab = gtab[[i]], grdevice = grdevice, path = path, ...)
    }
  } else {
    path <- file.path(paste(c(filepath, grdevice), collapse = "."))
    plotfun(gtab = gtab, grdevice = grdevice, path = path, ...)
  }
}

# Split list into smaller pices -------------------------------------------

chunkize <-  function(d, chunksize) split(d, ceiling(seq_along(d)/chunksize))

# get_pvalues function ---------------------------------------------------

get_pvalues <- function(restab, Accession = NULL){
  message(Accession)
  
  if(inherits(restab, "try-error")){
    return(NA)
  }
  
  rescols <- colnames(restab) %>% .[!is.na(.)] %>% .[!duplicated(.)]
  
  if(is.matrix(restab)){
    restab %<>% as.data.frame
  }
  
  # Fix colnames
  if(str_detect(tail(rescols, 1), "V[0-9]")){
    colnames(restab) <- c(tail(rescols, 1), setdiff(rescols, tail(rescols, 1)))
  }
  
  restab %<>% "["(rescols)
  restab <- restab[str_detect(colnames(restab), "base[Mm]ean|p(-)?val") & !str_detect(colnames(restab), "[Aa]dj|FDR|Corrected")]
  # restab %<>% select(matches("basemean|pval(ue)?(\\.[0-9])?$|^p-value"))
  
  if(ncol(restab) == 0 || !any(str_detect(colnames(restab), "p(-)?val"))){
    message('No column with p-values!')
    return(NULL)
  }
  
  if(any(str_detect(colnames(restab), "base[Mm]ean"))) {
    restab %<>% 
      mutate(basemean = rowMeans(select(., matches("basemean")))) %>% 
      filter(basemean > 1)
  }
  
  restab %<>% select(matches("p(-)?val")) %>% gather(na.rm = TRUE)
  
  if(class(restab$value) == "character"){
    restab %<>% mutate(value = readr::parse_number(value)) %>% filter(!is.na(value))
  }
  restab
}

# get_pvalues_basemean <- function(restab){
#   
#   if(inherits(restab, "try-error")){
#     return(NA)
#   }
#   
#   rescols <- colnames(restab) %>% .[!is.na(.)] %>% .[!duplicated(.)]
#   
#   if(is.matrix(restab)){
#     restab %<>% as.data.frame
#   }
#   
#   # Fix colnames
#   if(str_detect(tail(rescols, 1), "V[0-9]")){
#     colnames(restab) <- c(tail(rescols, 1), setdiff(rescols, tail(rescols, 1)))
#   }
# 
#   restab %<>% "["(rescols)
#   restab <- restab[str_detect(colnames(restab), "base[Mm]ean|^p(-)?val") & !str_detect(colnames(restab), "[Aa]dj|FDR|Corrected")]
#   
#   if(ncol(restab) == 0 || !any(str_detect(colnames(restab), "p(-)?val"))){
#     # message('No column with p-values!')
#     return(NULL)
#   }
#   
#   if(any(str_detect(colnames(restab), "base[Mm]ean"))) {
#     restab %<>% 
#       mutate(bmean = rowMeans(select(., matches("basemean"))))
#   }
#   
#   restab %<>% select(matches("bmean|p(-)?val"))
#   colnames(restab)[str_detect(colnames(restab),"p(-)?val")] <- "pvalue"
#   
#   if(any(str_detect(colnames(restab), "bmean"))){
#   colnames(restab)[str_detect(colnames(restab),"bmean")] <- "basemean"
#   } else {
#     restab$basemean <- NA
#   }
#   
#   return(restab)
# }

# # Destring function -------------------------------------------------------
# 
# destring <- function(x, keep = "0-9.-") {
#   return(as.numeric(gsub(paste("[^", keep, "]+", sep = ""), "", x)))
# }


# test integer function ---------------------------------------------------
testInteger <- function(x, id = NULL){
  message(id)
  if(!is.matrix(x)){
    message("Not matrix!")
    return(FALSE)
  }
  x <- x[sample(nrow(x), 100),]
  x <- sapply(x, function(z) as.integer(z) == z)
  all(sapply(x, all))
}

# Check for matrix rank ---------------------------------------------------

my_checkFullRank <- function(modelMatrix) {
  if (qr(modelMatrix)$rank < ncol(modelMatrix)) {
    if (any(apply(modelMatrix, 2, function(col) all(col == 0)))) {
      message("the model matrix is not full rank, so the model cannot be fit as specified.
           Levels or combinations of levels without any samples have resulted in
           column(s) of zeros in the model matrix.\n")
      return(TRUE)
    } else {
      message("the model matrix is not full rank, so the model cannot be fit as specified.
           One or more variables or interaction terms in the design formula are linear
           combinations of the others and must be removed.\n")
      return(TRUE)
    }
  }
  return(FALSE)
}

my_modelmatrix <- function(tmp){
  mm <- tmp %>% vapply(. %>% as.factor %>% as.numeric, numeric(nrow(tmp))) %>% as.data.frame
  design <- as.formula(paste("~", paste(colnames(mm), collapse = "+")))
  stats::model.matrix(design, mm)
}

# Infer exp design --------------------------------------------------------

remove_var <- function(tmp) {
  repeat {
    # remove col with max number of classes, when tie choose the last one (assuming that less important vars come later)
    tmp <- tmp[-max(which(vapply(tmp, n_distinct, integer(1)) == max(vapply(tmp, n_distinct, integer(1)))))]
    no_success <- group_by_(tmp, .dots=colnames(tmp)) %>% summarise(N=n()) %>% .$N %>% (function(x) any(x==1))
    # exit if the condition is met
    if (!no_success || ncol(tmp) <= 2) break
  }
  return(tmp)
}

is_chrnum <- function(x){
  all(!is.na(suppressWarnings(as.numeric(x))))
}

# metadata <- supptabs %>% filter(Accession=="GSE73857") %>% .$pdata %>% .[[1]]

get_design <- function(metadata, no.design = F){
  
  title <- metadata$title
  
  # Remove columns with NA-s
  nacols <- vapply(metadata, function(x) any(is.na(x)), logical(1))
  tmp <- metadata[, !nacols, drop = FALSE]
  
  # Remove columns with dates and extract protocols
  dates <- c("status", "submission_date", "last_update_date")
  protocols <- colnames(tmp)[str_detect(colnames(tmp), "extract_protocol")]
  tmp <- tmp[, setdiff(colnames(tmp), c(protocols, dates)), drop = FALSE]
  
  # Remove ftp & https links and publication data
  links <- vapply(tmp, function(z) all(str_detect(z, pattern = "ftp|http|Public")), logical(1))
  tmp <- tmp[, !links, drop = FALSE]
  
  # Replace plus and minus signs
  tmp <- tmp %>%
    vapply(str_replace, pattern = "\\+", replacement = "plus", character(nrow(tmp))) %>% as.data.frame() %>% 
    vapply(str_replace, pattern = "\\-", replacement = "minus", character(nrow(tmp))) %>% as.data.frame()
  
  # Some NAs are like "n/a" etc. or missing values like "", let's remove those too
  more_na <- vapply(tmp, function(x) any(str_detect(x, "\\bna|\\bn/a|NA|^$")), logical(1))
  tmp <- tmp[, !more_na, drop = FALSE]
  
  # Uninformative variables
  vars <- vapply(tmp, n_distinct, integer(1)) > 1 & vapply(tmp, n_distinct, integer(1)) != nrow(tmp)
  
  # Check if any of the variables remained and if not try to salvage some by chopping off 
  # two characters from the end of strings
  if(all(!vars)){
    tmp %<>% mutate_all(str_replace, pattern = ".{2}$", replacement = "")
    vars <- vapply(tmp, n_distinct, integer(1)) > 1 & vapply(tmp, n_distinct, integer(1)) != nrow(tmp)
    if(all(!vars)){
      (stop("Cannot identify independent variables in phenoData!"))
    }
  }
  
  # Remove vars with no variation and variables unique for each row
  tmp <- tmp[, vars, drop = FALSE]
  
  # Use title variable as rownames to match dependent variables with count data
  rownames(tmp) <- title
  
  # If no design 
  if(no.design){
    design <- as.formula("~ 1")
    return(list(colData = tmp, design = design))
  }
  
  # Test if we have only one variable left and can leave the function this time
  if(ncol(tmp)==1){
    message("Cool! One independent variable!")
    design <- as.formula(paste("~", colnames(tmp)))
    return(list(colData = tmp, design = design))
  }
  
  # if(ncol(tmp)>2){
  #   tmp <- remove_var(tmp)
  # }
  
  # Check matrix rank
  modelMatrix <- my_modelmatrix(tmp)
  while(qr(modelMatrix)$rank < ncol(modelMatrix)){
    trimm <- modelMatrix %>% "["(,-1) %>% cor %>% subselect::trim.matrix()
    vars <- setdiff(colnames(tmp), trimm$names.discarded)
    tmp <- tmp[, vars, drop = FALSE]
    modelMatrix <- my_modelmatrix(tmp)
  }
  
  # Drop highly correlated variables
  if(ncol(tmp)>1){
    cm <- modelMatrix %>% "["(,-1) %>% cor 
    cm[!lower.tri(cm)] <- 0 
    tmp <- tmp[,!apply(cm, 2, function(x) any(x > 0.8)), drop = FALSE]
  }
  
  design <- as.formula(paste("~", paste(colnames(tmp), collapse = "+")))
  list(colData = tmp, design = design)
}


# Arrange input data for model fitting ------------------------------------

model_input <- function(counts, metadata, id = NULL, ...){
  cat(id, "\n")
  
  colnames(counts) <- colnames(counts) %>% rm_punct_tolower()
  metadata$title <- metadata$title %>% rm_punct_tolower()
  countdata <- try(counts[, metadata$title, drop = FALSE], silent = T)
  
  if(inherits(countdata,"try-error")){
    message("No match between sample names and counts colnames!")
    return(data_frame())
  }
  
  # Convert strings in tables to NA and drop rows with NA-s
  countdata %<>% vapply(readr::parse_number, numeric(nrow(countdata))) %>% as_data_frame %>% filter(complete.cases(.))
  
  depvars <- try(get_design(metadata, ...), silent = T)
  if(inherits(depvars, "try-error")){
    message(depvars[1])
    return(data_frame())
  }
  
  c(list(countData = countdata), depvars)
}


# DESeq2 wrappers ---------------------------------------------------------

# Failsafe function to convert to deseq dataset  --------------------------
my_DESeqDataSetFromMatrix <- function(dataset, id=NULL, ...){
  message(id)
  argslist <- c(dataset, ...)
  dds <- try(do.call(DESeq2::DESeqDataSetFromMatrix, argslist), silent = T)
  if(inherits(dds, "try-error")){
    message(dds)
    return(dataset)
  }
  dds
}

my_DESeqDataSetFromMatrix2 <- function(countData, colData, design, id=NULL, ...) {
  message(id)
  dds <- try(DESeqDataSetFromMatrix(countData, colData, design, ...), silent = T)
  if(inherits(dds, "try-error")){
    warning(dds, immediate. = T)
    if(grepl("argument must be coercible", dds[1])){
      message("Fixing colnames.")
      rownames(countData) <- make.names(rownames(countData))
      dds <- DESeqDataSetFromMatrix(countData, colData, design, ...)
    }
  }
  dds
}

# Failsafe function to fit deseq model ------------------------------------
my_DESeq <- function(object, id=NULL, test=c("Wald","LRT"), reduced=~1, ...){
  message(id)
  
  args.list <- c(as.list(environment()), list(...))
  args.list["id"] <- NULL
  test <- match.arg(test)

  if(test=="Wald"){
    args.list["reduced"] <- NULL
  }
  
  try(do.call(DESeq2::DESeq, args.list))
}

# dds <- makeExampleDESeqDataSet()
# out <- my_DESeq(dds, id="fuck")
# results(out)

# Model input voom-limma --------------------------------------------------
#' @title Extract counts from supplementary table 
#' @description Match geo series matrix titles with supplementary table column names and keep this subset.
#' @param counts Data frame imported from geo series supplementary file.
#' @param samplenames Character vector. Geo series sample names obtained from series matrix file column 'title'.
#' @param id Character string. An optional supplementary file name.
#' @import dplyr
#' @import stringr
#' 
get_counts <- function(counts, samplenames, id = NULL){
  
  ## Output table id 
  message(id)
  
  ## Match metadata titles with counts column names
  colnames(counts) <- colnames(counts) %>% rm_punct_tolower()
  
  ## Prepend samplenames with 'x' when string starts with digit
  samplenames <-  rm_punct_tolower(samplenames) %>% 
    str_replace("(^[:digit:])", "x\\1")
  countdata <- try(counts[, samplenames, drop = FALSE], silent = T)
  
  ## When no match between series matrix titles and table colnames
  if(inherits(countdata,"try-error")){
    message("No match between sample names and counts colnames!")
    return(NULL)
  }
  
  ## Convert strings in tables to NA and drop rows with NA-s
  colclasses <- vapply(countdata, class, character(1)) %>% grepl("character", .)
  
  if(any(colclasses)){
    countdata <- vapply(countdata, readr::parse_number, numeric(nrow(countdata)))
  }
  
  ## If not matrix
  if(!is.matrix(countdata)){
    ## First column contains usually feature names
    rn <- counts[[1]]
    countdata <- as.matrix(countdata)
    row.names(countdata) <- rn
  }

  countdata %>% .[complete.cases(.),]
}

get_model <- function(metadata, moderated.vars=NULL, id=NULL){
  message(id)
  
  title <- metadata$title %>% rm_punct_tolower
  
  # Replace plus and minus signs
  tmp <- metadata %>%
    vapply(str_replace, pattern = "\\+", replacement = "plus", character(nrow(metadata))) %>% as.data.frame() %>% 
    vapply(str_replace, pattern = "\\-", replacement = "minus", character(nrow(metadata))) %>% as.data.frame()
  
  if(!is.null(moderated.vars)){
    message("Using provided variables.")
    tmp <- tmp[, moderated.vars, drop = FALSE]
    rownames(tmp) <- title
    return(tmp)
  }
  
  # Remove columns with NA-s
  nacols <- vapply(tmp, function(x) any(is.na(x)), logical(1))
  tmp <- tmp[, !nacols, drop = FALSE]
  
  # Some NAs are like "n/a" etc. or missing values like "", let's remove those too
  nacols <- vapply(tmp, function(x) any(str_detect(x, "\\bna|\\bn/a|NA|^$")), logical(1))
  tmp <- tmp[, !nacols, drop = FALSE]
  
  # Remove columns with dates 
  dates <- c("status", "submission_date", "last_update_date")
  tmp <- tmp[, setdiff(colnames(tmp), dates), drop = FALSE]
  
  # Remove ftp & https links and publication data
  links <- vapply(tmp, function(z) all(str_detect(z, pattern = "ftp|http|Public")), logical(1))
  tmp <- tmp[, !links, drop = FALSE]
  
  # Uninformative variables
  vars <- vapply(tmp, n_distinct, integer(1)) > 1 & vapply(tmp, n_distinct, integer(1)) != nrow(tmp)
  
  # # Remove variables denoting replicates but not unique...
  # tmp <- tmp[!vapply(tmp, function(x) all(grepl("\\<rep", x)), logical(1))]
  
  # Check if any of the variables remained and if not try to salvage some by chopping off 
  # two characters from the end of strings
  if(all(!vars)){
    tmp %<>% mutate_all(str_replace, pattern = ".{2}$", replacement = "")
    vars <- vapply(tmp, n_distinct, integer(1)) > 1 & vapply(tmp, n_distinct, integer(1)) != nrow(tmp)
    if(all(!vars)){
      warning("Cannot identify independent variables in phenoData!", immediate. = T)
      return(NULL)
    }
  }
  
  # Remove vars with no variation and variables unique for each row
  tmp <- tmp[, vars, drop = FALSE]
  
  # Use title variable as rownames to match dependent variables with count data
  rownames(tmp) <- title
  
  # Test if we have only one variable left and can leave the function this time
  if(ncol(tmp)==1){
    # design <- as.formula(paste("~", colnames(tmp)))
    message("One dependent variable!")
    return(tmp)
  }

  # Check matrix rank
  modelMatrix <- my_modelmatrix(tmp)
  while(qr(modelMatrix)$rank < ncol(modelMatrix)){
    trimm <- modelMatrix[,-1] %>% cor %>% subselect::trim.matrix()
    vars <- setdiff(colnames(tmp), trimm$names.discarded)
    tmp <- tmp[, vars, drop = FALSE]
    modelMatrix <- my_modelmatrix(tmp)
  }
  
  # Drop highly correlated variables
  if(ncol(tmp)>1){
    cm <- modelMatrix[,-1] %>% cor 
    cm[!lower.tri(cm)] <- 0 
    tmp <- tmp[, !apply(cm, 2, function(x) any(x > 0.8)), drop = FALSE]
  }
  
  # design <- as.formula(paste("~", paste(colnames(tmp), collapse = "+")))
  # list(object = design, data = tmp)
  tmp
}


# get_model_eset -----------------------------------------------------------

get_model_eset <- function(metadata, moderated.vars=NA, id=NULL){
  message(id)
  
  title <- metadata$title %>% rm_punct_tolower
  
  # Replace plus and minus signs
  tmp <- metadata %>%
    vapply(str_replace, pattern = "\\+", replacement = "plus", character(nrow(metadata))) %>% as.data.frame() %>% 
    vapply(str_replace, pattern = "\\-", replacement = "minus", character(nrow(metadata))) %>% as.data.frame()
  
  if(!is.na(moderated.vars)){
    message("Using provided variables.")
    tmp <- tmp[, moderated.vars, drop = FALSE]
    rownames(tmp) <- title
    return(tmp)
  }
  
  # Remove columns with NA-s
  nacols <- vapply(tmp, function(x) any(is.na(x)), logical(1))
  tmp <- tmp[, !nacols, drop = FALSE]
  
  # Some NAs are like "n/a" etc. or missing values like "", let's remove those too
  nacols <- vapply(tmp, function(x) any(str_detect(x, "\\bna|\\bn/a|NA|^$")), logical(1))
  tmp <- tmp[, !nacols, drop = FALSE]
  
  # Remove columns with dates 
  dates <- c("status", "submission_date", "last_update_date")
  tmp <- tmp[, setdiff(colnames(tmp), dates), drop = FALSE]
  
  # Remove ftp & https links and publication data and replicates columns
  links <- vapply(tmp, function(z) all(str_detect(z, pattern = "ftp|http|Public|rep(li)?|age|weight|cage")), logical(1))
  tmp <- tmp[, !links, drop = FALSE]
  
  # Uninformative variables
  vars <- vapply(tmp, n_distinct, integer(1)) > 1 & vapply(tmp, n_distinct, integer(1)) != nrow(tmp)
  
  # # Remove variables denoting replicates but not unique...
  # tmp <- tmp[!vapply(tmp, function(x) all(grepl("\\<rep", x)), logical(1))]
  
  # Check if any of the variables remained and if not try to salvage some by chopping off 
  # two characters from the end of strings
  if(all(!vars)){
    tmp %<>% mutate_all(str_replace, pattern = ".{2}$", replacement = "")
    vars <- vapply(tmp, n_distinct, integer(1)) > 1 & vapply(tmp, n_distinct, integer(1)) != nrow(tmp)
    if(all(!vars)){
      warning("Cannot identify independent variables in phenoData!", immediate. = T)
      return(NULL)
    }
  }
  
  # Remove vars with no variation and variables unique for each row
  tmp <- tmp[, vars, drop = FALSE]
  
  # Use title variable as rownames to match dependent variables with count data
  rownames(tmp) <- rownames(metadata)
  
  # Test if we have only one variable left and can leave the function this time
  if(ncol(tmp)==1){
    # design <- as.formula(paste("~", colnames(tmp)))
    message("One dependent variable!")
    return(tmp)
  }
  
  # Check matrix rank
  modelMatrix <- my_modelmatrix(tmp)
  while(qr(modelMatrix)$rank < ncol(modelMatrix)){
    trimm <- modelMatrix[,-1] %>% cor %>% subselect::trim.matrix()
    vars <- setdiff(colnames(tmp), trimm$names.discarded)
    tmp <- tmp[, vars, drop = FALSE]
    modelMatrix <- my_modelmatrix(tmp)
  }
  
  # Drop highly correlated variables
  if(ncol(tmp)>1){
    cm <- modelMatrix[,-1] %>% cor 
    cm[!lower.tri(cm)] <- 0 
    tmp <- tmp[, !apply(cm, 2, function(x) any(x > 0.67)), drop = FALSE]
  }
  
  # design <- as.formula(paste("~", paste(colnames(tmp), collapse = "+")))
  # list(object = design, data = tmp)
  tmp
}

