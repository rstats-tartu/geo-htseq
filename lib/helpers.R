

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
