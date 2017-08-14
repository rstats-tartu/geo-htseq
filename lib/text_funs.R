# Relevel function

relevel_fun <- function(x) {
  levels(x) <- rm_punct_tolower(levels(x))
  x
}

fix_levels <- function(x, id=NULL){
  message(id)
  x_rownames <- row.names(x)
  x <- mutate_all(x, relevel_fun)
  row.names(x) <- x_rownames
  x
}

# test integer ------------------------------------------------------------

test_int <- . %>% sapply(is.integer) %>% unlist %>% sum


# remove punctuation to lower ---------------------------------------------

rm_punct_tolower <- function(x){
  x <- stringr::str_replace_all(x, "\\+", "plus")
  x <- stringr::str_replace_all(x, "\\-", "minus")
  x <- stringr::str_replace_all(x, "[^[:alnum:]]", "")
  tolower(x)
}

# adjust colnames ---------------------------------------------------------

adj_colnames <- function(x) {
  newcolnames <- rm_punct_tolower(colnames(x))
  colnames(x) <- newcolnames
  return(x)
}

