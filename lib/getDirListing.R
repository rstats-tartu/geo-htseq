
getDirListing <- function(url) {
  message(url)
  # Takes a URL and returns a character vector of filenames
  a <- RCurl::getURL(url, dirlistonly = TRUE, ftp.use.epsv = FALSE, crlf = TRUE)
  ## Renaud Gaujoux reported problems behind firewall
  ## where the ftp index was converted to html content
  ## The IF statement here is his fix--harmless for the rest
  ## of us.
  if( grepl("^<HTML", a) ){ # process HTML content
    message("# Processing HTML result page (behind a proxy?) ... ", appendLF=FALSE)
    sa <- gsub('HREF', 'href', a, fixed = TRUE) # just not to depend on case change
    sa <- strsplit(sa, 'href', fixed = TRUE)[[1L]]
    pattern <- "^=\\s*[\"']/[^\"']+/([^/]+)[\"'].*"
    b <- as.matrix(gsub(pattern, "\\1", sa[grepl(pattern, sa)]))
    message('OK')
  } else { # standard processing of txt content
    tmpcon <- textConnection(a, "r")
    b <- read.table(tmpcon)
    close(tmpcon)
  }
  as.character(b[,ncol(b)])
}

getGEOSuppFileNames <- function(ftplink){
  url <- file.path(ftplink, "suppl/")
  fnames <- try(
    getDirListing(url)
    , silent = TRUE)
  if (inherits(fnames, "try-error")) {
    message(sprintf("No files found in %s.\nCheck URL manually if in doubt", url))
    return()
  }
  return(fnames)
}


# get dir listing tidy ----------------------------------------------------

#' Return current year when clock time is provided from entrez GEO ftp.
#' @description Entrez database shows clock time for files submitted on the current year.
#' @param x character string of clock time e.g. '6:21' or year.
#' @return character string with current year.
this_year <- function(x){
  if(str_detect(x, ":")) {
    x <- format(Sys.Date(), "%Y")
  }
  return(x)
}

#' Returns supplementary file names from Entrez GEO.
#' @description Replaces GEOquery packages getDirlisting function, seems to be faster.
#' @param r Object of class 'response'.
#' @return A data frame with columns date, size and suppfile. 
#' Date is the supplementary file date, size is filesize in bytes. 
#' Suppfile is the supplementary file name.
#' @examples \notrun{
#' 
#' ## Download supplementary file name
#' url <-  "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81555//suppl/"
#' r <- httr::GET(url)
#' suppfiles <- get_dirlist(r)
#' 
#' ## Download series matrix file name
#' url <-  "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81555//miniml/"
#' r <- httr::GET(url)
#' matrixfile <- get_dirlist(r)
#' 
#' }
#' @import dplyr
#' @import httr
#' @export
get_dirlist <- function(r){
  
  if(class(r) != "response") {
    warning("Not html response!")
    return()
  }
  
  cont <- content(r, as = "text", encoding = "UTF-8")
  tb <- readr::read_delim(cont, "\n", col_names = F)
  tb <- tidyr::separate(tb, "X1", paste0("C", 1:9), "[[:space:]]+")
  tb <- dplyr::mutate(tb, year = purrr::map_chr(C8, this_year),
                      month = which(stringr::str_detect(month.abb, C6)),
                      date = lubridate::dmy(paste(C7, month, year, sep = "-")))
  tb <- dplyr::select(tb, date, C5, C9)
  colnames(tb) <- c("date", "size", str_extract(r$url, "miniml|suppl"))
  return(tb)
}

# download supplementary files --------------------------------------------

geo_supp_dwnl <- function(ftplink, slug = c("suppl", "miniml"), filename, destdir) {
  
  message(filename)
  
  fp <- file.path(ftplink, slug, filename)
  dest <- file.path(destdir, filename)
  
  if(file.exists(dest)){
    message("File exists!")
    return()
  }
  
  out <- try(download.file(fp, dest))
  if(inherits(out, "try-error")){
    message("Download failed!")
  }
}


# matrix files ------------------------------------------------------------

getGEOSeriesFiles <- function(ftplink, subdir = "miniml/"){
  url <- file.path(ftplink, subdir)
  fname <- try(getDirListing(url), silent = TRUE)
  if (inherits(fname, "try-error")) {
    message(sprintf("No files found in %s.\nCheck URL manually if in doubt", url))
    return(NULL)
  }

  download.file(file.path(url, fname),
                destfile = file.path("data", subdir, fname),
                mode = 'wb')
}


