
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
    warning("No supplemental files found. Check URL manually if in doubt.")
    message(url)
    return()
  }
  return(fnames)
}


# download supplementary files --------------------------------------------

geo_supp_dwnl <- function(ftplink, suppfile, destdir) {
  message(suppfile)
  
  fp <- file.path(ftplink, "suppl", suppfile)
  dest <- file.path(destdir, suppfile)
  
  if(file.exists(dest)){
    message("File exists!")
    return()
  }
  
  out <- try(download.file(fp, dest))
  if(inherits(out, "try-error")){
    message("Ups! Something went wrong!")
  }
}


# matrix files ------------------------------------------------------------

getGEOSeriesFiles <- function(ftplink, subdir = "miniml/"){
  url <- file.path(ftplink, subdir)
  fname <- try(getDirListing(url), silent = TRUE)
  if (inherits(fname, "try-error")) {
    message("No files found.")
    message("Check URL manually if in doubt")
    message(url)
    return(NULL)
  }

  download.file(file.path(url, fname),
                destfile = file.path("data", subdir, fname),
                mode = 'wb')
}


