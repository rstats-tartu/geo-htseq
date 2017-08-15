
# get_GEO_Ids -------------------------------------------------------------

#' @title Get UID vector
#' @param query A GEO query string
#' @param database Entrez database, defaults to "gds" == GEO 
#' @param retmax maximum number of records to return, default is 500
#' @return A vector or GEO UIDs
#' @importFrom magrittr "%>%"
#' 
get_ids <- function(query, database = "gds", retmax = 500, ...){
  esearch <- "esearch.fcgi"
  base_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
  qres <- httr::GET(file.path(base_url, esearch), 
                    query = list(db = database, 
                                 term = query, 
                                 retmax = retmax, ...))
  rescont <- httr::content(qres)
  
  nids <- rescont %>% 
    xml2::xml_find_first(xpath = "//Count") %>% 
    xml2::xml_contents() %>% 
    xml2::xml_double()
  
  if(nids==0){
    stop("Query retrieved no results", call.=FALSE)
  } else {
    message(paste("Query retrieved", nids, "results"))
  }
  
  rescont %>% xml2::xml_find_all(xpath = "//Id") %>% xml2::xml_text()
}

# get_GEO_DocSums ---------------------------------------------------------

#' @param ids character vector of GEO UIDs
#' @param database Entrez database, defaults to "gds" == GEO 
#' @return A list of document summaries of class "xml_document" "xml_node"
#' 
get_docsums <- function(ids, database = "gds"){
  esummary <- "esummary.fcgi"
  base_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
  
  # Split UIDs into chunks of size max 500
  chunkize <-  function(d, chunksize) split(d, ceiling(seq_along(d)/chunksize))
  UID_chunks <- chunkize(ids, 500)
  
  get_qsums <- function(uid) {
    httr::GET(file.path(base_url, esummary), query = list(db = database, id = paste(uid, collapse = ",")))
  }
  
  qsums <- lapply(UID_chunks, get_qsums) 
  lapply(qsums, httr::content)
}


# extract_gds_docsums -----------------------------------------------------

#' @title Extract GEO DocSums into tibble
#' @param xmldoc A GEO query result contents, list of document summaries of class "xml_document" "xml_node"
#' @return A tibble of GEO document summaries
#' @importFrom magrittr "%>%"
#' 
extract_docsums <- function(xmldoc){
  
  d <- XML::xmlParse(xmldoc) 
  d <- XML::xmlRoot(d)
  
  items <- d[1]$DocSum %>% 
    XML::xmlSApply(XML::xmlGetAttr, name = "Name") %>% 
    unlist() %>% 
    c("Id",.) %>% 
    unname()
  
  d <- XML::xmlSApply(d, function(x) XML::getChildrenStrings(x)) 
  d <- dplyr::as_data_frame(t(d)) 
  colnames(d) <- items
  return(d)
}
