
#' Get text or attributes from xml.
#' @param node xml_node object.
#' @import xml2
get_nodedata <- function(node){
  
  ## Look for node text
  nodedata <- xml2::xml_text(node)
  nodename <- xml2::xml_name(node)
  names(nodedata) <- nodename
  
  ## If text is missing return attributes
  if(nchar(nodedata)==0){
    
    nodedata <- xml2::xml_attrs(node)
    
    if(nodename=="instrument") {
      
      names(nodedata) <- nodename
      
    } else {
      
      names(nodedata) <- paste0(nodename, '.', names(nodedata))
    }
  }
  
  return(nodedata)
}

#' Get text or attributes from xml.
#' @inheritParams get_nodedata
#' @import xml2
get_deepnode <- function(node){
  deepnode <- xml_child(node[[5]])
  deepnode_name <- xml_name(deepnode)
  deepnode_data <- xml_text(deepnode)
  # names(deepnode_data) <- paste0(deepnode_name, '.', names(deepnode_data))
  deepnode_data
}

#' Wrapper around get_nodedata
#' @inheritParams get_nodedata
#' @import xml2
get_nodedata_wrap <- function(node) {
  
  if(xml2::xml_length(node)>0) {
    nodename <- xml2::xml_name(node)
    node <- xml2::xml_contents(node)
    nodelist <- lapply(node, function(x) {
      if(xml2::xml_length(x)>0){
        return(get_deepnode(x))
      } else {
        get_nodedata(x)
      }
    })
    return(nodelist)
  } 
  
  get_nodedata(node)
}

#' Extracts data from SRA document summary ExpXml field into data frame.
#' @param xmlstring A character string with xml data. 
#' @import xml2
extract_expxml <- function(xmlstring){
  
  if(!is.character(xmlstring)&&!grepl("^<Summary><Title>", xmlstring)) stop("Not proper xml string!")
  
  ## Convert xml string into html
  xmlstr <- xml2::read_html(xmlstring)
  
  ## Extract html body with xml
  body <- xml2::xml_child(xmlstr, "body")
  body_contents <- xml2::xml_contents(body)
  body_cont_len <- xml2::xml_length(body_contents)
  
  ## Extract xml contents into list
  body_contents_list <- lapply(body_contents, get_nodedata_wrap)
  
  ## Convert nested list into named vector
  unlist(body_contents_list, use.names = TRUE)
}
