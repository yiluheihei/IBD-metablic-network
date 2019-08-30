#' convert ko id to gene id
#' 
#' @param ko ko id
#' @param species species
#' 
#' @return a three length list, containing ko id,  gene id and species
ko2gene <- function(ko, species = "hsa") {
  url <- file.path("http://rest.kegg.jp/link", species, ko)
  content <- tryCatch(readLines(url), error=function(e) NULL) %>% 
    strsplit("\t") %>% 
    unlist()
  
  res <- list(ko = content[1], 
    gene_id = gsub("hsa:", "", content[2]),
    species = "hsa"
  )
  
  return(res)
}
