#' generate topological properties tibble, including Betweenness centrality, Clustering coefficient, indegree, outdegree, degree
#' 
#' @param g igraph obejct, metablic network

generate_topos_tibble <- function(g ) {
  or <- vertex_attr(g, "OR")
  or_p <- vertex_attr(g, "OR_p")
  enrich_indx <- which(or > 2 & or_p < 0.05) # 15
  
  deplete_indx <-which(or < 0.5 & or_p < 0.05) # 31
  other_indx <- which((or >= 0.5 &  or <= 2) | or_p > 0.05)
  
  bc <- betweenness(g, normalized = TRUE)
  bc_tibble <- generate_s_topo_tibble(bc, g, "Betweenness centrality")
  
  cc <- transitivity(g, type ="weighted", vids = V(g), isolates=c("zero"))
  cc_tibble <- generate_s_topo_tibble(cc, g, "Clustering coefficient")
  
  indegree <- degree(g, mode = "in", normalized = TRUE)
  indegree_tibble <- generate_s_topo_tibble(indegree, g, "Indegree")
  
  outdegree <- degree(g, mode = "out", normalized = TRUE)
  outdegree_tibble <- generate_s_topo_tibble(outdegree, g, "Outdegree")
  
  degree <- degree(g, mode = "all", normalized = TRUE)
  degree_tibble <- generate_s_topo_tibble(degree, g, "Degree")
  
  
  
  return(
    bind_rows(bc_tibble, cc_tibble, indegree_tibble, outdegree_tibble, degree_tibble)
  )
}

#' generate a single topological property tibble, e.g. betweenness centrality
#' 
#' @param s_topo numerial vector value of topogical feature of all nodes
#' @param g igraph object
#' @param method character, specify the name of the topological property 
generate_s_topo_tibble <- function(s_topo, g, method) {
  or <- vertex_attr(g, "OR")
  or_p <- vertex_attr(g, "OR_p")
  enrich_indx <- which(or > 2 & or_p < 0.05) 
  
  deplete_indx <-which(or < 0.5 & or_p < 0.05) 
  other_indx <- which((or > 0.5 &  or < 2) | or_p >= 0.05)
  
  enrich_topo <- s_topo[enrich_indx]
  deplete_topo <- s_topo[deplete_indx]
  other_topo <- s_topo[other_indx]
  s_topo_tibble <- bind_rows(
    tibble(state = "Enrich", value = enrich_topo),
    tibble(state = "Deplete", value = deplete_topo),
    tibble(state = "Other", value = other_topo)
  )
  s_topo_tibble$method = method
  
  s_topo_tibble
}

#' index of the property of the KO, enrich, deplete or other
#' 
#' @param g igraph object
#' @return three length list
generate_associated_index <- function(g) {
  or <- vertex_attr(g, "OR")
  or_p <- vertex_attr(g, "OR_p")
  enrich_indx <- which(or > 2 & or_p < 0.05) 
  deplete_indx <-which(or < 0.5 & or_p < 0.05) 
  other_indx <- which((or > 0.5 &  or < 2) | or_p >= 0.05)
  
  res <- list(
    enrich_indx = enrich_indx,
    deplete_indx = deplete_indx,
    other_indx = other_indx
  )
  return(res)
  
}
