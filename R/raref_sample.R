#' network were constructed by considering an increasingly larger random subset of samples from each group
#'
#' @param ko_table ko abundance
#' @param state sample state
#' @param refnet reference metabolic network from kegg
#' 
#' @return two length list, health and disease network
raref_sample <- function(ko_table, state, refnet) {
  health_indx <- which(state == "healthy")
  health_count <- length(health_indx)
  health_table <- ko_table[, health_indx]
  n_sample <- length(state)
  disease_indx <- setdiff(1:n_sample, health_indx)
  disease_count <- n_sample - health_count
  disease_table <- ko_table[, disease_indx]
  
  raref_health_net <- map(1:health_count, ~ raref_net(.x, health_table, refnet))
  raref_disease_net <- map(1:disease_count, ~ raref_net(.x, disease_table, refnet))
  
  return(list(raref_health_net = raref_health_net, raref_disease_net = raref_disease_net))
}

#' metablic network of subset samples
raref_net <- function(indx, ko_table, refnet) {
  n <- ncol(ko_table)
  if (indx == 1) {
    kos <- row.names(ko_table)
    ko_abund <- ko_table[, sample(n, indx)]
    names(ko_abund) <- kos
    ko_abund <- ko_abund[ko_abund > 0]
    nodes <- intersect(V(refnet)$name, names(ko_abund))
  } else {
    selected_table <- ko_table[, sample(n, indx)]
    selected_table <- selected_table[rowSums(selected_table) > 0, ]
    nodes <- intersect(V(refnet)$name, row.names(selected_table))
  }
  
  res <- induced_subgraph(refnet, nodes) %>% sub_giant()
  return(res)
}
