#'  modularity of shuffled state-specific metabolic network
#' 
#' @param ko_table ko abundance table
#' @param state state of samples
#' @param refnet reference global metabolic network based on kegg pathway
#' 
#' @return two length list, healthy net modularity and disease net modularity
shuffle_modularity <- function(ko_table, state = state, refnet) {
  shuffled_net <- shuffle_label(ko_table, state = state, refnet)
  disease_net <- shuffled_net$disease_net
  health_net <- shuffled_net$health_net
  
  
  disease_modularity <- cluster_leading_eigen(as.undirected(disease_net), 
    options = list(maxiter = 1000000, ncv = 8)) %>% modularity()
  health_modularity <- cluster_leading_eigen(as.undirected(health_net), 
    options = list(maxiter = 1000000, ncv = 8)) %>% modularity()
  
  return(list(health_modularity = health_modularity, disease_modularity = disease_modularity))
}

shuffle_louvain_modularity <- function(ko_table, state = state, refnet) {
  shuffled_net <- shuffle_label(ko_table, state = state, refnet)
  disease_net <- shuffled_net$disease_net
  health_net <- shuffled_net$health_net
  
  
  disease_modularity <- cluster_louvain(as.undirected(disease_net)) %>% modularity()
  health_modularity <- cluster_louvain(as.undirected(health_net)) %>% modularity()
  
  return(list(health_modularity = health_modularity, disease_modularity = disease_modularity))
}

#' modularity of networks of shuffle edges 
shuffle_edge_modularity <- function(health_net, disease_net) {
  health_modularity <- randomize_net(health_net) %>%
    as.undirected() %>% 
    cluster_leading_eigen(options = list(maxiter = 1000000, ncv = 8)) %>% 
    modularity()
  disease_modularity <- randomize_net(disease_net) %>%
    as.undirected() %>% 
    cluster_leading_eigen(options = list(maxiter = 1000000, ncv = 8)) %>% 
    modularity()
  return(list(health_modularity = health_modularity, disease_modularity = disease_modularity))
}

#' modularity of networks of shuffle edges 
shuffle_edge_modularity2 <- function(net, health_kos, disease_kos) {
  random_net <- randomize_net(net) 
  
  health_kos <- intersect(V(random_net)$name, health_kos)
  disease_kos <- intersect(V(random_net)$name, disease_kos)
  health_net <- induced_subgraph(random_net, health_kos) %>% 
    sub_giant()
  disease_net <- induced_subgraph(random_net, disease_kos) %>% 
    sub_giant()
  
  health_modularity <- as.undirected(health_net) %>% 
    cluster_leading_eigen(options = list(maxiter = 1000000, ncv = 8)) %>% 
    modularity()
  disease_modularity <- as.undirected(disease_net) %>% 
    cluster_leading_eigen(options = list(maxiter = 1000000, ncv = 8)) %>% 
    modularity()
  
  return(list(health_modularity = health_modularity, disease_modularity = disease_modularity))
}

randomize_net <- function(net) {
  degree_out <- degree(net, mode = "out")
  degree_in <- degree(net, mode = "in")
  
  random_net <- sample_degseq(degree_out, degree_in, method = "simple")
  V(random_net)$name <- V(net)$name
  
  return(random_net)
}

#' shuffle sample labels
shuffle_label <- function(ko_table, state = state, refnet) {
  disease_count <- sum(state != "healthy")
  sample_count <- length(state)
  refnode <- V(refnet)$name
  
  disease_indx <- sample(1:sample_count, disease_count)
  health_indx <- setdiff(1:sample_count, disease_indx)
  disease_table <- ko_table[, disease_indx]
  # disease_table <- ko_table[, sample(sample_count, 50)]
  disease_table <- disease_table[rowSums(disease_table) > 0, ]
  health_table <- ko_table[, health_indx]
  #health_table <- ko_table[, sample(sample_count, 50)]
  health_table <- health_table[rowSums(health_table) > 0, ]
  disease_kos <- row.names(disease_table) %>% 
    intersect(refnode)
  health_kos <- row.names(health_table) %>% 
    intersect(refnode)
  disease_net <- induced_subgraph(refnet, disease_kos) %>% 
    sub_giant()
  health_net <- induced_subgraph(refnet, health_kos) %>% 
    sub_giant()
  
  return(list(health_net = health_net, disease_net = disease_net))
}

#' density of shuffled network
shuffle_density <- function(ko_table, state = state, refnet) {
  shuffled_net <- shuffle_label(ko_table, state = state, refnet)
  disease_net <- shuffled_net$disease_net
  health_net <- shuffled_net$health_net
  
  health_density <- edge_density(health_net)
  disease_density <- edge_density(disease_net)
  return(list(health_density = health_density, disease_density = disease_density))
}

#' shuffle edge density
shuffle_edge_density <- function(net, health_kos, disease_kos) {
  net <- rewire(net, keeping_degseq(niter = 20))
  
  health_kos <- intersect(V(net)$name, health_kos)
  disease_kos <- intersect(V(net)$name, disease_kos)
  health_net <- induced_subgraph(net, health_kos) %>% 
    sub_giant()
  disease_net <- induced_subgraph(net, disease_kos) %>% 
    sub_giant()
  
  health_density <- edge_density(health_net)
  disease_density <- edge_density(disease_net)
  return(list(health_density = health_density, disease_density = disease_density))
}

#' vertex count of shuffled network
shuffle_vcount <- function(ko_table, state, refnet) {
  shuffled_net <- shuffle_label(ko_table, state = state, refnet)
  disease_net <- shuffled_net$disease_net
  health_net <- shuffled_net$health_net
  
  health_vcount <- vcount(health_net)
  disease_vcount <- vcount(disease_net)
  return(list(health_vcount = health_vcount, disease_vcount = disease_vcount))
}
