#' natural connectivity of the network against attacks by sequentially removing hubs (nodes with high degree or betweenness centrality)
#' 
#' @param ig igraph object
#' @param method specify the feature to order the nodes
nc_attack <- function(ig, method = c("degree_betweenness", "degree", "betweenness", "random", "closeness")) {
  if (!method %in% c("degree_betweenness", "degree", "betweenness", "random", "closeness"))
    stop("method is one of degree_betweenness, degree, betweenness, random, closeness")
  oplan <- plan()
  on.exit(plan(oplan), add = TRUE)
  plan(multiprocess)
  if (method == "betweeness") {
    hubord <- order(rank(igraph::centr_betw(ig, directed = FALSE)$res), decreasing=TRUE)
  } else if (method == "closeness") {
    hubord <- order(rank(igraph::centr_clo(ig, mode = "all")$res), decreasing=TRUE)
  } else if (method == "degree_betweenness") {
    hubord <- order(rank(igraph::centr_betw(ig, directed = FALSE)$res), rank(igraph::degree(ig)), decreasing=TRUE)
  } else if (method == "degree") {
    hubord <- order(rank(igraph::degree(ig, normalized = TRUE)), decreasing = TRUE)
  } else {
    hubord <- sample(1:vcount(ig), vcount(ig))
  }
  net <- vector("list", round(vcount(ig)*0.8))
  res <- furrr::future_map(1:round(vcount(ig)*0.8), 
    function(i) {
      ind <- hubord[1:i]
      tmp <- delete_vertices(ig, V(ig)$name[ind])
      
      natcon(tmp)
    }
  )
  
  
  return(res)
}

#' natural connectivity
natcon <- function(ig) {
  N   <- vcount(ig)
  adj <- get.adjacency(ig)
  evals <- eigen(adj, only.values = TRUE)$value
  nc  <- log(mean(exp(evals)))
  nc / (N - log(N))
}

#' construct natural connectivity tibble for visualizaton using ggplot
make_nc_tibble <- function(nc, state = "HA", method = "degree_betweenness") {
  tibble(
    x = seq(0, 0.8,  len = length(nc)), 
    y = nc,
    state = rep(state, length(nc)),
    method = rep(method, length(nc))
  )
}