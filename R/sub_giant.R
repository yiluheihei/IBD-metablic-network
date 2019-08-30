#' sub the giant componets of a network
#' 
#' @param g igrpah object
sub_giant <- function(g) {
  if (!inherits(g, "igraph")) {
    stop("`g` must be a igraph object")
  }
  giant_comp <- components(g)
  max_indx <- which.max(giant_comp$csize)
  giant_nodes <- igraph::groups(giant_comp)[[max_indx]]
  gg <- induced_subgraph(g, giant_nodes)
  
  g
}
