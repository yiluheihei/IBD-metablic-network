#' construct integrated metabolic network
#' 
#' @param ssns list, a element represents a metabolic network of a given sample
#' @param sample.state character vector, state of all samples, e.g. IBD
#' @param diesase the disease state
diff_net <- function (ssns, sample.state, disease = "IBD") 
{
  abund <- lapply(ssns, function(x) {
    abund = get.vertex.attribute((x), name = "abundance")
    v.name = V(x)$name
    return(data.frame(v.name, abund))
  })
  abund <- Reduce(function(x, y) merge(x, y, by = "v.name", 
    all = TRUE), abund)
  abund[is.na(abund)] <- 0
  kos <- as.character(abund[, 1])
  ko.abund <- abund[, -1]
  if (length(sample.state) != ncol(ko.abund)) 
    stop("length of sample state must be equal to number of samples")
  state <- unique(sample.state)
  if (length(state) != 2) 
    stop("Differential analysis needed samples have two different state", 
      domain = NULL)
  index <- lapply(c(1, 2), function(x) which(match(sample.state,
    state) == x))
  names(index) <- state
  g <- induced_subgraph(RefDbcache$network, intersect(kos, 
    V(RefDbcache$network)$name))
  giant_comp <- components(g)
  giant_nodes <- groups(giant_comp)[[1]]
  g <- induced_subgraph(g, giant_nodes)
  #abund <- lapply(index, function(x) ko.abund[, x])
  kos_index <- match(V(g)$name, kos, nomatch = 0)
  ko.abund <- ko.abund[kos_index, ]
  row.names(ko.abund) <- V(g)$name
  or <- calculate_odds_ratio(ko.abund, sample.state)
  names(or) <- V(g)$name
  g <- set.vertex.attribute(g, name = "OR", value = or, index = V(g))
  
  # shuffled p value of odd ratio
  disease_count <- sum(sample.state == disease)
  shuffled_or <- replicate(1000, shuffle_or(ko.abund, disease_count = disease_count)) %>% 
    t() %>% 
    as_tibble()
  colnames(shuffled_or) <- V(g)$name
  or_h <- or[or >= 1]
  shuffled_or_h <- shuffled_or[, or >= 1]
  p_h <- map2_dbl(or_h, shuffled_or_h, ~ sum(.y > .x)/length(.y))
  or_l <- or[or < 1]
  shuffled_or_l <- shuffled_or[, or < 1]
  p_l <- map2_dbl(or_l, shuffled_or_l, ~ sum(.y < .x)/length(.y))
  shuffled_p <- c(p_h, p_l)
  shuffled_p <- shuffled_p[match(names(shuffled_or), names(shuffled_p))]
  g <- set.vertex.attribute(g, name = "OR_p", value = shuffled_p, index = V(g))
  
  return(g)
}
