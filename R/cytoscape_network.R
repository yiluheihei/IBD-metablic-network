# library(RCy3)
# 
# # node color
# or <- vertex_attr(diff_ssn, "OR")
# or_p <- vertex_attr(diff_ssn, "OR_p")
# enrich_indx <- which(or > 2 & or_p < 0.05) 
# deplete_indx <-which(or < 0.5 & or_p < 0.05) 
# enrich_ko <- V(diff_ssn)$name[enrich_indx]
# deplete_ko <- V(diff_ssn)$name[deplete_indx]
# other_indx <- setdiff(1:length(V(diff_ssn)), c(enrich_indx, deplete_indx))
# type <- vector(length(V(diff_ssn)), mode = "character")
# type[enrich_indx] <- "enrich"
# type[deplete_indx] <- "deplete"
# type[other_indx] <- "other"
# 
# diff_ssn <- set_vertex_attr(diff_ssn, "node_type", value = type)
# 
# # OR可能出现极值，导入cytoscape的时候出现错误
# # https://github.com/cytoscape/cyREST/issues/99
# cytoscape_net <- diff_ssn
# cytoscape_net <- delete_vertex_attr(cytoscape_net, "OR")
# 
# createNetworkFromIgraph(cytoscape_net, "network")
# 
# setVisualStyle("Solid", "network")
# 
# setNodeSizeMapping(
#   "name", mapping.type = "d",
#   default.size = 20,
#   style.name = "Solid",
#   network = "network"
# )
# setNodeColorMapping(
#   "node_type",
#   c("enrich", "deplete", "other"),
#   colors = c("#E64B35", "#00A087", "#ABABAB"),
#   mapping.type = "discrete",
#   style.name = "Solid",
#   network = "network"
# )
# 
# 
# # remove node and edge label
# deleteStyleMapping(style.name = "Solid", 
#   visual.prop = "NODE_LABEL")
# deleteStyleMapping(style.name = "Solid", 
#   visual.prop = "EDGE_LABEL")
# 
# exportImage("output/nework.svg", type = "svg",
#   network = "network")
# 
# col2rgb(c("#E64B35", "#00A087", "#ABABAB"))

  