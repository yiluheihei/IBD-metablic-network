#' Regroup uniref90 table to ko table
#'
#' @param gene_family_abundance uniref90 abundance
#' @param ko_uniref90 table of ko to uniref90
#' 
#' @return ko abundance table, the same result of humann2_regroup_table 
regroup_ko_table <- function(gene_family_abundance, ko_uniref90) {
  kos <- purrr::map_chr(ko_uniref90, ~ .x[1])
  uniref90s <- purrr::map(ko_uniref90, ~.x[-1])
  
  dplan <- plan()
  on.exit(plan(dplan))
  plan(multiprocess)
  kos_table <- future_map(
    uniref90s, 
    sum_ko_uniref90, gene_family_abundance = gene_family_abundance,
    .progress = TRUE
  ) %>%
    do.call(rbind, .)
  
  row.names(kos_table) <- kos
  kos_indx <- which(rowSums(kos_table) != 0)
  kos_table <- kos_table[kos_indx, ]
  
  kos_table
}

sum_ko_uniref90 <- function(gene_family_abundance, uniref90) {
  features <- row.names(gene_family_abundance)
  indx <- features %in% uniref90
  
  ko_table <- gene_family_abundance[indx, ] 
  
  if (is.matrix(ko_table)) {
    colSums(ko_table)
  } else {
    ko_table
  }
}

