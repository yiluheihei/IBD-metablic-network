#' calculate the odds ratio of each ko of metagenomic samples
#' 
#' @param kos_table a data frame, ko abundance of all samples, each variable represents a sample while each row represents a ko.
#' @param state a vector whose length must be equal to sample count, specifiy the state fo the sample
#' @param disease character, disease state
#' 
#' @return vector, whose length is equal to  sample count
calculate_odds_ratio <- function(kos_table, state, disease = "IBD") {
  disease_table <- kos_table[, state == disease]
  health_table <- kos_table[, state != disease]
  
  disease_ratio <- get_ratio(disease_table, feature_row = TRUE)
  health_ratio <- get_ratio(health_table, feature_row = TRUE)
  
  disease_ratio/health_ratio
}

#' shuffle the smaple lables and calculate the odds ratio
#' 
#' @param ko_table a data frame, ko abundance of all samples, each variable represents a sample while each row represents a ko.
#' @param disease_count integer, the number of samples in disease state
#' 
#' @return vector, whose length is equal to  sample count
shuffle_or <- function(ko_table, disease_count) {
  
  sample_count <- ncol(ko_table)
  disease_indx <- sample(1:sample_count, disease_count)
  health_indx <- setdiff(1:sample_count, disease_indx)
  
  disease_table <- ko_table[, disease_indx]
  health_table <- ko_table[, health_indx]
  
  disease_ratio <- get_ratio(disease_table, feature_row = TRUE)
  health_ratio <- get_ratio(health_table, feature_row = TRUE)
  
  or <- disease_ratio/health_ratio
  names(or) <- row.names(ko_table)
  
  
  or
}

#' ratio
get_ratio <- function(abundance, feature_row = FALSE) {
  if (feature_row) {
    ko_abd <- rowSums(abundance)
  } else {
    ko_abd <- colSums(abundance)
  }
  
  no_ko_abd <- sum(abundance) - ko_abd
  
  ko_abd/no_ko_abd
} 

