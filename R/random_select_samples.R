#' generate random sample index 
#' 
#' @param state sample state
#' @param n the number of samples to select
#' 
#' @return a two length list, containing the index of random selected sample
random_select_samples <- function(state, n = 50) {
  n_sample <- length(state)
  health_indx <- which(state == 'healthy')
  health_n <- sample(1:length(health_indx), n)
  selected_health <- health_indx[health_n]
  
  disease_indx <- which(state != 'healthy')
  disease_n <- sample(1:length(disease_indx), n)
  selected_disease <- disease_indx[disease_n]
  
  return(list(health_indx = selected_health, disease_indx = selected_disease))
}
