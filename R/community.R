#' Generate a microbial community 
#' 
#' Generate relative abundances of taxa from a beta distribution. Useful for illustrating and simulations.
#' 
#' @param CC Species richness of the community
#' @param alpha_parameter the alpha parameter for the beta distribution
#' @param beta_parameter the alpha parameter for the beta distribution
#' 
#' @export
make_community <- function(CC, alpha_parameter, beta_parameter) {
  relative_proportions <- rbeta(CC, alpha_parameter, beta_parameter) 
  data.frame("taxon" = paste("taxon", 1:CC, sep = ""),
             "abundances" = relative_proportions / sum(relative_proportions))
}

#' @export
community_richnness <- function(comm) {
  stopifnot(abs(comm$abundances %>% sum - 1) < 1e-8)
  nrow(comm)
}
