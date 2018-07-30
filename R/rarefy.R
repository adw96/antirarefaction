#' Rarefy a sample to a given level
#'
#' @param smpl See sample_richness for details.
#' 
#' @export
rarefied_sample <- function(smpl, level) {
  smpl[1:level, ]
}