#' Sample species from a population 
#' 
#' @param population The names of the species
#' @param sample_depth The number of individuals to observe (the read depth)
#' 
#' @return The species labels of the sampled species
#' @export
sample_from_population <- function(population, 
                                   sample_depth) {
  
  data.frame("sample" = 1:sample_depth,
             "taxon" = sample(x = population$taxon,
                              size = sample_depth, 
                              prob = population$abundances, replace = T) %>%
               as.character)
  
}

#' @export
richness_from_sample <- function(smpl) {
  breakaway::sample_richness((smpl$taxon %>% 
                                table %>% 
                                data.frame)$Freq)$estimate
}

