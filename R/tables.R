#' @import magrittr
#' @export
make_table <- function(Ca, Cb, aa, ab, ba, bb, na1, na2, nb1, nb2, replicates = 1) {
  
  population1 <- make_community(Ca, alpha_parameter = aa, beta_parameter = ab)
  population2 <- make_community(Cb, alpha_parameter = ba, beta_parameter = bb)
  
  # observations and truth
  truth <- rbind(data.frame("reads" = 0, 
                            "richness" = length(population1), 
                            "richness_max" = NA,
                            "shannon" = breakaway::shannon(population1), 
                            "shannon_max" = NA,
                            "ecosystem" =  "ecosystem1", 
                            "type" = "truth"),
                 data.frame("reads" = 0, 
                            "richness" = length(population2), 
                            "richness_max" = NA,
                            "shannon" = breakaway::shannon(population2), 
                            "shannon_max" = NA,
                            "ecosystem" =  "ecosystem2", 
                            "type" = "truth"))
  
  observations <- rbind(get_point(population1, na1, 1, "population1", replicates = replicates),
                        get_point(population1, na2, 2, "population1", replicates = replicates),
                        get_point(population2, nb1, 1, "population2", replicates = replicates),
                        get_point(population2, nb2, 2, "population2", replicates = replicates))
  
  rarefaction_level <- observations %>% 
    dplyr::filter(type == "observation") %>% 
    dplyr::select(reads) %>% 
    min
  
  rarefaction_table <- rbind(observations, truth)
  rarefaction_table %>% 
    dplyr::select(ecosystem) %>% 
    unlist %>% 
    gsub("[^[:digit:]]","",.) %>% 
    as.numeric -> rarefaction_table$ecosystem_number
  
  rarefaction_table
}

#' Estimate the probability of making the correct decision for different approaches
#' 
#' Accepts a data structure and returns the estimated error rates for different approaches
#' 
#' @param repeats Number of times to run the simulation
#' 
#' @import breakaway
#' @export
get_error_rates <- function(repeats, 
                            replicates,
                            Ca, Cb, 
                            aa, ab, ba, bb, 
                            na1, na2, nb1, nb2) {
  
  populationA <- make_community(Ca, alpha_parameter = aa, beta_parameter = ab)
  populationB <- make_community(Cb, alpha_parameter = ba, beta_parameter = bb)
  
  results <- data.frame()
  
  pvalues <- matrix(NA, nrow = 3, ncol = repeats)
  correct <- matrix(TRUE, nrow = 3, ncol = repeats)
  
  ## TODO: run in parallel
  for (i in 1:repeats) {
    
    # generate samples
    rr <- replicates
    read_counts <- round(c(seq(from = na1, to = na2, length.out = rr),
                           seq(from = nb1, to = nb2, length.out = rr)))
    samples <- list()
    for (k in 1:(2*rr)) {
      my_name <- get(paste("population", ifelse(k <= rr, "A", "B"), sep=""))
      my_sample <- sample_from_population(my_name, read_counts[k])
      
      samples[[k]] <- my_sample
    }
    
    # estimate richness
    cs <- lapply(samples, sample_richness) %>% unlist
    cs_rarefied <- lapply(samples, rarefied_richness, level = min(read_counts)) %>% unlist
    cc2 <- lapply(samples, estimate_richness_breakaway) %>% data.frame
    covariate <- rep(c("A", "B"), each = rr)
    cs_summary <- lm(cs ~ covariate) %>% summary
    cs_r_summary <- lm(cs_rarefied ~ covariate) %>% summary
    
    pvalues[1, i] <- cs_summary$coef[2,4]
    pvalues[2, i] <- cs_r_summary$coef[2,4]
    
    bta <- betta(cc2[1, ], cc2[2, ], 
                 X = cbind("Intercept" = 1, 
                           "CovariateB" = c(rep(0, rr), rep(1, rr) )))
    pvalues[3, i] <- bta$table[2,3]
    
    if (Ca != Cb) {
      correct[1, i] <- sign(cs_summary$coef[2,1]) == sign(Cb - Ca)
      correct[2, i] <- sign(cs_r_summary$coef[2,1]) == sign(Cb - Ca)
      
      ## TODO: is this right??
      # bta$table[1,3] is the reduction in richness due to covariate B
      # If Cb < Ca, should be negative
      # i.e. if sign(Cb - Ca) is -1, so should sign(bta$table[1,3])
      # Looks correct to me
      # TODO: actually should be [2,1]
      correct[3, i] <- sign(bta$table[1,3]) == sign(Cb - Ca)
    }
    
  }
  
  ## Measuring error rate
  ## If equal, H0 true; if p-value small, reject H0
  if(Ca == Cb) {
    ## Correctly detect no difference
    threshold <- pvalues > 0.05
  } else {
    ## Correctly detect difference
    threshold <- pvalues < 0.05
  }
  
  # To be correct overall, need to make correct decision in correct direction
  overall_correct <- threshold & correct
  correct <- apply(overall_correct, 1, mean)
  
  names(correct) <- c("Raw", "Rarefied", "Corrected")
  
  correct %<>% t %>% data.frame
  correct
}