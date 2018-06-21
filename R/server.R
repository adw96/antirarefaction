#' Generate a microbial community 
#' 
#' @export
make_community <- function(CC, alpha_parameter, beta_parameter) {
  relative_proportions <- rbeta(CC, alpha_parameter, beta_parameter) 
  raw_abundances <- relative_proportions/sum(relative_proportions)
  raw_abundances
}

sample_richness <- function(community) {
  community %>% unique %>% length
}

rarefied_richness <- function(my_sample, level) {
  my_sample %>% sample(size = level, replace = FALSE) -> rarefied_sample
  rarefied_sample %>% unique %>% length
}

estimate_richness_breakaway <- function(my_sample) {
  my_sample %>% table %>% make_frequency_count_table -> my_count_table
  result <- try(breakaway::breakaway(my_count_table, output=F, plot = F, answers = T), 
                silent = T)
  if (class(result) == "list")  {
    c(result$est, result$seest)
  } else {
    print(my_count_table)
    c(NA, NA)
  }
}

draw_rarefaction <- function(population, pts = NA, replicates = 5, ecosystem = NA) {
  rarefaction_table <- data.frame()
  if (all(is.na(pts))) pts <- c(round(seq(from = 5, to = length(population$names)/100, length = 100)),
                                round(seq(from = length(population$names)/100, to = length(population$names)/10, length = 100)),
                                round(seq(from = length(population$names)/10, to = length(population$names)/4, length = 100)),
                                round(seq(from = length(population$names)/4, to = length(population$names), length = 10)))
  who <- replicate(replicates, sample(population$names, length(population$names), replace = F))
  for (i in pts) {
    subsample_diversity <- rep(NA, replicates)
    subsample_shannon <- rep(NA, replicates)
    for (j in 1:replicates) {
      subsample <- who[1:i, j]
      subsample_diversity[j] <- length(unique(subsample))
      community <- c(as.matrix(table(subsample)))
      proportions <- community/sum(community)
      subsample_shannon[j] <- breakaway::shannon(proportions)
    }
    rarefaction_table %<>% rbind(data.frame("reads" = i,
                                            "richness" = min(subsample_diversity), "richness_max" = max(subsample_diversity),
                                            "shannon" = min(subsample_shannon), "shannon_max" = max(subsample_shannon),
                                            "ecosystem" = ecosystem, "type" = "rcurve"))
  }
  rarefaction_table
}

get_point <- function(ecosystem, read, label) {
  read <- round(read)
  my_population <- get(paste("population", ecosystem, sep=""))
  ecosystem_name <- paste("ecosystem", ecosystem, sep="")
  subsample <- sample(x = c(1:length(my_population$proportions)),
                      size = read, prob = my_population$proportions, replace = T)
  
  point <- data.frame("reads" = read,
                      "richness" = length(unique(subsample)), "richness_max" = NA,
                      "shannon" = breakaway::shannon(c(as.matrix(table(subsample)))/sum(c(as.matrix(table(subsample))))), "shannon_max" = NA,
                      "ecosystem" = ecosystem_name, "type" = "observation")
  subsample1 <- list()
  subsample1$names <- subsample
  subsample1$proportions <- c(as.matrix(table(subsample)))/sum(c(as.matrix(table(subsample))))
  rarefied_sample <- draw_rarefaction(subsample1, replicates = 1, pts = seq(from = 1, to = read, length = 100),
                                      ecosystem = ecosystem_name)
  rarefied_sample$type <- paste("rarefied_sample", label, sep="")
  xx <- subsample %>% table %>% make_frequency_count_table
  
  estimate <- xx %>% breakaway(answers = T, plot = F, output=F) %$% est
  error <- xx %>% breakaway(answers = T, plot = F, output=F) %$% seest
  richness_estimate <- c(estimate - 2*error, estimate + 2*error)
  
  estimate <- data.frame("reads" = read,
                         "richness" = richness_estimate[1], "richness_max" = richness_estimate[2],
                         "shannon" = 0, "shannon_max" = 100000,
                         "ecosystem" = ecosystem_name, "type" = paste("estimate", label, sep=""))
  rbind(point, rarefied_sample, estimate)
}

make_table <- function(Ca, Cb, aa, ab, ba, bb, na1, na2, nb1, nb2) {
  
  population1 <- make_community(Ca, alpha_parameter = aa, beta_parameter = ab)
  population2 <- make_community(Cb, alpha_parameter = ba, beta_parameter = bb)
  
  # observations and truth
  truth <- rbind(data.frame("reads" = 0, 
                            "richness" = length(population1$proportions), 
                            "richness_max" = NA,
                            "shannon" = breakaway::shannon(population1$proportions), 
                            "shannon_max" = NA,
                            "ecosystem" =  "ecosystem1", 
                            "type" = "truth"),
                 data.frame("reads" = 0, 
                            "richness" = length(population2$proportions), 
                            "richness_max" = NA,
                            "shannon" = breakaway::shannon(population2$proportions), 
                            "shannon_max" = NA,
                            "ecosystem" =  "ecosystem2", 
                            "type" = "truth"))
  
  
  
  observations <- rbind(get_point(1, na1, 1),
                        get_point(1, na2, 2),
                        get_point(2, nb1, 1),
                        get_point(2, nb2, 2))
  
  rarefaction_level <- observations %>% filter(type == "observation") %>% select(reads) %>% min
  
  rarefaction_table <- rbind(observations, truth)
  rarefaction_table %>% 
    select(ecosystem) %>% 
    unlist %>% 
    substring(10, 11) %>% 
    as.numeric -> rarefaction_table$ecosystem_number
  
  rarefaction_table_unequal <- rarefaction_table
  rarefaction_table_unequal
}

make_plot <- function(rarefaction_table_unequal,
                      Ca, Cb, aa, ab, ba, bb, na1, na2, nb1, nb2) {
  ## TODO: define Ca, Cb, na1, na2...
  plot_base <- ggplot(data = rarefaction_table_unequal, aes(col = ecosystem)) +
    theme_bw() +
    geom_hline(data = rarefaction_table_unequal %>% filter(type == "truth"),
               aes(yintercept=richness, col = ecosystem, lty = ecosystem)) +
    scale_linetype_manual(guide = "none", values=c("twodash", "dotted")) +
    labs(x="No. reads", y="Taxonomic richness",col = "") +
    coord_cartesian(ylim=c(0, max(Ca, Cb)*1.1))
  
  plot_a <- plot_base +
    xlim(0, 1.5*max(na1, na2, nb1, nb2)) +
    geom_ribbon(data = rarefaction_table_unequal %>% filter(type == "rcurve"),
                aes(x=reads, ymin=richness, ymax=richness_max, col = ecosystem, fill = ecosystem), alpha = 0, linetype = 0) +
    geom_point(data = rarefaction_table_unequal %>% filter(type == "observation"),
               aes(x=reads, y=richness, col = ecosystem))+
    geom_line(data = rarefaction_table_unequal %>% filter(type == "rarefied_sample1"),
              aes(x=reads, y=richness, col = ecosystem)) +
    geom_line(data = rarefaction_table_unequal %>% filter(type == "rarefied_sample2"),
              aes(x=reads, y=richness, col = ecosystem)) +
    theme(legend.position = c(0.7, 0.2), legend.text=element_text(size=7))  +
    scale_color_discrete(name="",
                         breaks=c("ecosystem1", "ecosystem2"),
                         labels=c("Environment A", "Environment B")) +
    guides(fill=FALSE)+
    geom_text(data = data.frame("x"=(nb1), "y"= max(Ca, Cb)*1.05), aes(x,y),
              col = "black", label = "True taxonomic richnesses", cex = 3) +
    ggtitle("Rarefaction curve for alpha diversity")
  
  plot_b <- plot_base +
    geom_point(data = rarefaction_table_unequal %>% filter(type == "observation"),
               aes(x=ecosystem_number+c(-0.1, 0.1), y=richness, col = ecosystem)) +
    scale_x_continuous(name="Environment", breaks=c(1,2), labels=c("A", "B"), limits = c(0.5, 2.5)) +
    theme(legend.position="none") +
    ggtitle("Raw alpha diversity")
  
  
  plot_c <- plot_base +
    scale_x_continuous(name="Environment", breaks=c(1,2), labels=c("A", "B"), limits = c(0.5, 2.5)) +
    geom_point(data = rarefaction_table_unequal %>% filter(type %in% c("rarefied_sample1", "rarefied_sample2"))  %>%
                 group_by(ecosystem, ecosystem_number, type)  %>%
                 filter(reads <= rarefaction_level) %>%
                 filter(reads == max(reads)),
               aes(x=ecosystem_number+c(-0.1, 0.1), y=richness, col = ecosystem))   +
    theme(legend.position="none") +
    ggtitle("Rarefield alpha diversity")
  
  plot_d <- plot_base +
    scale_x_continuous(name="Environment", breaks=c(1,2), labels=c("A", "B"), limits = c(0.5, 2.5)) +
    geom_linerange(data = rarefaction_table_unequal %>% filter(type == "estimate1"),
                   aes(x = ecosystem_number-0.1, ymin=richness, ymax = richness_max), lty = 1) +
    geom_linerange(data = rarefaction_table_unequal %>% filter(type == "estimate2"),
                   aes(x = ecosystem_number+0.1, ymin=richness, ymax = richness_max), lty = 1) +
    theme(legend.position="none") +
    ggtitle("Bias + variance correction")
  
  grid.arrange(grobs = list(plot_a, plot_b, plot_c, plot_d),
               layout_matrix = matrix(c(1,1,2,3,4), nrow = 1))
}

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
  error <- matrix(NA, nrow = 1, ncol = 3)
  
  pvalues <- matrix(NA, nrow = 3, ncol = repeats)
  correct <- matrix(TRUE, nrow = 3, ncol = repeats)
  
  for (i in 1:repeats) {
    
    # generate samples
    rr <- replicates
    read_counts <- round(c(seq(from = na1, to = na2, length.out = rr),
                           seq(from = nb1, to = nb2, length.out = rr)))
    samples <- list()
    for (k in 1:(2*rr)) {
      my_name <- get(paste("population", ifelse(k <= rr, "A", "B"), sep=""))
      my_sample <- sample(1:length(my_name$proportions),
                          size = read_counts[k], prob = my_name$proportions, replace = T)
      samples[[k]] <- my_sample
    }
    
    # estimate richness
    cs <- lapply(samples, sample_richness) %>% unlist
    cs_rarefied <- lapply(samples, rarefied_richness, level = min(read_counts))  %>% unlist
    cc2 <- lapply(samples, estimate_richness_breakaway) %>% data.frame
    covariate <- rep(c("A", "B"), each = rr)
    cs_summary <- lm(cs ~ covariate) %>% summary
    cs_r_summary <- lm(cs_rarefied ~ covariate) %>% summary
    
    pvalues[1, i] <- cs_summary$coef[2,4]
    pvalues[2, i] <- cs_r_summary$coef[2,4]
    
    bta <- betta(cc2[1, ], cc2[2, ], X = cbind("Intercept" = 1, "CovariateB" = c(rep(0, rr), rep(1, rr) )))
    pvalues[3, i] <- bta$table[2,3]
    
    if (Ca != Cb) {
      correct[1, i] <- sign(cs_summary$coef[2,1]) == sign(Cb - Ca)
      correct[2, i] <- sign(cs_r_summary$coef[2,1]) == sign(Cb - Ca)
      
      ## TODO: is this right??
      correct[3, i] <- sign(bta$table[1,3]) == sign(Cb - Ca)
    }
    
  }
  
  if(Ca == Cb) {
    threshold <- pvalues > 0.05
  } else {
    threshold <- pvalues < 0.05
  }
  
  overall_correct <- threshold & correct
  error <- apply(overall_correct, 1, mean)
  names(error) <- c("Raw", "Rarefied", "Corrected")
  t(error)
  
}
