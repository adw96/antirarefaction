#' Subsample species to draw a rarefaction curve
#' 
#' @param population The population to draw from
#' @param pts A vector of depths to subsample to
#' @param replicates The number of times to repeat (for confidence intervals)
#' @param ecosystem The ecosystem name (useful for comparing ecosystems in plots)
#' 
#' @import magrittr
#' @import dplyr
#' 
#' @return A data frame with each row corresponding to a sampling depth (number of rows = length of pts)
#' 
#' @export
draw_rarefaction <- function(population, pts = NA, replicates = 5, ecosystem = NA) {
  rarefaction_table <- data.frame()
  
  CC <- length(population)
  pts %<>% round
  pts %<>% unique
  
  if (all(is.na(pts))) pts <- c(round(seq(from = 5, to = CC/100, length = 100)),
                                round(seq(from = CC/100, to = CC/10, length = 100)),
                                round(seq(from = CC/10, to = CC/4, length = 100)),
                                round(seq(from = CC/4, to = CC, length = 10)))
  who <- replicate(replicates, 
                   sample_from_population(population, 
                                          max(pts)))
  
  # terrible way of doing this
  # maybe not so bad, since need nested
  # definitely could be faster
  # TODO fix up use of shannon
  for (i in pts) {
    subsample_diversity <- rep(NA, replicates)
    subsample_shannon <- rep(NA, replicates)
    for (j in 1:replicates) {
      subsample <- who[1:i, j]
      subsample_diversity[j] <- sample_richness(subsample)
      community <- c(as.matrix(table(subsample)))
      proportions <- community/sum(community)
      subsample_shannon[j] <- breakaway::shannon(proportions)
      
      rarefaction_table %<>% rbind(data.frame("reads" = i,
                                              "richness" = min(subsample_diversity), 
                                              "richness_max" = max(subsample_diversity),
                                              "shannon" = min(subsample_shannon), 
                                              "shannon_max" = max(subsample_shannon),
                                              "ecosystem" = ecosystem, 
                                              "type" = "rcurve"))
    }
  }
  rarefaction_table
}

#' Points for plotting a rarefaction curved based on a sample
#' 
#' @param my_population The population distribution
#' @param read Sample depth
#' @param label label
#' @param ecosystem_name ecosystem_name
#' 
#' @import magrittr
#' 
#' @export
get_point <- function(my_population, read, label, ecosystem_name, 
                      replicates = 1,
                      pts) {
  
  read %<>% round
  subsample <- sample_from_population(my_population, read)
  
  point <- data.frame("reads" = read,
                      "richness" = sample_richness(subsample), 
                      "richness_max" = NA,
                      "shannon" = breakaway::shannon(c(as.matrix(table(subsample)))/sum(c(as.matrix(table(subsample))))), 
                      "shannon_max" = NA,
                      "ecosystem" = ecosystem_name, 
                      "type" = "observation")
  
  subsample1 <- subsample %>% table %>% to_proportions(type = "column")
  ### draw_rarefaction props, not labels
  rarefied_sample <- draw_rarefaction(subsample1, 
                                      replicates = replicates, 
                                      pts = pts,
                                      ecosystem = ecosystem_name)
  rarefied_sample$type <- paste("rarefied_sample", label, sep="")
  xx <- subsample %>% table %>% make_frequency_count_table
  
  estimate <- xx %>% breakaway(answers = T, plot = F, output=F) %$% est
  error <- xx %>% breakaway(answers = T, plot = F, output=F) %$% seest
  richness_estimate <- c(estimate - 2*error, estimate + 2*error)
  
  estimate <- data.frame("reads" = read,
                         "richness" = richness_estimate[1], 
                         "richness_max" = richness_estimate[2],
                         "shannon" = 0, 
                         "shannon_max" = 100000,
                         "ecosystem" = ecosystem_name, 
                         "type" = paste("estimate", label, sep=""))
  rbind(point, rarefied_sample, estimate)
}


#' make_plot makes plots
#' 
#' @import dplyr
#' @import ggplot2
#' 
#' @export
plot_a <- function(rarefaction_table_unequal,
                   Ca, Cb, aa, ab, ba, bb, na1, na2, nb1, nb2) {
  ## TODO: define Ca, Cb, na1, na2...
  plot_base <- ggplot(data = rarefaction_table_unequal, 
                      aes_string(col = "ecosystem")) +
    theme_bw() +
    geom_hline(data = rarefaction_table_unequal %>% 
                 filter(type == "truth"),
               aes_string(yintercept="richness", 
                          col = "ecosystem", lty = "ecosystem")) +
    scale_linetype_manual(guide = "none", values=c("twodash", "dotted")) +
    labs(x="No. reads", y="Taxonomic richness",col = "") +
    coord_cartesian(ylim=c(0, max(Ca, Cb)*1.1))
  
  plot_a <- plot_base +
    xlim(0, 1.5*max(na1, na2, nb1, nb2)) +
    geom_ribbon(data = rarefaction_table_unequal %>% filter(type == "rcurve"),
                aes_string(x="reads", ymin="richness", ymax="richness_max", 
                           col = "ecosystem", fill = "ecosystem"), alpha = 0, linetype = 0) +
    geom_point(data = rarefaction_table_unequal %>% filter(type == "observation"),
               aes_string(x="reads", y="richness", col = "ecosystem"))+
    geom_line(data = rarefaction_table_unequal %>% filter(type == "rarefied_sample1"),
              aes_string(x="reads", y="richness", col = "ecosystem")) +
    geom_line(data = rarefaction_table_unequal %>% filter(type == "rarefied_sample2"),
              aes_string(x="reads", y="richness", col = "ecosystem")) +
    theme(legend.position = c(0.7, 0.2), legend.text=element_text(size=7))  +
    scale_color_discrete(name="",
                         breaks=c("ecosystem1", "ecosystem2"),
                         labels=c("Environment A", "Environment B")) +
    guides(fill=FALSE)+
    # geom_text(data = data.frame("x"=(nb1), 
    #                             "y"= max(Ca, Cb)*1.05), 
    #           aes_string("x", "y"),
    #           col = "black", label = "True taxonomic richnesses", cex = 3) +
    # ggtitle("Rarefaction curve for alpha diversity") +
    NULL
  
  plot_a
  
}

plot_b <- function(rarefaction_table_unequal,
                   Ca, Cb, aa, ab, ba, bb, na1, na2, nb1, nb2) {
  ## TODO: define Ca, Cb, na1, na2...
  plot_base <- ggplot(data = rarefaction_table_unequal, 
                      aes_string(col = "ecosystem")) +
    theme_bw() +
    geom_hline(data = rarefaction_table_unequal %>% 
                 filter(type == "truth"),
               aes_string(yintercept="richness", col = "ecosystem", lty = "ecosystem")) +
    scale_linetype_manual(guide = "none", values=c("twodash", "dotted")) +
    labs(x="No. reads", y="Taxonomic richness", col = "") +
    coord_cartesian(ylim=c(0, max(Ca, Cb)*1.1))
  
  plot_base +
    geom_point(data = rarefaction_table_unequal %>% filter(type == "observation"),
               aes(x=ecosystem_number+c(-0.1, 0.1), y=richness, col = ecosystem)) +
    scale_x_continuous(name="Environment", breaks=c(1,2), labels=c("A", "B"), limits = c(0.5, 2.5)) +
    theme(legend.position="none") +
    ggtitle("Raw alpha diversity")
  
}

#' plot_c makes plots
#' 
#' @import dplyr
#' @import ggplot2
#' 
#' @export
plot_c <- function() {
  
  plot_base +
    scale_x_continuous(name="Environment", breaks=c(1,2), labels=c("A", "B"), limits = c(0.5, 2.5)) +
    geom_point(data = rarefaction_table_unequal %>% filter(type %in% c("rarefied_sample1", "rarefied_sample2"))  %>%
                 group_by(ecosystem, ecosystem_number, type)  %>%
                 filter(reads <= rarefaction_level) %>%
                 filter(reads == max(reads)),
               aes_string(x=ecosystem_number+c(-0.1, 0.1), y=richness, col = ecosystem))   +
    theme(legend.position="none") +
    ggtitle("Rarefield alpha diversity")
}

#' plot_d makes plots
#' 
#' @import dplyr
#' @import ggplot2
#' 
#' @export
plot_d <- function() {
  plot_base +
    scale_x_continuous(name="Environment", breaks=c(1,2), labels=c("A", "B"), limits = c(0.5, 2.5)) +
    geom_linerange(data = rarefaction_table_unequal %>% filter(type == "estimate1"),
                   aes_string(x = ecosystem_number-0.1, ymin=richness, ymax = richness_max), lty = 1) +
    geom_linerange(data = rarefaction_table_unequal %>% filter(type == "estimate2"),
                   aes_string(x = ecosystem_number+0.1, ymin=richness, ymax = richness_max), lty = 1) +
    theme(legend.position="none") +
    ggtitle("Bias + variance correction")
  
}

