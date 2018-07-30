library(antirarefaction)
context("Test output")


test_that("make_community behaves OK", {
  
  # make_community
  expect_silent(popn <- make_community(20, 2, 2))
  expect_equal(popn %>% nrow, 20)
  expect_equal(popn$abundances %>% sum, 1)
  
  # sample_from_population
  expect_equal(sample_from_population(popn, 4) %>% nrow, 4)
  
  # breakaway::sample_richness
  expect_silent(smpl <- sample_from_population(popn, 500))
  expect_equal(richness_from_sample(smpl), 
               20)
  
  expect_lt(rarefied_sample(smpl, 10) %>% richness_from_sample, richness_from_sample(smpl))
  
})

test_that("plotting", {
  Ca = 1400
  Cb = 1000
  aa = 2
  ba = 1
  ab = 2
  bb = 2
  na = 3000
  nb = 7000
  
  popnA <- make_community(Ca, aa, ab)
  popnB <- make_community(Cb, ba, bb)
  
  smplA <- sample_from_population(popnA, na)
  smplB <- sample_from_population(popnB, nb)
  
  popnA %>% community_richnness
  
  truth <- data.frame("sample" = c("A", "B"),
                      "richness" = c(popnA %>% community_richnness, 
                                     popnA %>% community_richnness),
                      "depth" = NA,
                      "type" = "true")
  truth
  
  obsd <- data.frame("sample" = c("A", "B"),
                     "richness" = c(smplA %>% richness_from_sample, 
                                    smplB %>% richness_from_sample),
                     "depth" = c(na, nb),
                     "type" = "observed")
  
  # subsample from larger sample to level of smaller sample
  rs <- rarefied_sample(smpl = if(na < nb) { smplB } else {smplA},
                        level = ifelse(na < nb, na, nb))
  
  rrd <- data.frame("sample" = ifelse(na < nb, "B", "A"),
                    "richness" = rs %>% richness_from_sample, 
                    "depth" = ifelse(na < nb, na, nb),
                    "type" = "rarefied")
  
  ## rarefy
  
  # nlower <- ifelse(na < nb, na, nb) + 1
  # nupper <- ifelse(na < nb, nb, na) - 1
  # 
  # rrds <- sapply(nlower:nupper, function(n) {
  #   rarefied_sample(smpl = if(na < nb) { smplB } else {smplA},
  #                   level = n) %>% richness_from_sample
  # })
  # 
  # rrds_df <- data.frame("sample" = ifelse(na < nb, "B", "A"),
  #                       "richness" = rrds, 
  #                       "depth" = nlower:nupper,
  #                       "type" = "rarefied_seq")
  # 
  # all <- rbind(truth, obsd, rrd, rrds_df)
  
  all <- rbind(truth, obsd, rrd)
  
  gg <- all %>%
    ggplot(aes(x = depth, y = richness, col = sample)) +
    theme_bw() +
    geom_point(data = all %>% dplyr::filter(type == "observed")) +
    geom_point(data = all %>% dplyr::filter(type == "rarefied"), pch = 2) +
    # geom_point(data = all %>% dplyr::filter(type == "rarefied_seq"), cex = 0.2, alpha = 0.2) +
    geom_hline(data = all %>% dplyr::filter(type == "true"),
               aes_string(yintercept="richness", lty = "sample", col = "sample"))  +
    scale_linetype_manual(guide = "none", values=c("twodash", "dotted")) +
    labs(x="No. reads", y="Taxonomic richness",col = "") +
    coord_cartesian(ylim=c(0, max(Ca, Cb)*1.1)) +
    xlim(0, 1.5*max(na, nb))
  
  expect_is(gg, "ggplot")
  
})


# 
# test_that("get_point behave OK", {
#   popn <- make_community(20, 2, 2)
#   
#   expect_is(get_point(popn, 
#                       read = 10, label = 1, ecosystem_name = "popn"), 
#             "data.frame")
#   
# })
# 
# test_that("draw_rarefaction behave OK", {
#   popn <- make_community(10, 2, 2)
#   
#   expect_is(draw_rarefaction(population = popn, pts = 7), 
#             "data.frame")
#   
#   expect_is(draw_rarefaction(population = popn, pts = 7, replicates = 1), 
#             "data.frame")
#   
#   expect_is(draw_rarefaction(population = popn, pts = 7, replicates = 2), 
#             "data.frame")
#   expect_is(draw_rarefaction(population = popn, pts = 3:10, replicates = 2), 
#             "data.frame")
#   
#   expect_is(draw_rarefaction(population = popn, pts = c(5,7), replicates = 2), 
#             "data.frame")
#   
# })
# 
# 
# test_that("make_table behaves ok", {
#   
#   expect_is(mt <- make_table(Ca = 140, 
#                              Cb = 140, 
#                              aa = 2, 
#                              ba = 1.8, 
#                              ab = 2, 
#                              bb = 2, 
#                              na1 = 70, 
#                              nb1 = 80, 
#                              na2 = 10, 
#                              nb2 = 20), 
#             "data.frame")
#   
#   
#   expect_is(plot_a(mt, 
#                    Ca = 140, 
#                    Cb = 140, 
#                    aa = 2, 
#                    ba = 1.8, 
#                    ab = 2, 
#                    bb = 2, 
#                    na1 = 70, 
#                    nb1 = 80, 
#                    na2 = 10, 
#                    nb2 = 20), 
#             "ggplot")
# })
# 
# 
# test_that("get_error_rates", {
#   
#   expect_is(get_error_rates(repeats = 5, 
#                             replicates = 10,
#                             Ca = 1400, 
#                             Cb = 1400, 
#                             aa = 2, 
#                             ba = 1.8, 
#                             ab = 2, 
#                             bb = 2, 
#                             na1 = 7000, 
#                             nb1 = 7000, 
#                             na2 = 7000, 
#                             nb2 = 700) , "data.frame")
#   
#   ## construct a setting where raw should be fine
#   
#   ## construct a setting where rarefied should be fine
#   
#   ## construct a setting where betta should be best
#   
#   
# })

