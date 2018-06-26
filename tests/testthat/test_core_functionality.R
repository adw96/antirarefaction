library(antirarefaction)
context("Test output")

test_that("make_community behaves OK", {
  
  expect_equal(make_community(1400, 2, 2) %>% length, 1400)
  expect_equal(make_community(700, 1, 0.5) %>% sum, 1)
  
})

test_that("sample_from_population behaves OK", {
  popn <- make_community(10, 2, 2)
  expect_equal(sample_from_population(popn, 4) %>% length, 4)
})

test_that("sample_richness behaves OK", {
  
  expect_equal(sample_richness(paste("species", 1:10)), 10)
  
})

test_that("get_point and draw_rarefaction behave OK", {
  
  expect_equal(sample_from_population(popn, 10) %>% length, 10)
  
  expect_is(draw_rarefaction(population = popn, pts = 7), 
            "data.frame")
  expect_is(draw_rarefaction(population = popn, pts = 7, replicates = 1), 
            "data.frame")
  expect_is(draw_rarefaction(population = popn, pts = 7, replicates = 2), 
            "data.frame")
  expect_is(draw_rarefaction(population = popn, pts = 3:10, replicates = 2), 
            "data.frame")
  
  expect_is(draw_rarefaction(population = popn, pts = c(5,7), replicates = 2), 
            "data.frame")
  
  popn <- make_community(2000, 2, 2)
  expect_is(draw_rarefaction(population = popn), 
            "data.frame")
  
  popn <- make_community(20, 2, 2)
  sample_from_population(popn)
  get_point(popn, 
            read = 10, label = 1, ecosystem_name = "popn")
  
  popn <- make_community(200, 2, 2)
  sample_from_population(popn)
  get_point(popn, 
            read = 10, label = 1, ecosystem_name = "popn")
  
})


test_that("make_table behaves ok", {
  
  mt <- make_table(Ca = 1400, 
                   Cb = 1400, 
                   aa = 2, 
                   ba = 1.8, 
                   ab = 2, 
                   bb = 2, 
                   na1 = 7000, 
                   nb1 = 7000, 
                   na2 = 7000, 
                   nb2 = 700)
  make_plot(mt, 
            Ca = 1400, 
            Cb = 1400, 
            aa = 2, 
            ba = 1.8, 
            ab = 2, 
            bb = 2, 
            na1 = 7000, 
            nb1 = 7000, 
            na2 = 7000, 
            nb2 = 700)
})


test_that("get_error_rates", {
  
  ger <- get_error_rates(repeats = 5, 
                         replicates = 10,
                         Ca = 1400, 
                         Cb = 1400, 
                         aa = 2, 
                         ba = 1.8, 
                         ab = 2, 
                         bb = 2, 
                         na1 = 7000, 
                         nb1 = 7000, 
                         na2 = 7000, 
                         nb2 = 700) 
  
  ger
  
  
  expect_is(ger, "data.frame")
  
  ## construct a setting where raw should be fine
  
  ## construct a setting where rarefied should be fine
  
  ## construct a setting where betta should be best
  
  
})

Ca = 1400 
Cb = 1400 
aa = 2 
ba = 1.8 
ab = 2 
bb = 2 
na1 = 7000 
nb1 = 7000 
na2 = 7000 
nb2 = 700
