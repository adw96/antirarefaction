library(antirarefaction)
context("Test output")

test_that("make_community", {
  
  expect_equal(make_community(1400, 2, 2) %>% length, 1400)
  expect_equal(make_community(700, 1, 0.5) %>% sum, 1)
  
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
  
  
  expect_is(ger, "data.frame")

  ## construct a setting where raw should be fine

  ## construct a setting where rarefied should be fine
  
  ## construct a setting where betta should be best
  
  
})

