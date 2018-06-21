library(antirarefaction)
context("Test output")

# set.seed(1)
# n <- 6
# n_taxa <- 5
# my_counts <- matrix(rpois(n*n_taxa, lambda=10), nrow = n)
# my_discrete_covariate <- cbind(1, rep(c(0,1), each = n/2), rep(c(0,1), n/2))
# my_continuous_covariate <- rnorm(n)
# 
# test_that("errors are thrown", {
#   expect_error(divnet(matrix(c(10, 20, 10, 1), nrow= 1), tuning="test"))  
#   expect_error(divnet(matrix(c(10, 20, 10, 1), ncol=2), tuning="test"))
#   expect_error(divnet(matrix(c(10, 20, 10, 1, 50, 0), ncol=2), tuning="test"))
#   expect_is(divnet(matrix(c(10, 20, 11, 1, 50, 1), nrow=2), 
#                    tuning="test"), "diversityEstimates")
# })
