### How to use this package

library(antirarefaction)

Ca <- 1400
Cb <- 1400
na1 <- 7000
nb1 <- 7000
na2 <- 7000
nb2 <- 7000
aa <- 2
ba <- 1.8
ab <- 2
bb <- 2
replicates <- 2
repeats <- 5

get_error_rates(repeats, 
                replicates,
                Ca, Cb, 
                aa, ab, ba, bb, 
                na1, na2, nb1, nb2) 