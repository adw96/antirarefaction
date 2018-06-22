### How to use this package

library(devtools)
build("/Users/amy/Documents/software/antirarefaction/")
library(roxygen2)
roxygenise("/Users/amy/Documents/software/antirarefaction/")
install("/Users/amy/Documents/software/antirarefaction/", quick = T)
library(antirarefaction)
devtools::test()

