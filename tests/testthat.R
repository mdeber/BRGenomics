Sys.setenv("R_TESTS" = "")
library(testthat)
library(BRGenomics)

test_check("BRGenomics")
