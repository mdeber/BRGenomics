context("counting signal by region")
library(BRGenomics)

data("PROseq")
data("txs_dm6_chr4")

test_counts <- getCountsByRegions(PROseq, txs_dm6_chr4)

test_that("signal counts correct returns correct vector", {
    expect_is(test_counts, "numeric")
    expect_equal(length(test_counts), length(txs_dm6_chr4))
    expect_equivalent(test_counts[1:4], c(1, 58, 13, 125))
})

ps_rename <- PROseq
names(mcols(ps_rename)) <- "signal"

test_that("signal counts works with alternative metadata fields", {
    expect_error(getCountsByRegions(ps_rename, txs_dm6_chr4))
    expect_equivalent(test_counts, getCountsByRegions(ps_rename,
                                                      txs_dm6_chr4,
                                                      field = "signal"))
})

test_that("signal counts work over multiple metadata fields", {
    ps_rename$posnum <- seq_along(ps_rename)
    counts_multiple <- getCountsByRegions(ps_rename, txs_dm6_chr4,
                                          field = c("signal", "posnum"))
    expect_is(counts_multiple, "data.frame")
    expect_equivalent(counts_multiple$signal, test_counts)
})



# getCountsByPositions ----------------------------------------------------
context("Counting signal at positions within regions")




# Calculating pause indices -----------------------------------------------
context("Calculating pausing indices")


