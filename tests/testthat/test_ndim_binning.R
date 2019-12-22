context("N-dimensional binning")
library(BRGenomics)

data("PROseq")
data("txs_dm6_chr4")

txs_pr <- promoters(txs_dm6_chr4, 0, 100)
txs_gb <- flank(txs_pr, 300, start = FALSE)

counts_pr <- getCountsByRegions(PROseq, txs_pr)
counts_gb <- getCountsByRegions(PROseq, txs_gb)

idx <- union(which(counts_pr > 0), which(counts_gb > 0))
counts_pr <- counts_pr[idx]
counts_gb <- counts_gb[idx]

bin_deciles <- binNdimensions(counts_pr, counts_gb, quantiles = 10)

test_that("ndim binning on numeric vectors", {
    expect_is(bin_deciles, "data.frame")
    expect_equivalent(names(bin_deciles), c("bin.counts_pr",
                                            "bin.counts_gb"))
})

test_that("names given in function call are used in output", {
    bin_alt <- binNdimensions(pr = counts_pr, gb = counts_gb,
                              quantiles = 10)
    expect_equivalent(names(bin_alt), c("bin.pr", "bin.gb"))
})

test_that("ndim binning the same on a dataframe", {
    expect_equivalent(bin_deciles,
                      binNdimensions(data.frame(counts_pr, counts_gb),
                                     quantiles = 10))
})

test_that("error if >1 unnamed args given when dataframe given", {
    expect_error(binNdimensions(data.frame(counts_pr, counts_gb),
                                counts_gb, quantiles = 10))
})

test_that("ndim binning on list, and mixed list/vector inputs", {
    expect_equivalent(bin_deciles,
                      binNdimensions(list(counts_pr), list(counts_gb),
                                     quantiles = 10))
    expect_equivalent(bin_deciles,
                      binNdimensions(list(counts_pr), counts_gb,
                                     quantiles = 10))
})

pidx <- counts_pr / counts_gb

test_that("can use different quantiles for different dimensions", {
    expect_is(binNdimensions(counts_pr, counts_gb, pidx), "data.frame")
    expect_error(binNdimensions(counts_pr, counts_gb, pidx,
                                quantiles = c(10, 15)))
    bin_3d <- binNdimensions(counts_pr, counts_gb, pidx,
                             quantiles = c(10, 15, 20))
    expect_is(bin_3d, "data.frame")
    expect_equivalent(sapply(bin_3d, max), c(10, 15, 20))
})

test_that("messages produced on attempts to coerce unexpected classes", {
    expect_warning(binNdimensions(cbind(counts_pr, counts_gb)))
    expect_error(binNdimensions(rle(counts_pr), rle(counts_gb)))
    expect_warning(binNdimensions(counts_pr = Rle(counts_pr),
                                  counts_gb = Rle(counts_pr), quantiles = 10))
})




