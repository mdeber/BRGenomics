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

df <- data.frame(pr = counts_pr, gb = counts_gb)

test_that("names given in input are used in output", {
    bin_deciles <- binNdimensions(df, nbins = 10)
    expect_equivalent(names(bin_deciles), c("bin.pr", "bin.gb"))
})

pidx <- counts_pr / counts_gb

test_that("can use different nbins for different dimensions", {
    df3 <- data.frame(counts_pr, counts_gb, pidx)

    expect_is(binNdimensions(df3), "data.frame")
    expect_error(binNdimensions(df3, nbins = c(10, 15)))
    bin_3d <- binNdimensions(df3, nbins = c(10, 15, 20))
    expect_is(bin_3d, "data.frame")
    expect_equivalent(sapply(bin_3d, max), c(10, 15, 20))
})


# Test aggregation --------------------------------------------------------


test_that("can aggregate in ndimensional bins", {
    mean_bins <- aggregateByNdimensionalBins(seq_len(nrow(df)), df, ncores = 2)

    expect_is(mean_bins, "data.frame")
    expect_equal(names(mean_bins), c("bin.pr", "bin.gb", "value"))
    expect_equal(round(mean_bins$value[1]), 139)
    expect_equal(sum(is.na(mean_bins$value)), 21)
})


test_that("can modify function and default output values", {
    sum_bins <- aggregateByNdimensionalBins(seq_len(nrow(df)), df, FUN = sum,
                                            ncores = 2)
    expect_equal(sum_bins$value[1], 18540)
    expect_true(is.na(sum_bins$value[3]))

    sum_bins2 <- aggregateByNdimensionalBins(seq_len(nrow(df)), df, FUN = sum,
                                             empty = 0, ncores = 2)
    expect_equal(sum_bins2$value[1], 18540)
    expect_equal(sum_bins2$value[3], 0)
})


# Test density ------------------------------------------------------------

test_that("can find density in ndimensional bins", {
    bin_density <- densityInNdimensionalBins(df, ncores = 2)
    expect_equivalent(bin_density$value[1:3], c(133, 9, 0))
})


