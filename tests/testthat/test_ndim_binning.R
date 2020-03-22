context("N-dimensional binning")
library(BRGenomics)


# Test binning ------------------------------------------------------------

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

bin_deciles <- binNdimensions(df, nbins = 10, ncores = 1)

test_that("names given in input are used in output", {
    expect_equivalent(names(bin_deciles), c("bin.pr", "bin.gb"))
})

pidx <- counts_pr / counts_gb

test_that("can use different nbins for different dimensions", {
    df3 <- data.frame(counts_pr, counts_gb, pidx)

    expect_is(binNdimensions(df3, ncores = 1), "data.frame")
    expect_error(binNdimensions(df3, nbins = c(10, 15), ncores = 1))
    bin_3d <- binNdimensions(df3, nbins = c(10, 15, 20), ncores = 1)
    expect_is(bin_3d, "data.frame")
    expect_equivalent(sapply(bin_3d, max), c(10, 15, 21))
})

test_that("use bin centers, not numbers", {
    bin_deciles <- binNdimensions(df, nbins = 10, ncores = 1)
    bindec_centers <- binNdimensions(df, nbins = 10, use_bin_numbers = FALSE,
                                     ncores = 1)

    # remake bins in test, and use bin_deciles to index them
    # check for bin.pr
    bins.pr <- seq(min(df$pr), max(df$pr), length.out = 10+1) # make bins
    centers.pr <- sapply(1:( length(bins.pr)-1 ),
                         function(i) median(bins.pr[i:(i+1)]))

    expect_equivalent(bindec_centers$bin.pr,
                      centers.pr[bin_deciles$bin.pr])

    # check for bin.gb
    bins.gb <- seq(min(df$gb), max(df$gb), length.out = 10+1) # make bins
    centers.gb <- sapply(1:( length(bins.gb)-1 ),
                         function(i) median(bins.gb[i:(i+1)]))

    expect_equivalent(bindec_centers$bin.gb,
                      centers.gb[bin_deciles$bin.gb])
})


# Test aggregation --------------------------------------------------------


test_that("can aggregate in ndimensional bins", {
    mean_bins <- aggregateByNdimensionalBins(seq_len(nrow(df)), df, ncores = 1)

    expect_is(mean_bins, "data.frame")
    expect_equal(names(mean_bins), c("bin.pr", "bin.gb", "value"))
    expect_equal(round(mean_bins$value[1]), 141)
    expect_equal(sum(is.na(mean_bins$value)), 21)

    # test using x = column name
    dfplus <- df
    dfplus$value <- seq_len(nrow(df))
    mbinsp <- aggregateByNdimensionalBins("value", dfplus, ncores = 1)
    expect_identical(mean_bins, mbinsp)
})

test_that("can modify function and default output values", {
    sum_bins <- aggregateByNdimensionalBins(seq_len(nrow(df)), df, FUN = sum,
                                            ncores = 1)
    expect_equal(sum_bins$value[1], 18521)
    expect_true(is.na(sum_bins$value[4]))

    sum_bins2 <- aggregateByNdimensionalBins(seq_len(nrow(df)), df, FUN = sum,
                                             empty = 0, ncores = 1)
    expect_equal(sum_bins2$value[1], 18521)
    expect_equal(sum_bins2$value[4], 0)
})

test_that("can separate 'empty NAs' from 'output NAs'", {
    # the key here is not dropping empty bins; not setting empty bins to NA; and
    # not ignoring NA values returned from FUN -> there's unique method dispatch
    # here because we need to keep NAs returned by aggregate
    sd_bins <- aggregateByNdimensionalBins(seq_len(nrow(df)), df, FUN = sd,
                                           ignore.na = FALSE, empty = 3.14159,
                                           ncores = 1)
    sd_bins2 <- aggregateByNdimensionalBins(seq_len(nrow(df)), df, FUN = sd,
                                            ignore.na = FALSE, empty = 0,
                                            ncores = 1)
    idx_empty <- which(sd_bins$value == 3.14159)
    expect_true(all(sd_bins2$value[idx_empty] == 0))
    expect_equivalent(sd_bins$value[-idx_empty], sd_bins2$value[-idx_empty])
})

# Test density ------------------------------------------------------------

test_that("can find density in ndimensional bins", {
    bin_density <- densityInNdimensionalBins(df, ncores = 1)
    expect_equivalent(bin_density$value[1:3], c(131, 6, 3))
})
