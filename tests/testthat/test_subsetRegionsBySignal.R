context("Subsetting genelists by signal")
library(BRGenomics)

data("PROseq")
data("txs_dm6_chr4")

test_regions <- txs_dm6_chr4[1:100]

out_high <- subsetRegionsBySignal(test_regions,
                                  dataset.gr = PROseq,
                                  quantiles = c(0.5, 1))

out_low <- subsetRegionsBySignal(test_regions,
                                 dataset.gr = PROseq,
                                 quantiles = c(0, 0.5))

test_that("can output subset genelist", {
    expect_is(out_high, "GRanges")
    expect_is(out_low, "GRanges")
    expect_equivalent(names(mcols(test_regions)), names(mcols(out_high)))
    expect_equivalent(names(mcols(test_regions)), names(mcols(out_low)))
    expect_true(length(out_high) < length(test_regions))
    expect_true(length(out_low) < length(test_regions))
})

test_that("subsetted genelist has expected length", {
    expect_equal(0.5*length(test_regions), length(out_high))
    expect_equal(0.5*length(test_regions), length(out_low))

    out <- subsetRegionsBySignal(test_regions, PROseq, c(0.05, 0.95))
    expect_equal(0.9*length(test_regions), length(out))
})

test_that("subsetted genelist has expected breaks", {
    test_counts <- getCountsByRegions(PROseq, test_regions)
    counts_high <- getCountsByRegions(PROseq, out_high)
    counts_low <- getCountsByRegions(PROseq, out_low)

    expect_equal(max(test_counts), max(counts_high))
    expect_equal(min(test_counts), min(counts_low))
    expect_true(min(counts_high) >= max(counts_low))
})

test_that("subsetting by density", {
    density_high <- subsetRegionsBySignal(test_regions,
                                          dataset.gr = PROseq,
                                          quantiles = c(0.5, 1),
                                          density = TRUE)
    expect_false(all(width(out_high) == width(density_high)))
})

test_that("ordering by rank", {
    ordered_high <- subsetRegionsBySignal(test_regions,
                                          dataset.gr = PROseq,
                                          quantiles = c(0.5, 1),
                                          order.by.rank = TRUE)
    expect_true(all(order(out_high) == seq_along(ordered_high)))
    expect_false(all(order(ordered_high) == seq_along(ordered_high)))
})

test_that("subsetting from quantiles (0,1) returns full genelist", {
    out <- subsetRegionsBySignal(test_regions,
                                 dataset.gr = PROseq,
                                 quantiles = c(0, 1))
    expect_equal(length(out), length(test_regions))
    expect_equivalent(out, sort(test_regions))
})

test_that("lower_quantile = 1 | upper_quantile = 0 returns empty genelist", {
    out <- subsetRegionsBySignal(test_regions,
                                 dataset.gr = PROseq,
                                 quantiles = c(1, 1))
    expect_is(out, "GRanges")
    expect_equal(length(out), 0)

    out <- subsetRegionsBySignal(test_regions,
                                 dataset.gr = PROseq,
                                 quantiles = c(0, 0))
    expect_is(out, "GRanges")
    expect_equal(length(out), 0)
})

test_that("error returned if lower quantile not less than upper quantile", {
    expect_error(subsetRegionsBySignal(test_regions, PROseq, c(0.6, 0.4)))
})

test_that("error returned if quantiles outside (0,1)", {
    expect_error(subsetRegionsBySignal(test_regions, PROseq, c(0.5, 1.01)))
    expect_error(subsetRegionsBySignal(test_regions, PROseq, c(-0.01, 0.5)))
})


