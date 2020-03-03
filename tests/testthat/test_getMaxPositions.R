context("Finding max-signal positions within regions of interest")
library(BRGenomics)

data(PROseq)
data("txs_dm6_chr4")

test_that("Empty input regions returns empty GRanges", {
    out <- getMaxPositionsBySignal(dataset.gr = PROseq,
                                   regions.gr = txs_dm6_chr4[0])

    expect_is(out, "GRanges")
    expect_equal(length(out), 0)
    expect_equivalent(names(mcols(out)), names(mcols(txs_dm6_chr4)))
})

test_that("keep.signal with empty input regions is properly formatted", {
    # should add metadata for MaxSiteSignal
    out <- getMaxPositionsBySignal(dataset.gr = PROseq,
                                   regions.gr = txs_dm6_chr4[0],
                                   keep.signal = TRUE)
    expect_equal(ncol(mcols(out)), ncol(mcols(txs_dm6_chr4)) + 1)
})

test_that("No signal overlapping regions returns formatted, empty GRanges", {
    # 3 regions without overlapping signal
    out <- getMaxPositionsBySignal(dataset.gr = PROseq,
                                   regions.gr = txs_dm6_chr4[c(24, 60, 62)])
    expect_is(out, "GRanges")
    expect_equal(length(out), 0)
    expect_equivalent(names(mcols(out)), names(mcols(txs_dm6_chr4)))
})

test_regions <- unique(txs_dm6_chr4)[1:10]
test_regions_promoter <- promoters(test_regions, 0, 300)

test_that("Max sites found for regions with signal", {
    out <- getMaxPositionsBySignal(dataset.gr = PROseq,
                                   regions.gr = test_regions_promoter)
    expect_equal(length(out), length(subsetByOverlaps(test_regions_promoter,
                                                      PROseq)))
    expect_equivalent(start(out)[1:4], c(42798, 44800, 69401, 69401))
    expect_true(all(width(out) == 1))
})

test_that("Correct signal metadata added for regions with signal", {
    out <- getMaxPositionsBySignal(dataset.gr = PROseq,
                                   regions.gr = test_regions_promoter,
                                   keep.signal = TRUE)
    # two repeated tests
    expect_equivalent(start(out)[1:4], c(42798, 44800, 69401, 69401))
    expect_equal(ncol(mcols(out)), ncol(mcols(txs_dm6_chr4)) + 1)

    score_col <- which(!names(mcols(out)) %in% names(mcols(txs_dm6_chr4)))
    expect_equivalent(mcols(out)[1:4, score_col], c(3, 1, 122, 122))
})

test_that("Max sites found with multi-width input", {
    out <- getMaxPositionsBySignal(dataset.gr = PROseq,
                                   regions.gr = test_regions)
    expect_is(out, "GRanges")
    expect_equal(length(out), length(test_regions)) # all have signal
    expect_equivalent(start(out)[1:4], c(1295, 42798, 45915, 60309))
})

test_that("Max signals found in larger bins", {
    out <- getMaxPositionsBySignal(dataset.gr = PROseq,
                                   regions.gr = test_regions,
                                   binsize = 10, bin.centers = FALSE)
    expect_equal(length(out), length(test_regions)) # all have signal
    expect_equivalent(start(out)[1:4], c(1289, 42794, 45914, 60907))
})

test_that("Max signals found in larger bins, bin.centers = TRUE", {
    out <- getMaxPositionsBySignal(dataset.gr = PROseq,
                                   regions.gr = test_regions,
                                   binsize = 10, bin.centers = TRUE)
    expect_equal(length(out), length(test_regions)) # all have signal
    expect_equivalent(start(out)[1:4], c(1293, 42798, 45918, 60911))
    out2 <- getMaxPositionsBySignal(dataset.gr = PROseq,
                                    regions.gr = test_regions,
                                    binsize = 11, bin.centers = TRUE)
    expect_equivalent(start(out2)[1:4], c(1291, 43285, 45912, 60913))
})





