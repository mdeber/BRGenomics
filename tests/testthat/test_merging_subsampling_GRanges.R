context("Random subsampling GRanges data")
library(BRGenomics)
data("PROseq")

set.seed(11)
ps_tenth <- subsampleGRanges(PROseq, prop = 0.1)

test_that("Subsampling produces GRanges", {
    expect_is(ps_tenth, "GRanges")
})

test_that("Subsampling gives correct signal count", {
    expect_equal(floor(sum(score(PROseq))), sum(score(ps_tenth)))
})

test_that("Subsamplied reproduced using set.seed", {
    expect_equal(length(ps_tenth), 6397)
})


# Merging GRanges ---------------------------------------------------------

context("Merging GRanges data")
ps_breaks <- seq(1, length(PROseq), floor(0.25*length(PROseq)))
PROseq$quartile <- findInterval(seq_along(PROseq),
                                ps_breaks,
                                rightmost.closed = TRUE)

test_that("Merging single dataset produces warning", {
    expect_warning(mergeGRangesData(PROseq))
})

test_that("Merging is non-destructive", {
    expect_identical(PROseq, suppressWarnings(mergeGRangesData(PROseq)))
    expect_identical(PROseq, mergeGRangesData(subset(PROseq, quartile == 1),
                                              subset(PROseq, quartile == 2),
                                              subset(PROseq, quartile == 3),
                                              subset(PROseq, quartile == 4)))
})

ps_10ranges <- PROseq[seq(1, 100, 10)] # ensure disjoint when resized below

test_that("Merging non-single-width GRanges produces warning", {
    expect_warning(mergeGRangesData(PROseq, resize(ps_ranges, 2)))
})

test_that("Merging sums signals in overlaps", {
    merge_10ranges <- mergeGRangesData(PROseq, ps_10ranges)
    expect_equal(length(PROseq), length(merge_10ranges))

    overlap_scores_combined <- score(merge_10ranges)[seq(1, 100, 10)]
    overlap_scores_original <- score(PROseq)[seq(1, 100, 10)]

    expect_equivalent(overlap_scores_combined, 2*overlap_scores_original)
    expect_equal(sum(score(merge_10ranges)),
                 sum(score(PROseq), overlap_scores_original))
})

test_that("Merging concatenates non-overlapping ranges", {
    merge_nonoverlap <- mergeGRangesData(ps_10ranges, shift(ps_10ranges, 1))
    expect_equal(length(merge_nonoverlap), 2*length(ps_10ranges))
    expect_true(all( merge_nonoverlap %in%
                         c(ps_10ranges, shift(ps_10ranges, 1)) ))
})



