context("Random subsampling GRanges data")
library(BRGenomics)
data("PROseq")

set.seed(11)
ps_tenth <- subsampleGRanges(PROseq, prop = 0.1)

test_that("Subsampling produces GRanges", {
    expect_is(ps_tenth, "GRanges")
})

test_that("Subsampling gives correct signal count", {
    expect_equal(floor(0.1*sum(score(PROseq))), sum(score(ps_tenth)))
})

test_that("Subsamplied reproduced using set.seed", {
    expect_equal(length(ps_tenth), 6397)
})

test_that("can subsample with field = NULL", {
    reads <- rep(PROseq, times = PROseq$score)
    score(reads) <- NULL
    reads_tenth <- subsampleGRanges(reads, prop = 0.1, field = NULL)
    expect_is(reads_tenth, "GRanges")
    expect_equivalent(mcols(reads), mcols(reads_tenth))
    expect_equal(length(reads_tenth), floor(0.1*length(reads)))
    expect_error(subsampleGRanges(reads, prop = 0.1))
})

test_that("can subsample when simple normalization factor was applied", {
    norm_ps <- PROseq
    score(norm_ps) <- 0.89 * score(norm_ps)
    expect_warning(subsampleGRanges(norm_ps, prop = 0.1))

    set.seed(11)
    norm_tenth <- suppressWarnings(subsampleGRanges(norm_ps, prop = 0.1))

    expect_equal(length(norm_tenth), length(ps_tenth))
    expect_equivalent(0.89*ps_tenth$score, norm_tenth$score)
})

test_that("error when incorrect input arguments", {
    expect_error(subsampleGRanges(PROseq))
    expect_error(subsampleGRanges(PROseq, prop = 0.1, n = 7415))
})

test_that("error when non-simple normalization", {
    breakpoint <- floor(0.5*length(PROseq))
    norm_1 <- PROseq[1:breakpoint]
    norm_2 <- PROseq[(breakpoint+1):length(PROseq)]
    score(norm_1) <- 0.75*norm_1$score
    score(norm_2) <- 0.50*norm_2$score
    norm_mix <- c(norm_1, norm_2)

    expect_error(subsampleGRanges(norm_mix, prop = 0.1))
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
    expect_warning(mergeGRangesData(PROseq, resize(ps_10ranges, 2)))
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
    # %in% operator fails depending on environment?
    # expect_true(all( merge_nonoverlap %in%
    #                      c(ps_10ranges, shift(ps_10ranges, 1)) ))
})



