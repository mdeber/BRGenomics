# Merging GRanges ---------------------------------------------------------

context("Merging GRanges data")
library(BRGenomics)
data("PROseq")

ps_breaks <- seq(1, length(PROseq), floor(0.25*length(PROseq)))
PROseq$quartile <- findInterval(seq_along(PROseq),
                                ps_breaks,
                                rightmost.closed = TRUE)

test_that("Merging is non-destructive", {
    expect_identical(subset(PROseq, select = score),
                     suppressWarnings(mergeGRangesData(PROseq, ncores = 1)))
    expect_identical(subset(PROseq, select = score),
                     mergeGRangesData(subset(PROseq, quartile == 1),
                                      subset(PROseq, quartile == 2),
                                      subset(PROseq, quartile == 3),
                                      subset(PROseq, quartile == 4),
                                      ncores = 1))
})

ps_10ranges <- PROseq[seq(1, 100, 10)] # ensure disjoint when resized below

test_that("Merging non-single-width GRanges produces warning", {
    expect_warning(mergeGRangesData(PROseq, resize(ps_10ranges, 2), ncores = 1))
})

test_that("Merging sums signals in overlaps", {
    merge_10ranges <- mergeGRangesData(PROseq, ps_10ranges, ncores = 1)
    expect_equal(length(PROseq), length(merge_10ranges))

    overlap_scores_combined <- score(merge_10ranges)[seq(1, 100, 10)]
    overlap_scores_original <- score(PROseq)[seq(1, 100, 10)]

    expect_equivalent(overlap_scores_combined, 2*overlap_scores_original)
    expect_equal(sum(score(merge_10ranges)),
                 sum(score(PROseq), overlap_scores_original))
})

test_that("Merging concatenates non-overlapping ranges", {
    merge_nonoverlap <- mergeGRangesData(ps_10ranges, shift(ps_10ranges, 1),
                                         ncores = 1)
    expect_equal(length(merge_nonoverlap), 2*length(ps_10ranges))
    # %in% operator fails depending on environment?
    # expect_true(all( merge_nonoverlap %in%
    #                      c(ps_10ranges, shift(ps_10ranges, 1)) ))
})

gr1 <- PROseq[10:13]
gr2 <- PROseq[12:15]
gr_multi <- mergeGRangesData(gr1, gr2, multiplex = TRUE, ncores = 1)

test_that("Multiplexed merging works with input GRanges", {
    expect_equivalent(names(mcols(gr_multi)), c("gr1", "gr2"))
    expect_equivalent(sapply(mcols(gr_multi), sum),
                      c(sum(gr1$score), sum(gr2$score)))
    expect_equivalent(subset(gr_multi, gr1 > 0, select = gr1)$gr1,
                      gr1$score)
})

test_that("Multiplexed merging works with listed input", {
    expect_equivalent(gr_multi, mergeGRangesData(list(gr1 = gr1, gr2 = gr2),
                                                 multiplex = TRUE, ncores = 1))
    expect_equivalent(gr_multi, mergeGRangesData(list(gr1 = gr1),
                                                 list(gr2 = gr2),
                                                 multiplex = TRUE, ncores = 1))
})

# Sampling GRanges --------------------------------------------------------

context("Random subsampling GRanges data")

set.seed(11)
ps_tenth <- subsampleGRanges(PROseq, prop = 0.1)

test_that("Subsampling produces GRanges", {
    expect_is(ps_tenth, "GRanges")
})

test_that("Subsampling gives correct signal count", {
    expect_equal(round(0.1*sum(score(PROseq))), sum(score(ps_tenth)))
})

test_that("Subsamplied reproduced using set.seed", {
    expect_equal(length(ps_tenth), 6330)
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

test_that("can subsample multiplexed GRanges", {
    ssmulti <- subsampleGRanges(gr_multi, prop = 0.1, field = c("gr1", "gr2"),
                                ncores = 1)
    expect_is(ssmulti, "GRanges")
    expect_equivalent(names(mcols(ssmulti)), c("gr1", "gr2"))
    expect_equal(sum(ssmulti$gr1), round(0.1*sum(gr_multi$gr1)))
    expect_equal(sum(ssmulti$gr2), round(0.1*sum(gr_multi$gr2)))
})

test_that("can subsample when simple normalization factor was applied", {
    norm_ps <- PROseq
    score(norm_ps) <- 0.89 * score(norm_ps)
    expect_warning(subsampleGRanges(norm_ps, prop = 0.1))

    set.seed(11)
    norm_tenth <- suppressWarnings(subsampleGRanges(norm_ps, prop = 0.1))
    signal_ps_tenth <- sum(score(ps_tenth))
    signal_norm_tenth <- sum(score(norm_tenth))
    expect_true(abs(0.89*signal_ps_tenth - signal_norm_tenth) < 1)
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
