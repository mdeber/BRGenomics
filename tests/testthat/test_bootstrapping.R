context("bootstrapping signal counts")
library(BRGenomics)

data("PROseq")
data("txs_dm6_chr4")

txs_pr <- promoters(txs_dm6_chr4, 0, 100)
countsmat <- getCountsByPositions(PROseq, txs_pr)

set.seed(11)
df <- metaSubsampleMatrix(countsmat, first.output.xval = -10, ncores = 1)

test_that("can get bootstrapped means from counts matrix", {
    expect_is(df, "data.frame")
    expect_equal(nrow(df), 100)
    expect_equivalent(names(df),
                      c("x", "mean", "lower", "upper", "sample.name"))
    expect_equivalent(df$x, seq_len(100) - 11)
    expect_equivalent(round(df$mean*100)[1:5], c(3, 0, 6, 6, 0))
})

test_that("sampling GRanges same as sampling matrix", {
    set.seed(11)
    expect_equivalent(df, metaSubsample(PROseq, txs_pr,
                                        first.output.xval = -10,
                                        sample.name = "countsmat", ncores = 1))
})

test_that("error produced for variable-width regions", {
    expect_error(metaSubsample(PROseq, txs_dm6_chr4, ncores = 1))
})

test_that("error produced for too few iterations", {
    expect_silent(metaSubsample(PROseq, txs_pr, n.iter = 5, ncores = 1))
    expect_error(metaSubsample(PROseq, txs_pr, n.iter = 4, ncores = 1))
})

test_that("message produced when n.iter = 1", {
    expect_message(metaSubsample(PROseq, txs_pr, n.iter = 1, ncores = 1))
})

set.seed(11)
df_bin <- metaSubsample(PROseq, txs_pr, binsize = 10,
                        first.output.xval = -10, ncores = 1)

test_that("binning correct with bootstrapping", {
    expect_equivalent(nrow(df_bin), 10)
    expect_equivalent(df_bin$x, seq(-5.5, 84.5, 10))
})

test_that("binning correct with bootstrapping counts matrix", {
    set.seed(11)
    df_binmat <- metaSubsampleMatrix(countsmat, binsize = 10,
                                     first.output.xval = -10,
                                     sample.name = "PROseq", ncores = 1)
    expect_equivalent(df_bin$x, df_binmat$x)
})

test_that("multicore bootstrapping successful", {
    set.seed(11)
    expect_equivalent(df, metaSubsample(PROseq, txs_pr,
                                        first.output.xval = -10,
                                        sample.name = "countsmat",
                                        ncores = 2))
})

ps_rename <- PROseq
names(mcols(ps_rename)) <- "signal"
ps_rename$posnum <- seq_along(ps_rename)

test_that("bootstrapping successful over several fields", {
    set.seed(11)
    df_alt <- metaSubsample(ps_rename, txs_pr, field = "signal",
                            first.output.xval = -10, sample.name = "countsmat",
                            ncores = 2)
    expect_equivalent(df, df_alt)

    set.seed(11)
    df_m <- metaSubsample(ps_rename, txs_pr, field = c("signal", "posnum"),
                          first.output.xval = -10,
                          sample.name = c("countsmat", "posnum"),
                          ncores = 2)
    expect_is(df_m, "data.frame")
    expect_equal(nrow(df_m), 2 * nrow(df))
    expect_equivalent(as.character(unique(df_m$sample.name)),
                      c("countsmat", "posnum"))
    expect_equivalent(subset(df_m, sample.name == "countsmat"), df)

    expect_warning(metaSubsample(ps_rename, txs_pr,
                                 field = c("signal", "posnum"),
                                 sample.name = "countsmat",
                                 remove.empty = TRUE, ncores = 1))
})



