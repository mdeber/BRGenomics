context("counting signal by region")
library(BRGenomics)

data("PROseq")
data("txs_dm6_chr4")

test_counts <- getCountsByRegions(PROseq, txs_dm6_chr4)

test_that("signal counts correct returns correct vector", {
    expect_is(test_counts, "integer")
    expect_equal(length(test_counts), length(txs_dm6_chr4))
    expect_equivalent(test_counts[1:4], c(1, 58, 13, 125))
})

ps_rename <- PROseq
names(mcols(ps_rename)) <- "signal"

test_that("signal counts works with alternative metadata fields", {
    expect_message(getCountsByRegions(ps_rename, txs_dm6_chr4))
    expect_equivalent(test_counts, getCountsByRegions(ps_rename, txs_dm6_chr4,
                                                      field = "signal"))
})

ps_rename$posnum <- seq_along(ps_rename)

test_that("signal counts work over multiple metadata fields", {
    counts_multiple <- getCountsByRegions(ps_rename, txs_dm6_chr4,
                                          field = c("signal", "posnum"),
                                          ncores = 2)
    expect_is(counts_multiple, "data.frame")
    expect_equivalent(counts_multiple$signal, test_counts)
})

# getCountsByPositions ----------------------------------------------------

context("Counting signal at positions within regions")

txs_pr <- promoters(txs_dm6_chr4, 0, 100)
countsmat <- getCountsByPositions(PROseq, txs_pr)

test_that("simple counts matrix correctly calculated", {
    expect_is(countsmat, "matrix")
    expect_equivalent(dim(countsmat), c(length(txs_dm6_chr4), 100))
    expect_equivalent(rowSums(countsmat), getCountsByRegions(PROseq, txs_pr))
    # with normalization factor
    expect_equivalent(getCountsByPositions(PROseq, txs_pr, NF = 0.5),
                      countsmat * 0.5)
})

test_that("simple counts matrix works over multiple metadata fields", {
    fieldcounts <- getCountsByPositions(ps_rename, txs_pr,
                                        field = c("signal", "posnum"),
                                        ncores = 2)
    expect_is(fieldcounts, "list")
    expect_equivalent(names(fieldcounts), c("signal", "posnum"))
    expect_is(fieldcounts[[1]], "matrix")
    expect_is(fieldcounts[[2]], "matrix")
    expect_equivalent(fieldcounts[[1]], countsmat)
    # with normalization factor
    expect_equivalent(fieldcounts[[2]] * 0.5,
                      getCountsByPositions(ps_rename, txs_pr,
                                           field = c("signal", "posnum"),
                                           NF = c(1, 0.5), ncores = 2)[[2]])
})

test_that("error if multi-width is not explicit", {
    expect_error(getCountsByPositions(PROseq, txs_dm6_chr4))
})

countslist <- getCountsByPositions(PROseq, txs_dm6_chr4,
                                   simplify.multi.widths = "list")

test_that("can get list from multi-width counts matrix", {
    expect_is(countslist, "list")
    expect_equal(length(countslist), nrow(countsmat))
    expect_equivalent(width(txs_dm6_chr4), lengths(countslist))

    idx <- which(width(txs_dm6_chr4) >= 100)
    expect_equivalent(countsmat[idx, ],
                      t(sapply(countslist[idx], "[", 1:100)))

    # with normalization factor
    expect_equivalent(countslist[[1]] * 0.5,
                      getCountsByPositions(PROseq, txs_dm6_chr4,
                                           simplify.multi.widths = "list",
                                           NF = 0.5)[[1]])
})

test_that("can get 0-padded counts matrix for multi-width regions", {
    padmat <- getCountsByPositions(PROseq, txs_dm6_chr4,
                                   simplify.multi.widths = "pad 0")
    expect_is(padmat, "matrix")
    expect_equivalent(dim(padmat),
                      c(length(txs_dm6_chr4), max(width(txs_dm6_chr4))))
    expect_equivalent(sapply(countslist, sum), rowSums(padmat))
})

test_that("can get NA-padded counts matrix for multi-width regions", {
    padmat <- getCountsByPositions(PROseq, txs_dm6_chr4,
                                   simplify.multi.widths = "pad NA")
    expect_is(padmat, "matrix")
    expect_equal(sum(is.na(rowSums(padmat))), length(txs_dm6_chr4) - 1)
})

test_that("melting option works", {
    countsdf <- getCountsByPositions(PROseq, txs_pr, melt = TRUE)
    expect_equivalent(as.vector(t(countsmat)), countsdf$signal)
    expect_equal(ncol(countsdf), 3)

    # for multiple fields
    fieldcountsdf <- getCountsByPositions(ps_rename, txs_pr,
                                          field = c("signal", "posnum"),
                                          melt = TRUE, ncores = 2)
    expect_is(fieldcountsdf, "data.frame")
    expect_equal(ncol(fieldcountsdf), 4) # has sample names now
    expect_equivalent(unique(fieldcountsdf[,4]), c("signal", "posnum"))

    # for multi-width regions
    meltlist <- getCountsByPositions(PROseq, txs_dm6_chr4,
                                     simplify.multi.widths = "list",
                                     melt = TRUE)
    expect_is(meltlist, "data.frame")
    expect_equal(ncol(meltlist), 3)
    expect_equivalent(unlist(countslist), meltlist[, 3])
})

test_that("error on incorrect simplify argument", {
    expect_error(getCountsByPositions(PROseq, txs_dm6_chr4,
                                      simplify.multi.widths = "notright"))
})

test_that("can perform arbitrary binning operations on count matrix", {
    binsums <- getCountsByPositions(PROseq, txs_pr, binsize = 10, FUN = sum)
    expect_is(binsums, "matrix")
    expect_equivalent(rowSums(binsums), rowSums(countsmat))
    expect_equivalent(dim(binsums), c(length(txs_dm6_chr4), 10))

    binmeans <- getCountsByPositions(PROseq, txs_pr, binsize = 10, FUN = mean)
    expect_is(binmeans, "matrix")
    expect_equivalent(which(binsums == 0), which(binmeans == 0))
    expect_true(sum(binsums) != sum(binmeans))

    binmedians <- getCountsByPositions(PROseq, txs_pr, binsize = 10,
                                       FUN = function(x) quantile(x, 0.5))
    expect_is(binmedians, "matrix")
    expect_equivalent(dim(binmedians), dim(binsums))
    expect_false(all(rowSums(binmedians) == 0))
    expect_false(all((binmedians == 0) == (binmeans == 0)))
})

# Calculating pause indices -----------------------------------------------

context("Calculating pausing indices")

txs_gb <- flank(txs_pr, 500, start = FALSE)
pidx <- getPausingIndices(PROseq, txs_pr, txs_gb)

test_that("can calculate pausing indices", {
    expect_is(pidx, "numeric")
    expect_equivalent(round(pidx[1:4]), c(0, 2, 3, 0))
})

test_that("error if regions.pr not matching regions.gb", {
    expect_error(getPausingIndices(PROseq, txs_pr, txs_gb[1:10]))
})

test_that("length.normalize works (check when FALSE)", {
    pidx_norm <- getPausingIndices(PROseq, txs_pr, txs_gb,
                                   length.normalize = FALSE)
    expect_true(!all(pidx == pidx_norm))
    expect_equivalent(round(pidx_norm)[1:4], c(0, 0, 1, 0))
})

pidx_re <- getPausingIndices(PROseq, txs_pr, txs_gb, remove.empty = TRUE)

test_that("remove.empty works", {
    expect_true(length(pidx_re) < length(pidx))
    expect_equal(length(pidx_re),
                 length( which(getCountsByRegions(PROseq, txs_pr) != 0) ))
})

test_that("can get pause indices over multiple fields", {
    pidx_multi <- getPausingIndices(ps_rename, txs_pr, txs_gb,
                                    field = c("signal", "posnum"),
                                    ncores = 2)
    expect_is(pidx_multi, "data.frame")
    expect_equivalent(names(pidx_multi), c("signal", "posnum"))
    expect_equivalent(pidx_multi[, 1], pidx)

    pidx_multi_re <- getPausingIndices(ps_rename, txs_pr, txs_gb,
                                       field = c("signal", "posnum"),
                                       remove.empty = TRUE, ncores = 2)
    expect_is(pidx_multi, "data.frame")
    expect_equivalent(names(pidx_multi), c("signal", "posnum"))
    expect_equivalent(pidx_multi_re[, 1], pidx_re)
})
