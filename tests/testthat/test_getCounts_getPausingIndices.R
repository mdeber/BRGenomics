context("counting signal by region")
library(BRGenomics)

data("PROseq")
data("txs_dm6_chr4")

# for testing list inputs
ps_list <- list(A_rep1 = PROseq[seq(1, length(PROseq), 6)],
                A_rep2 = PROseq[seq(4, length(PROseq), 6)],
                B_rep1 = PROseq[seq(2, length(PROseq), 6)],
                B_rep2 = PROseq[seq(5, length(PROseq), 6)],
                C_rep1 = PROseq[seq(3, length(PROseq), 6)],
                C_rep2 = PROseq[seq(6, length(PROseq), 6)])


# getCountsByRegions ------------------------------------------------------

test_counts <- getCountsByRegions(PROseq, txs_dm6_chr4)

test_that("signal counts correct returns correct vector", {
    expect_is(test_counts, "integer")
    expect_equal(length(test_counts), length(txs_dm6_chr4))
    expect_equivalent(test_counts[1:4], c(1, 59, 13, 126))
})

test_that("melt option works for single sample", {
    mcounts <- getCountsByRegions(PROseq, txs_dm6_chr4, melt = TRUE, ncores = 1)
    expect_is(mcounts, "data.frame")
    expect_equivalent(names(mcounts), c("region", "signal"))
    expect_equivalent(mcounts$signal, test_counts)
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
                                          ncores = 1)
    expect_is(counts_multiple, "data.frame")
    expect_equivalent(counts_multiple$signal, test_counts)
})

countsl <- getCountsByRegions(ps_list, txs_dm6_chr4, ncores = 1)

test_that("signal counting works over a list", {
    expect_is(countsl, "data.frame")
    expect_equal(names(ps_list), names(countsl))
    expect_equal(nrow(countsl), length(txs_dm6_chr4))
    expect_equal(sum(sapply(countsl, sum)), sum(test_counts))

    # with nfs
    countslnf <- getCountsByRegions(ps_list, txs_dm6_chr4, NF=1:6, ncores=2)
    expect_equivalent(countslnf, as.data.frame(Map("*", countsl, 1:6)))
})

test_that("melt option works for list", {
    mcountsl <- getCountsByRegions(ps_list, txs_dm6_chr4, melt = TRUE,
                                   ncores = 1)
    expect_is(mcountsl, "data.frame")
    expect_equivalent(names(mcountsl), c("region", "signal", "sample"))
    expect_equivalent(subset(mcountsl, sample == "B_rep1")$signal,
                      countsl[, "B_rep1"])
})

test_that("error if NFs given don't match input", {
    expect_error(getCountsByRegions(ps_list, txs_dm6_chr4, NF = 1:5,
                                    ncores = 1))
    expect_error(getCountsByRegions(ps_rename, txs_dm6_chr4,
                                    field = c("signal", "posnum"),
                                    NF = 1, ncores = 1))
})

bl <- txs_dm6_chr4[2]

test_that("blacklisting works", {
    blcounts <- getCountsByRegions(PROseq, txs_dm6_chr4, blacklist = bl)

    test_bl <- test_counts
    test_bl[2] <- 0
    expect_equivalent(blcounts, test_bl)
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

test_that("melt option melts simple counts matrix", {
    countsmelt <- getCountsByPositions(PROseq, txs_pr, melt = TRUE)
    expect_is(countsmelt, "data.frame")
    expect_equal(ncol(countsmelt), 3)
    expect_equal(length(countsmat), nrow(countsmelt))
    expect_equivalent(as.vector(t(countsmat)), countsmelt[,3])
})

test_that("simple counts matrix works over multiple metadata fields", {
    fieldcounts <- getCountsByPositions(ps_rename, txs_pr,
                                        field = c("signal", "posnum"),
                                        ncores = 1)
    expect_is(fieldcounts, "list")
    expect_equivalent(names(fieldcounts), c("signal", "posnum"))
    expect_is(fieldcounts[[1]], "matrix")
    expect_is(fieldcounts[[2]], "matrix")
    expect_equivalent(fieldcounts[[1]], countsmat)
    # with normalization factor
    expect_equivalent(fieldcounts[[2]] * 0.5,
                      getCountsByPositions(ps_rename, txs_pr,
                                           field = c("signal", "posnum"),
                                           NF = c(1, 0.5), ncores = 1)[[2]])
})

test_that("simple counts matrix works for list input", {
    countsl <- getCountsByPositions(ps_list, txs_pr, ncores = 1)
    expect_is(countsl, "list")
    expect_equivalent(names(countsl), names(ps_list))
    expect_is(countsl[[1]], "matrix")
    expect_equivalent(dim(countsmat), dim(countsl[[1]]))

    arr3d <- simplify2array(countsl)
    expect_equivalent(apply(arr3d, 1:2, sum), countsmat)
})

test_that("melting for list input returns single dataframe", {
    countslmelt <- getCountsByPositions(ps_list, txs_pr, melt=TRUE, ncores = 1)
    expect_is(countslmelt, "data.frame")
    expect_equal(ncol(countslmelt), 4)
    expect_equivalent(unique(countslmelt[,4]), names(ps_list))
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
                                          melt = TRUE, ncores = 1)
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

    # can't set argument if not multiwidth (not my ideal...)
    expect_error(getCountsByPositions(PROseq, txs_pr,
                                      simplify.multi.widths = "list"))
})

test_that("error on empty regions", {
    expect_error(getCountsByPositions(PROseq, txs_pr[0], ncores = 1))
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

test_that("blacklisting works", {
    blcbp <- getCountsByPositions(PROseq, txs_pr, blacklist = bl)

    test_bl <- countsmat
    test_bl[2,] <- 0
    expect_equivalent(blcbp, test_bl)

    blcbpna <- getCountsByPositions(PROseq, txs_pr, blacklist = bl,
                                    NA_blacklisted = TRUE)
    expect_true(all(is.na(blcbpna[2,])))
})


test_that("multiwidth, blacklisting over list works", {
    blcl <- getCountsByPositions(PROseq, txs_dm6_chr4, blacklist = bl,
                                 simplify.multi.widths = "list")
    expect_is(blcl, "list")
    expect_true(all(blcl[[2]] == 0))

    blclna <- getCountsByPositions(PROseq, txs_dm6_chr4, blacklist = bl,
                                   simplify.multi.widths = "list",
                                   NA_blacklisted = TRUE)
    expect_true(all(is.na(blclna[[2]])))
    expect_equivalent(blcl[[1]], blclna[[1]])
    expect_equivalent(blcl[[3]], blclna[[3]])
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

test_that("melt option works for single dataset", {
    mpidx <- getPausingIndices(PROseq, txs_pr, txs_gb, melt = TRUE)
    expect_is(mpidx, "data.frame")
    expect_equivalent(names(mpidx), c("region", "pauseIndex"))
    expect_equivalent(mpidx$pauseIndex, pidx)
})

test_that("can get pause indices over multiple fields", {
    pidx_multi <- getPausingIndices(ps_rename, txs_pr, txs_gb,
                                    field = c("signal", "posnum"),
                                    ncores = 1)
    expect_is(pidx_multi, "data.frame")
    expect_equivalent(names(pidx_multi), c("signal", "posnum"))
    expect_equivalent(pidx_multi[, 1], pidx)

    # with remove.empty
    pidx_multi_re <- getPausingIndices(ps_rename, txs_pr, txs_gb,
                                       field = c("signal", "posnum"),
                                       remove.empty = TRUE, ncores = 1)
    expect_is(pidx_multi, "data.frame")
    expect_equivalent(names(pidx_multi), c("signal", "posnum"))
    expect_equivalent(pidx_multi_re[, 1], pidx_re)
})

pidxl <- getPausingIndices(ps_list, txs_pr, txs_gb, ncores = 1)

test_that("can get pause indices for a list", {
    expect_is(pidxl, "data.frame")
    expect_equivalent(colnames(pidxl), names(ps_list))
    expect_equal(nrow(pidxl), length(txs_pr))

    # with remove.empty
    pidxl_re <- getPausingIndices(ps_list, txs_pr, txs_gb, remove.empty = TRUE,
                                  ncores = 1)
    expect_equal(nrow(pidxl_re), 189)
})

test_that("melt option works for a list", {
    mpidxl <- getPausingIndices(ps_list, txs_pr, txs_gb, melt = TRUE,
                                ncores = 1)
    expect_is(mpidxl, "data.frame")
    expect_equivalent(names(mpidxl), c("region", "pauseIndex", "sample"))
    expect_equivalent(subset(mpidxl, sample == "B_rep1")$pauseIndex,
                      pidxl[, "B_rep1"])
})

test_that("blacklisting works", {
    blpidx <- getPausingIndices(PROseq, txs_pr, txs_gb, blacklist = bl,
                                ncores = 1)
    expect_true(is.finite(pidx[2]))
    expect_true(!is.finite(blpidx[2]))

    # with a list
    blpidxl <- getPausingIndices(ps_list, txs_pr, txs_gb, blacklist = bl,
                                 ncores = 1)
    expect_true(identical(blpidxl[1,], pidxl[1,]))
    expect_true(!identical(blpidxl[2,], pidxl[2,]))
})
