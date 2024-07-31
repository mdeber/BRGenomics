context("Counting signal by region")
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

# for testing collapsed ranges (i.e. "stranded coverage" or similar)
pscov <- getStrandedCoverage(PROseq, ncores = 1)

pscov_list <- list(A_rep1 = pscov[seq(1, length(pscov), 6)],
                   A_rep2 = pscov[seq(4, length(pscov), 6)],
                   B_rep1 = pscov[seq(2, length(pscov), 6)],
                   B_rep2 = pscov[seq(5, length(pscov), 6)],
                   C_rep1 = pscov[seq(3, length(pscov), 6)],
                   C_rep2 = pscov[seq(6, length(pscov), 6)])

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
    countslnf <- getCountsByRegions(ps_list, txs_dm6_chr4, NF=1:6, ncores=1)
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



# CBR with Range Expansion ------------------------------------------------

context("signal counting with range expansion")

test_that("expanding ranges gives correct counts", {
    test_counts_exp <- getCountsByRegions(pscov, txs_dm6_chr4,
                                          expand_ranges = TRUE)
    expect_identical(test_counts, test_counts_exp)
})

test_that("failure to expand gives wrong counts", {
    expect_equivalent(test_counts[2], 59)
    expect_equivalent(getCountsByRegions(pscov, txs_dm6_chr4[2]), 54)
})

test_that("melt option works for single sample", {
    mcounts <- getCountsByRegions(PROseq, txs_dm6_chr4, melt = TRUE, ncores = 1)
    mcountsexp <- getCountsByRegions(pscov, txs_dm6_chr4, expand_ranges = TRUE,
                                     melt = TRUE, ncores = 1)
    expect_identical(mcountsexp, mcounts)
})

pscov_rename <- pscov
names(mcols(pscov_rename)) <- "signal"

test_that("expand_ranges with alternative metadata fields", {
    expect_message(getCountsByRegions(pscov_rename, txs_dm6_chr4,
                                      expand_ranges = TRUE))
    expect_equivalent(test_counts, getCountsByRegions(pscov_rename,
                                                      txs_dm6_chr4,
                                                      expand_ranges = TRUE,
                                                      field = "signal"))
})

pscov_rename$posnum <- seq_along(pscov_rename)

test_that("expand ranges counts over multiple metadata fields", {
    counts_multiple <- getCountsByRegions(pscov_rename, txs_dm6_chr4,
                                          field = c("signal", "posnum"),
                                          expand_ranges = TRUE,
                                          ncores = 1)
    expect_is(counts_multiple, "data.frame")
    expect_equivalent(counts_multiple$signal, test_counts)
})


countsl_exp <- getCountsByRegions(pscov_list, txs_dm6_chr4,
                                  expand_ranges = TRUE, ncores = 1)

test_that("expanded ranges counting over a list", {
    expect_is(countsl_exp, "data.frame")
    expect_equal(names(pscov_list), names(countsl_exp))
    expect_equal(nrow(countsl_exp), length(txs_dm6_chr4))
    expect_equal(sum(sapply(countsl_exp, sum)), sum(test_counts))

    # with nfs
    countsl_exp_nf <- getCountsByRegions(pscov_list, txs_dm6_chr4, NF = 1:6,
                                         expand_ranges = TRUE, ncores = 1)
    expect_equivalent(countsl_exp_nf, as.data.frame(Map("*", countsl_exp, 1:6)))
})

test_that("melt option works for list from expanded ranges", {
    mcountsl_exp <- getCountsByRegions(pscov_list, txs_dm6_chr4, melt = TRUE,
                                       expand_ranges = TRUE, ncores = 1)
    expect_is(mcountsl_exp, "data.frame")
    expect_equivalent(names(mcountsl_exp), c("region", "signal", "sample"))
    expect_equivalent(subset(mcountsl_exp, sample == "B_rep1")$signal,
                      countsl_exp[, "B_rep1"])
})

test_that("blacklisting works with expanded ranges", {
    blcounts <- getCountsByRegions(pscov, txs_dm6_chr4, expand_ranges = TRUE,
                                   blacklist = bl)

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

# melt:
countsdf <- getCountsByPositions(PROseq, txs_pr, melt = TRUE)
# melt + multiple fields:
fieldcountsdf <- getCountsByPositions(ps_rename, txs_pr,
                                      field = c("signal", "posnum"),
                                      melt = TRUE, ncores = 1)
# melt + multi-width regions:
meltlist <- getCountsByPositions(PROseq, txs_dm6_chr4,
                                 simplify.multi.widths = "list",
                                 melt = TRUE)

test_that("melting option works", {
    expect_equivalent(as.vector(t(countsmat)), countsdf$signal)
    expect_equal(ncol(countsdf), 3)

    # for multiple fields
    expect_is(fieldcountsdf, "data.frame")
    expect_equal(ncol(fieldcountsdf), 4) # has sample names now
    expect_equivalent(unique(fieldcountsdf[,4]), c("signal", "posnum"))

    # for multi-width regions
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



# CBP with Range Expansion ------------------------------------------------

context("Counting signal at positions with expanded ranges")

test_that("error if expand_ranges not set with collapsed input", {
    expect_error(getCountsByPositions(pscov, txs_pr))
})

countsmat_exp <- getCountsByPositions(pscov, txs_pr, expand_ranges = TRUE)

test_that("counts matrix correct with range expansion", {
    expect_identical(countsmat, countsmat_exp)
})

test_that("melt option for countsmatrix with range expansion", {
    countsmelt <- getCountsByPositions(pscov, txs_pr, expand_ranges = TRUE,
                                       melt = TRUE)
    expect_is(countsmelt, "data.frame")
    expect_equal(ncol(countsmelt), 3)
    expect_equal(length(countsmat), nrow(countsmelt))
    expect_equivalent(as.vector(t(countsmat)), countsmelt[,3])
})

test_that("counts matrix with range expansion & multiple metadata fields", {
    fieldcounts <- getCountsByPositions(pscov_rename, txs_pr,
                                        field = c("signal", "posnum"),
                                        expand_ranges = TRUE, ncores = 1)
    expect_is(fieldcounts, "list")
    expect_equivalent(names(fieldcounts), c("signal", "posnum"))
    expect_is(fieldcounts[[1]], "matrix")
    expect_is(fieldcounts[[2]], "matrix")
    expect_equivalent(fieldcounts[[1]], countsmat)
    # with normalization factor
    expect_equivalent(fieldcounts[[2]] * 0.5,
                      getCountsByPositions(pscov_rename, txs_pr,
                                           field = c("signal", "posnum"),
                                           expand_ranges = TRUE, NF = c(1, 0.5),
                                           ncores = 1)[[2]])
})

# cutting some tests to avoid bioc build timeouts

# test_that("counts matrix with list input & range expansion", {
#     countsl <- getCountsByPositions(pscov_list, txs_pr, expand_ranges = TRUE,
#                                     ncores = 1)
#     expect_is(countsl, "list")
#     expect_equivalent(names(countsl), names(pscov_list))
#     expect_is(countsl[[1]], "matrix")
#     expect_equivalent(dim(countsmat), dim(countsl[[1]]))
# 
#     arr3d <- simplify2array(countsl)
#     expect_equivalent(apply(arr3d, 1:2, sum), countsmat)
# })

# test_that("with range expansion, melt list input to single dataframe", {
#     countslmelt <- getCountsByPositions(pscov_list, txs_pr,
#                                         expand_ranges = TRUE, melt = TRUE,
#                                         ncores = 1)
#     expect_is(countslmelt, "data.frame")
#     expect_equal(ncol(countslmelt), 4)
#     expect_equivalent(unique(countslmelt[,4]), names(pscov_list))
# })

test_that("expanded ranges, error if multi-width is not explicit", {
    expect_error(getCountsByPositions(pscov, txs_dm6_chr4,
                                      expand_ranges = TRUE))
})

# countslist_exp <- getCountsByPositions(pscov, txs_dm6_chr4,
#                                        expand_ranges = TRUE,
#                                        simplify.multi.widths = "list")
# 
# test_that("can get list from multi-width counts matrix", {
#     expect_identical(countslist, countslist_exp)
# 
#     # with normalization factor
#     expect_equivalent(countslist_exp[[1]] * 0.5,
#                       getCountsByPositions(pscov, txs_dm6_chr4,
#                                            expand_ranges = TRUE,
#                                            simplify.multi.widths = "list",
#                                            NF = 0.5)[[1]])
# })

test_that("expand ranges, get 0-padded counts matrix for multi-width", {
    padmat <- getCountsByPositions(pscov, txs_dm6_chr4, expand_ranges = TRUE,
                                   simplify.multi.widths = "pad 0")
    expect_is(padmat, "matrix")
    expect_equivalent(dim(padmat),
                      c(length(txs_dm6_chr4), max(width(txs_dm6_chr4))))
    expect_equivalent(sapply(countslist, sum), rowSums(padmat))
})

# test_that("expand ranges, NA-padded counts matrix for multi-width", {
#     padmat <- getCountsByPositions(pscov, txs_dm6_chr4, expand_ranges = TRUE,
#                                    simplify.multi.widths = "pad NA")
#     expect_is(padmat, "matrix")
#     expect_equal(sum(is.na(rowSums(padmat))), length(txs_dm6_chr4) - 1)
# })

test_that("expand_ranges, melting option works", {
    expect_identical(countsdf, getCountsByPositions(pscov, txs_pr,
                                                    expand_ranges = TRUE,
                                                    melt = TRUE))
    # for multiple fields
    fieldcountsdfexp <- getCountsByPositions(pscov_rename, txs_pr,
                                             field = c("signal", "posnum"),
                                             expand_ranges = TRUE, melt = TRUE,
                                             ncores = 1)
    expect_is(fieldcountsdfexp, "data.frame")
    expect_equal(ncol(fieldcountsdf), 4) # has sample names now
    expect_equivalent(unique(fieldcountsdfexp[,4]), c("signal", "posnum"))
    expect_identical(subset(fieldcountsdf, sample == "signal")$signal,
                     subset(fieldcountsdfexp, sample == "signal")$signal)

    # for multi-width regions
    expect_identical(meltlist,
                     getCountsByPositions(pscov, txs_dm6_chr4,
                                     simplify.multi.widths = "list",
                                     expand_ranges = TRUE, melt = TRUE))
})

test_that("expand ranges, error on incorrect simplify argument", {
    expect_error(getCountsByPositions(pscov, txs_dm6_chr4, expand_ranges = TRUE,
                                      simplify.multi.widths = "notright"))

    # can't set argument if not multiwidth (not my ideal...)
    expect_error(getCountsByPositions(pscov, txs_pr, epand_ranges = TRUE,
                                      simplify.multi.widths = "list"))
})

test_that("expand ranges, error on empty regions", {
    expect_error(getCountsByPositions(pscov, txs_pr[0], expand_ranges = TRUE,
                                      ncores = 1))
})

test_that("expand ranges, arbitrary binning operations", {
    binsums <- getCountsByPositions(pscov, txs_pr, binsize = 10, FUN = sum,
                                    expand_ranges = TRUE)
    expect_is(binsums, "matrix")
    expect_equivalent(rowSums(binsums), rowSums(countsmat))
    expect_equivalent(dim(binsums), c(length(txs_dm6_chr4), 10))

    binmeans <- getCountsByPositions(pscov, txs_pr, binsize = 10, FUN = mean,
                                     expand_ranges = TRUE)
    expect_is(binmeans, "matrix")
    expect_equivalent(which(binsums == 0), which(binmeans == 0))
    expect_true(sum(binsums) != sum(binmeans))

    binmedians <- getCountsByPositions(pscov, txs_pr, binsize = 10,
                                       FUN = function(x) quantile(x, 0.5),
                                       expand_ranges = TRUE)
    expect_is(binmedians, "matrix")
    expect_equivalent(dim(binmedians), dim(binsums))
    expect_false(all(rowSums(binmedians) == 0))
    expect_false(all((binmedians == 0) == (binmeans == 0)))
})

test_that("expand ranges with blacklisting", {
    blcbp <- getCountsByPositions(pscov, txs_pr, expand_ranges = TRUE,
                                  blacklist = bl)

    test_bl <- countsmat
    test_bl[2,] <- 0
    expect_equivalent(blcbp, test_bl)

    blcbpna <- getCountsByPositions(pscov, txs_pr, expand_ranges = TRUE,
                                    blacklist = bl, NA_blacklisted = TRUE)
    expect_true(all(is.na(blcbpna[2,])))
})


# test_that("expand ranges + multiwidth, blacklisting over list works", {
#     blcl <- getCountsByPositions(pscov, txs_dm6_chr4, expand_ranges = TRUE,
#                                  blacklist = bl, simplify.multi.widths = "list")
#     expect_is(blcl, "list")
#     expect_true(all(blcl[[2]] == 0))
# 
#     blclna <- getCountsByPositions(pscov, txs_dm6_chr4, expand_ranges = TRUE,
#                                    blacklist = bl, NA_blacklisted = TRUE,
#                                    simplify.multi.widths = "list")
#     expect_true(all(is.na(blclna[[2]])))
#     expect_equivalent(blcl[[1]], blclna[[1]])
#     expect_equivalent(blcl[[3]], blclna[[3]])
# })



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



# Test Expand Ranges (Work in Progress) -----------------------------------

context("Signal counting by region with range expansion (coverage signal)")

# Get PRO-seq coverage (collapse adjacent positions);
#  subset for collapsed reads (ranges of interest for these comparisons)
cvg <- subset(getStrandedCoverage(PROseq, ncores = 1), width > 1)
cvg <- subset(cvg, width > 1)

# Get overlapping signal from the expanded data (basepair-resolution GRanges)
bpr <- subsetByOverlaps(PROseq, cvg)

# expanded- and unexpanded-counts on compressed (cvg) and expanded (bpr) data
bpr_noexp <- getCountsByRegions(bpr, txs_dm6_chr4)
bpr_exp <- getCountsByRegions(bpr, txs_dm6_chr4, expand_ranges = TRUE)
cvg_noexp <- getCountsByRegions(cvg, txs_dm6_chr4)
cvg_exp <- getCountsByRegions(cvg, txs_dm6_chr4, expand_ranges = TRUE)


test_that("expansion doesn't affect basepair-resolution GRanges", {
    # (b/c disjoint, single-width)
    expect_equivalent(bpr_exp, bpr_noexp)
})

test_that("expanding coverage data increases signal", {
    expect_true(all(cvg_noexp <= bpr_noexp))
    expect_true(all(cvg_exp >= cvg_noexp))
})

context("Signal counting by position with range expansion (coverage signal)")

# expanded- and unexpanded-counts on compressed (cvg) and expanded (bpr) data[*]
mat_bpr_noexp <- getCountsByPositions(bpr, promoters(txs_dm6_chr4, 0, 100))
mat_bpr_exp <- getCountsByPositions(bpr, promoters(txs_dm6_chr4, 0, 100),
                                    expand_ranges = TRUE)
# (this currently just here to check if there's an error)
mat_cvg_exp <- getCountsByPositions(cvg, promoters(txs_dm6_chr4, 0, 100),
                                    expand_ranges = TRUE)

# [*] but can't get unexpanded counts on expanded data
test_that("error on unexpanded counts on expanded data", {
    expect_error(getCountsByPositions(cvg, promoters(txs_dm6_chr4, 0, 100)))
})

# (nothing affected - must write test cases)
test_that("expansion doesn't affect countsmat when no multi-width ranges", {
    expect_equivalent(mat_bpr_exp, mat_bpr_noexp)
})

