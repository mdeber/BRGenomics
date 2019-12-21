context("genebodies function")
library(BRGenomics)

data("txs_dm6_chr4")

test_that("Function returns GRanges", {
    expect_is(genebodies(txs_dm6_chr4), "GRanges")
})

test_that("Length filtering correct for fixes = start, end", {
    gb_list <- list(genebodies(txs_dm6_chr4, 0, 0, min_window = 500),
                    genebodies(txs_dm6_chr4, 500, 0, min_window = 0),
                    genebodies(txs_dm6_chr4, 0, -500, min_window = 0),
                    genebodies(txs_dm6_chr4, 250, -250, min_window = 0),
                    genebodies(txs_dm6_chr4, 125, -125, min_window = 250))
    lengths <- vapply(gb_list, length, FUN.VALUE = integer(1))
    expect_equal(length(unique(lengths)), 1)
})

test_that("Region sizes correct for fixes = start, end", {
    gb <- genebodies(txs_dm6_chr4, 300, -300, min_window = 500)
    expect_equivalent(width(gb[c(1, 151, 301)]), c(3561, 39277, 6866))
})

test_that("Length filtering & region sizes correct for fixes = start, start", {
    gb_list <- list(genebodies(txs_dm6_chr4,
                               0, 500, fix_end = "start", min_window = 0),
                    genebodies(txs_dm6_chr4,
                               250, 500, fix_end = "start", min_window = 0),
                    genebodies(txs_dm6_chr4,
                               0, 250, fix_end = "start", min_window = 250),
                    genebodies(txs_dm6_chr4,
                               -250, 250, fix_end = "start", min_window = 250),
                    genebodies(txs_dm6_chr4,
                               -750, 250, fix_end = "start", min_window = 250),
                    genebodies(txs_dm6_chr4,
                               -750, -250, fix_end = "start", min_window = 750))

    lengths <- vapply(gb_list, length, FUN.VALUE = integer(1))
    expect_equal(length(unique(lengths)), 1)

    widths <- vapply(gb_list, function(x) length(unique(width(x))),
                     FUN.VALUE = integer(1))
    expect_true(all(widths == 1))
})

test_that("Length filtering & region sizes correct for fixes = end, end", {
    gb_list <- list(genebodies(txs_dm6_chr4,
                               0, 500, fix_start = "end", min_window = 500),
                    genebodies(txs_dm6_chr4,
                               250, 500, fix_start = "end", min_window = 750),
                    genebodies(txs_dm6_chr4,
                               -250, 250, fix_start = "end", min_window = 250),
                    genebodies(txs_dm6_chr4,
                               -500, 250, fix_start = "end", min_window = 0),
                    genebodies(txs_dm6_chr4,
                               -250, 750, fix_start = "end", min_window = 250))

    lengths <- vapply(gb_list, length, FUN.VALUE = integer(1))
    expect_equal(length(unique(lengths)), 1)

    widths <- vapply(gb_list, function(x) length(unique(width(x))),
                     FUN.VALUE = integer(1))
    expect_true(all(widths == 1))
})
