context("genebodies function")
library(BRGenomics)

data("txs_dm6_chr4")

test_that("Function returns GRanges", {
    expect_is(genebodies(txs_dm6_chr4), "GRanges")
})

test_that("Length filtering correct for fixes = start, end", {
    gb_list <- list(genebodies(txs_dm6_chr4, 0, 0, min.window = 500),
                    genebodies(txs_dm6_chr4, 500, 0, min.window = 0),
                    genebodies(txs_dm6_chr4, 0, -500, min.window = 0),
                    genebodies(txs_dm6_chr4, 250, -250, min.window = 0),
                    genebodies(txs_dm6_chr4, 125, -125, min.window = 250))
    lengths <- vapply(gb_list, length, FUN.VALUE = integer(1))
    expect_equal(length(unique(lengths)), 1)
})

test_that("Region sizes correct for fixes = start, end", {
    gb <- genebodies(txs_dm6_chr4, 300, -300, min.window = 500)
    expect_equivalent(width(gb[c(1, 151, 301)]), c(3561, 39277, 6866))
})

test_that("Length filtering & region sizes correct for fixes = start, start", {
    gb_list <- list(genebodies(txs_dm6_chr4,
                               0, 500, fix.end = "start", min.window = 500),
                    genebodies(txs_dm6_chr4,
                               0, 250, fix.end = "start", min.window = 500),
                    genebodies(txs_dm6_chr4,
                               -250, 250, fix.end = "start", min.window = 500))

    lengths <- vapply(gb_list, length, FUN.VALUE = integer(1))
    expect_equal(length(unique(lengths)), 1)

    widths <- vapply(gb_list, function(x) length(unique(width(x))),
                     FUN.VALUE = integer(1))
    expect_true(all(widths == 1))
})

test_that("Length filtering & region sizes correct for fixes = end, end", {
    gb_list <- list(genebodies(txs_dm6_chr4,
                               0, 500, fix.start = "end", min.window = 500),
                    genebodies(txs_dm6_chr4,
                               250, 500, fix.start = "end", min.window = 500),
                    genebodies(txs_dm6_chr4,
                               -250, 250, fix.start = "end", min.window = 500))

    lengths <- vapply(gb_list, length, FUN.VALUE = integer(1))
    expect_equal(length(unique(lengths)), 1)

    widths <- vapply(gb_list, function(x) length(unique(width(x))),
                     FUN.VALUE = integer(1))
    expect_true(all(widths == 1))
})

test_that("genebodies() identical to promoters() when expected", {
    expect_equivalent(genebodies(txs_dm6_chr4, -50, 150, fix.end = "start"),
                      promoters(txs_dm6_chr4, 50, 150))
})

test_that("errors on invalid inputs", {
    expect_error(genebodies(txs_dm6_chr4, 0, 500, fix.start = "typo"))
    expect_error(genebodies(txs_dm6_chr4, 0, 500,
                 fix.start = "end", fix.end = "start"))
    expect_error(genebodies(txs_dm6_chr4, 200, 100, fix.end = "start"))
})

test_that("unstranded ranges filtered, and warning given", {
    unst_range <- txs_dm6_chr4[1:10]
    strand(unst_range) <- "*"
    expect_warning(genebodies(c(txs_dm6_chr4, unst_range)))
    expect_equivalent(genebodies(txs_dm6_chr4),
                      suppressWarnings(genebodies(c(txs_dm6_chr4, unst_range))))
})





