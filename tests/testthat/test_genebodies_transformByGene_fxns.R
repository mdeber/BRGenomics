
# genebodies fxn ----------------------------------------------------------

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


# reduceByGene fxn --------------------------------------------------------

context("reduceByGene function")

rbg_gr <- GRanges(seqnames = "chr1",
                  ranges = IRanges(start = c(1, 3, 7, 10),
                                   end = c(4, 5, 9, 11)),
                  strand = "+")

# just change the gene names with gn; make a list of expected ranges like so:
# exp = list(a = 1:5, a = 7:9, b = 10:11)
testfun <- function(gn, disjoin = FALSE, exp) {
    out <- reduceByGene(rbg_gr, gn, disjoin)
    ex.r <- Reduce(c, (lapply(exp, function(x) IRanges(start = min(x),
                                                       end = max(x)))))
    all(names(out) == names(exp)) & all(ranges(out) == ex.r)
}

test_that("with disjoin = FALSE", {
    expect_true(testfun(c("a", "a", "a", "a"), exp = list(a = 1:5,
                                                          a = 7:11)))
    expect_true(testfun(c("a", "a", "a", "b"), exp = list(a = 1:5,
                                                          a = 7:9,
                                                          b = 10:11)))
    expect_true(testfun(c("a", "a", "b", "b"), exp = list(a = 1:5,
                                                          b = 7:11)))
    expect_true(testfun(c("a", "b", "b", "b"), exp = list(a = 1:4,
                                                          b = 3:5,
                                                          b = 7:11)))
})

test_that("with disjoin = TRUE", {
    disjoin <- TRUE
    expect_true(testfun(c("a", "a", "a", "a"), disjoin, exp = list(a = 1:5,
                                                                   a = 7:11)))
    expect_true(testfun(c("a", "a", "a", "b"), disjoin, exp = list(a = 1:5,
                                                                   a = 7:9,
                                                                   b = 10:11)))
    expect_true(testfun(c("a", "b", "b", "b"), disjoin, exp = list(a = 1:2,
                                                                   b = 5,
                                                                   b = 7:11)))
    expect_true(testfun(c("b", "a", "a", "a"), disjoin, exp = list(b = 1:2,
                                                                   a = 5,
                                                                   a = 7:11)))
})

end(rbg_gr)[1] <- 10
testfun <- testfun

test_that("all unambiguous segments kept", {
    expect_true(testfun(c("a", "b", "b", "b"),
                        exp = list(a = 1:10, b = 3:5, b = 7:11)))
    expect_true(testfun(c("a", "b", "b", "b"), disjoin = TRUE,
                        list(a = 1:2, a = 6, b = 11)))
})

test_that("list/GRangesList inputs", {
    test_list <- function(gn, dj) {
        rbg_grl <- lapply(unique(gn), function(n) rbg_gr[which(gn == n)])
        redux <- reduceByGene(rbg_gr, gene_names = gn, disjoin = dj)
        lredux <- reduceByGene(rbg_grl, gene_names = c("a", "b"), disjoin = dj)
        expect_identical(redux, lredux)
    }
    test_list(c("a", "b", "b", "b"), FALSE)
    test_list(c("a", "a", "b", "b"), FALSE)
    test_list(c("a", "a", "a", "b"), FALSE)
    test_list(c("a", "b", "b", "b"), TRUE)
    test_list(c("a", "a", "b", "b"), TRUE)
    test_list(c("a", "a", "a", "b"), TRUE)
})


# intersectByGene fxn -----------------------------------------------------

context("intersectByGene")

grint <- GRanges(seqnames = "chr1",
                 ranges = IRanges(start = c(1, 3, 6, 7),
                                  end = c(4, 8, 9, 8)),
                 strand = "+")

testfun <- function(gn, exp) {
    out <- intersectByGene(grint, gn)
    ex.r <- Reduce(c, (lapply(exp, function(x) IRanges(start = min(x),
                                                       end = max(x)))))
    all(names(out) == names(exp)) & all(ranges(out) == ex.r)
}

test_that("intersections returned", {
    expect_true(testfun(c("a", "b", "b", "b"), list(a = 1:4, b = 7:8)))
    expect_true(testfun(c("a", "a", "b", "b"), list(a = 3:4, b = 7:8)))
    expect_true(testfun(c("a", "a", "a", "b"), list(b = 7:8)))
    expect_true(testfun(c("a", "b", "a", "b"), list(b = 7:8)))
    expect_true(testfun(c("a", "b", "b", "a"), list(b = 6:8)))
    expect_true(length(intersectByGene(grint, c("a", "a", "a", "a"))) == 0)
})

test_that("list/GRangesList input", {
    test_list <- function(gn) {
        rbg_grl <- lapply(unique(gn), function(n) rbg_gr[which(gn == n)])
        gint <- intersectByGene(rbg_gr, gene_names = gn)
        lgint <- intersectByGene(rbg_grl, gene_names = c("a", "b"))
        expect_identical(gint, lgint)
    }
    test_list(c("a", "b", "b", "b"))
    test_list(c("a", "a", "b", "b"))
    test_list(c("a", "a", "a", "b"))
    test_list(c("a", "b", "b", "b"))
    test_list(c("a", "a", "b", "b"))
    test_list(c("a", "a", "a", "b"))
})




