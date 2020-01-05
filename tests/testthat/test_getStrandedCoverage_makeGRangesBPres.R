# Testing getStrandedCoverage ---------------------------------------------

context("getting stranded coverage")
library(BRGenomics)

data("PROseq_paired")
paired_3p_ends <- resize(PROseq_paired, 1, fix = "end")
paired_3p_cov <- getStrandedCoverage(paired_3p_ends)


test_that("Stranded coverage makes stranded GRanges", {
    expect_is(paired_3p_cov, "GRanges")
    expect_is(getStrandedCoverage(PROseq_paired), "GRanges")
    expect_true(all(as.character(strand(paired_3p_cov)) %in% c("+", "-")))
})

test_that("Stranded coverage output not single-width", {
    expect_true(all( width(paired_3p_ends) == 1 ))
    expect_false(all( width(paired_3p_cov) == 1 ))
})

test_that("Stranded coverage is disjoint", {
    expect_false(isDisjoint(paired_3p_ends))
    expect_true(isDisjoint(paired_3p_cov))
})

# signal weighted by width represents total reads counted
sum_3pcov <- sum(score(paired_3p_cov))
sum_3pcov_by_width <- sum( score(paired_3p_cov) * width(paired_3p_cov) )

test_that("Stranded coverage correctly calculated", {
    expect_equal(length(paired_3p_ends), 53179)
    expect_equal(length(paired_3p_cov), 43385)

    expect_equal(sum(score(paired_3p_ends)), 73887)
    expect_equal(sum_3pcov, 69349)
    expect_equal(sum_3pcov_by_width, sum(score(paired_3p_ends)))
})

test_that("stranded coverage correct without any weighting", {
    reads_gr <- rep(paired_3p_ends, paired_3p_ends$score)
    mcols(reads_gr) <- NULL
    reads_cov <- getStrandedCoverage(reads_gr, field = NULL)

    expect_equivalent(reads_cov, paired_3p_cov)
    expect_error(getStrandedCoverage(reads_gr))
})

test_that("Empty input returns GRanges", {
    null_cov <- getStrandedCoverage(paired_3p_ends[0])
    expect_is(null_cov, "GRanges")
    expect_equal(length(null_cov), 0)
    expect_equivalent(seqinfo(paired_3p_ends), seqinfo(null_cov))
    expect_equivalent(mcols(paired_3p_ends[0]), mcols(null_cov))
})

test_that("error on invalid field", {
    expect_error(getStrandedCoverage(paired_3p_ends, "wrongname"))
})


# Test makeGRangesBRG -----------------------------------------------------

context("making GRanges BP-resolution")

test_that("disJoint input for makeGRangesBRG makes error", {
    expect_error(makeGRangesBRG(PROseq_paired))
    expect_error(makeGRangesBRG(paired_3p_ends))
})

test_that("Making GRangesBPres doesn't affect BPres GRanges object", {
    data("PROseq")
    expect_identical(makeGRangesBRG(PROseq), PROseq)
})

paired_3p_bpres <- makeGRangesBRG(paired_3p_cov)

test_that("Can make single-width GRanges", {
    expect_is(paired_3p_bpres, "GRanges")
    expect_true(isDisjoint(paired_3p_bpres))
    expect_true(all(width(paired_3p_bpres) == 1))
})

test_that("BPresGRanges preserves input information", {
    expect_equal(length(paired_3p_bpres), sum(width(paired_3p_cov)))
    expect_equal(sum_3pcov_by_width, sum(score(paired_3p_bpres)))
})








