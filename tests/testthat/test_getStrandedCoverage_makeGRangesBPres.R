Context("Functions on dataset GRanges")
library(BRGenomics)

data("PROseq_paired")
paired_3p_ends <- resize(PROseq_paired, 1, fix = "end")
paired_3p_cov <- getStrandedCoverage(paired_3p_ends)


# Testing getStrandedCoverage ---------------------------------------------

testthat("Stranded coverage makes stranded GRanges", {
    expect_is(paired_3p_cov, "GRanges")
    expect_is(getStrandedCoverage(PROseq_paired), "GRanges")
    expect_true(all(strand(paired_3p_cov) %in% c("+", "-")))
})

testthat("Stranded coverage output not single-width", {
    expect_true(all( width(paired_3p_ends) == 1 ))
    expect_false(all( width(paired_3p_cov) == 1 ))
})

testthat("Stranded coverage is disjoint", {
    expect_false(isDisjoint(paired_3p_ends))
    expect_true(isDisjoint(paired_3p_cov))
})

# signal weighted by width represents total reads counted
sum_3pcov <- sum(score(paired_3p_cov))
sum_3pcov_by_width <- sum( score(paired_3p_cov) * width(paired_3p_cov) )

testthat("Stranded coverage correctly calculated", {
    expect_equal(length(paired_3p_ends), 52464)
    expect_equal(length(paired_3p_cov), 42916)

    expect_equal(sum(score(paired_3p_ends)), 73011)
    expect_equal(sum_3pcov, 68562)
    expect_equal(sum_3pcov_by_width, 73011)
})



# Testing makeGRangesBPres ------------------------------------------------

testthat("disJoint input for makeGRangesBPres makes error", {
    expect_error(makeGRangesBPres(PROseq_paired))
    expect_error(makeGRangesBPres(paired_3p_ends))
})

testthat("Making GRangesBPres doesn't affect BPres GRanges object", {
    data("PROseq")
    expect_identical(makeGRangesBPres(PROseq), PROseq)
})

paired_3p_bpres <- makeGRangesBPres(paired_3p_cov)

testthat("Can make single-width GRanges", {
    expect_is(paired_3p_bpres, "GRanges")
    expect_true(isDisjoint(paired_3p_bpres))
    expect_true(all(width(paired_3p_bpres) == 1))
})

testthat("BPresGRanges preserves input information", {
    expect_equal(length(paired_3p_bpres), sum(width(paired_3p_cov)))
    expect_equal(sum_3pcov_by_width, sum(score(paired_3p_bpres)))
})








