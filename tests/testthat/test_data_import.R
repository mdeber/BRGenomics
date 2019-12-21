context("Data Import")
library(BRGenomics)

# PRO-seq import

ps_p_file <- system.file("extdata", "PROseq_dm6_chr4_plus.bw",
                         package = "BRGenomics")
ps_m_file <- system.file("extdata", "PROseq_dm6_chr4_minus.bw",
                         package = "BRGenomics")

test_that("PROseq files are found", {
    expect_true(length(ps_p_file) > 1)
    expect_true(length(ps_m_file) > 1)
})

ps <- import.PROseq(ps_p_file, ps_m_file, "dm6")

test_that("PROseq files import", {
    expect_is(ps, "GRanges")
})

test_that("Imported PROseq formatted correctly", {
    expect_equal(genome(ps), "dm6")
    expect_true(all(width(ps) == 1))
    expect_equal(length(ps), 47533)
    expect_equal(sum(score(ps)), 74157)
    expect_is(ps$score, "integer")
    expect_equal(unique(seqnames(ps)), "chr4")
    expect_true(all(strand(ps) %in% c("+", "-")))
})

# CoPRO import

paired_p_file <- system.file("extdata", "PROseq_dm6_chr4_plus.bedGraph",
                             package = "BRGenomics")
paired_m_file <- system.file("extdata", "PROseq_dm6_chr4_minus.bedGraph",
                             package = "BRGenomics")

test_that("Paired PROseq files are found", {
    expect_true(length(paired_p_file) > 1)
    expect_true(length(paired_m_file) > 1)
})

ps_paired <- import.CoPRO(paired_p_file, paired_m_file, "dm6")

test_that("Paired PROseq files import", {
    expect_is(ps_paired, "GRanges")
})


test_that("Paired PROseq formatted correctly", {
    expect_equal(genome(ps_paired), "dm6")
    expect_equal(length(ps_paired), 52464)
    expect_equal(sum(score(ps_paired)), 73011)
    expect_is(ps_paired$score, "integer")
    expect_equal(unique(seqnames(ps)), "chr4")
    all(strand(ps_paired) %in% c("+", "-"))
})





