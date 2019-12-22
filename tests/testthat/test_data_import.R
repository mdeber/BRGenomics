context("Data Import")
library(BRGenomics)



ps_p_file <- system.file("extdata", "PROseq_dm6_chr4_plus.bw",
                         package = "BRGenomics")
ps_m_file <- system.file("extdata", "PROseq_dm6_chr4_minus.bw",
                         package = "BRGenomics")


# Check raw data files ----------------------------------------------------


test_that("PROseq files are found", {
    expect_is(ps_p_file, "character")
    expect_is(ps_m_file, "character")
    expect_true(nchar(ps_p_file) > 1)
    expect_true(nchar(ps_m_file) > 1)
})


# Generic bigWig import ---------------------------------------------------


test_that("Single bigWig can be imported", {
    bw <- import.bw_trim(ps_p_file, "dm6")
    expect_equal(c(genome(bw), use.names = FALSE), "dm6")
    expect_equal(length(bw), 16587)
    expect_is(bw$score, "integer")
    expect_true(all(as.character(strand(ps)) == "*"))
})


# PRO-seq import ----------------------------------------------------------


ps <- import.PROseq(ps_p_file, ps_m_file, "dm6")

test_that("PROseq files import", {
    expect_is(ps, "GRanges")
})

test_that("Imported PROseq formatted correctly", {
    expect_equal(c(genome(ps), use.names = FALSE), "dm6")
    expect_true(all(width(ps) == 1))
    expect_equal(length(ps), 47533)
    expect_equal(sum(score(ps)), 74157)
    expect_is(ps$score, "integer")
    expect_equal(as.character(unique(seqnames(ps))), "chr4")
    expect_true(all(as.character(strand(ps)) %in% c("+", "-")))
})


# Paired ends import ------------------------------------------------------


paired_p_file <- system.file("extdata", "PROseq_dm6_chr4_plus.bedGraph",
                             package = "BRGenomics")
paired_m_file <- system.file("extdata", "PROseq_dm6_chr4_minus.bedGraph",
                             package = "BRGenomics")

test_that("Paired PROseq files are found", {
    expect_is(paired_p_file, "character")
    expect_is(paired_m_file, "character")
    expect_true(nchar(paired_p_file) > 1)
    expect_true(nchar(paired_m_file) > 1)
})

ps_paired <- import.CoPRO(paired_p_file, paired_m_file, "dm6")

test_that("Paired PROseq files import", {
    expect_is(ps_paired, "GRanges")
})


test_that("Paired PROseq formatted correctly", {
    expect_equal(c(genome(ps_paired), use.names = FALSE), "dm6")
    expect_equal(length(ps_paired), 52464)
    expect_equal(sum(score(ps_paired)), 73011)
    expect_is(ps_paired$score, "integer")
    expect_equal(as.character(unique(seqnames(ps_paired))), "chr4")
    expect_true(all(as.character(strand(ps_paired)) %in% c("+", "-")))
})





