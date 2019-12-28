context("Data import")
library(BRGenomics)

ps_p_file <- system.file("extdata", "PROseq_dm6_chr4_plus.bw",
                         package = "BRGenomics")
ps_m_file <- system.file("extdata", "PROseq_dm6_chr4_minus.bw",
                         package = "BRGenomics")


# Check internal datasets -------------------------------------------------

test_that("Internal RData present", {
    data("PROseq")
    expect_is(PROseq, "GRanges")
    data("PROseq_paired")
    expect_is(PROseq_paired, "GRanges")
})

# Check raw data files ----------------------------------------------------

test_that("PROseq files are found", {
    expect_is(ps_p_file, "character")
    expect_is(ps_m_file, "character")
    expect_true(nchar(ps_p_file) > 1)
    expect_true(nchar(ps_m_file) > 1)
})


# Check tidyChromosome function -------------------------------------------

test_that("tidyChromosome works", {
    bw <- import.bw(ps_p_file)
    seqlevels(bw) <- c("chr4", "chr2L", "chrM", "chrX", "chrY", "unassigned")
    genome(bw) <- "dm6"
    seqnames(bw)[1:5] <- c("chr2L", "chrM", "chrX", "chrY", "unassigned")

    trim_x <- tidyChromosomes(bw, keep.X = FALSE, keep.M = TRUE,
                              keep.nonstandard = TRUE)
    expect_equal(length(trim_x), length(bw)-1)
    expect_equal(seqinfo(trim_x),
                 sortSeqlevels(dropSeqlevels(seqinfo(bw), "chrX")))

    trim_y <- tidyChromosomes(bw, keep.Y = FALSE, keep.M = TRUE,
                              keep.nonstandard = TRUE)
    expect_equal(length(trim_y), length(bw)-1)
    expect_equal(seqinfo(trim_y),
                 sortSeqlevels(dropSeqlevels(seqinfo(bw), "chrY")))

    # default keeps sex, trims M and non-standard
    trim_default <- tidyChromosomes(bw, keep.X = TRUE, keep.Y = TRUE,
                                    keep.M = FALSE, keep.nonstandard = FALSE)
    expect_equal(length(trim_default), length(bw)-2)
    expect_equal(seqinfo(trim_default),
                 sortSeqlevels(dropSeqlevels(seqinfo(bw),
                                             c("chrM", "unassigned"))))
})


# Generic bigWig import ---------------------------------------------------

test_that("Single bigWig can be imported", {
    bw <- import.bw_trim(ps_p_file, "dm6")
    expect_equal(c(genome(bw), use.names = FALSE), "dm6")
    expect_equal(length(bw), 16587)
    expect_is(bw$score, "integer")
    expect_true(all(as.character(strand(bw)) == "*"))
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





