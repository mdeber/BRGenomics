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
    bw <- rtracklayer::import.bw(ps_p_file)
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


# Stranded bigWig import --------------------------------------------------

ps <- import_bigWig(ps_p_file, ps_m_file, "dm6")

test_that("PROseq files import", {
    expect_is(ps, "GRanges")
})

test_that("Imported PROseq formatted correctly", {
    expect_equal(c(genome(ps), use.names = FALSE), "dm6")
    expect_true(all(width(ps) == 1))
    expect_equal(length(ps), 47380)
    expect_equal(sum(score(ps)), 73887)
    expect_is(ps$score, "integer")
    expect_equal(as.character(unique(seqnames(ps))), "chr4")
    expect_true(all(as.character(strand(ps)) %in% c("+", "-")))
})


# Stranded bedGraph import ------------------------------------------------

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

ps_paired <- import_bedGraph(paired_p_file, paired_m_file, "dm6")

test_that("Paired PROseq files import", {
    expect_is(ps_paired, "GRanges")
})

test_that("Paired PROseq formatted correctly", {
    expect_equal(c(genome(ps_paired), use.names = FALSE), "dm6")
    expect_equal(length(ps_paired), 53179)
    expect_equal(sum(score(ps_paired)), 73887)
    expect_is(ps_paired$score, "integer")
    expect_equal(as.character(unique(seqnames(ps_paired))), "chr4")
    expect_true(all(as.character(strand(ps_paired)) %in% c("+", "-")))
})



# Bam file import ---------------------------------------------------------

ps_bam <- system.file("extdata", "PROseq_dm6_chr4.bam", package = "BRGenomics")

test_that("Bam file is found", {
    expect_is(ps_bam, "character")
    expect_true(nchar(ps_bam) > 1)
})

test_that("Bam file imports", {
    expect_is(import_bam(ps_bam, paired_end = FALSE), "GRanges")
})

test_that("Options applied", {
    ps_reads <- import_bam(ps_bam, revcomp = TRUE, paired_end = FALSE)
    expect_equal(unique(width(import_bam(ps_bam, trim.to = "5p",
                                         paired_end = FALSE))), 1)
    expect_equal(unique(width(import_bam(ps_bam, trim.to = "3p",
                                         paired_end = FALSE))), 1)
    expect_equal(unique(width(import_bam(ps_bam, trim.to = "center",
                                         paired_end = FALSE))), 1)
    expect_true(isDisjoint(import_bam(ps_bam, trim.to = "5p",
                                      paired_end = FALSE)))
    expect_equal(sum(as.character(strand(ps_reads)) == "+"),
                 sum(as.character(strand(
                     import_bam(ps_bam, paired_end = FALSE))) == "-"))
    ps_nostrand <- import_bam(ps_bam, ignore.strand = TRUE, paired_end = FALSE)
    expect_true(length(ps_nostrand) < length(ps_reads))
    expect_equal(unique(as.character(strand(ps_nostrand))), "*")
    expect_equal(sum(score(ps_nostrand)), sum(score(ps_reads)))
    expect_equal(sum(score(ps_reads)),
                 sum(score(import_bam(ps_bam, trim.to = "3p",
                                      paired_end = FALSE))))
})

