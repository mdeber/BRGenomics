library(BRGenomics)

# Check included data -----------------------------------------------------

context("Included data files")

test_that("Internal RData present", {
    data("PROseq")
    expect_is(PROseq, "GRanges")
    data("PROseq_paired")
    expect_is(PROseq_paired, "GRanges")
})

ps_p_file <- system.file("extdata", "PROseq_dm6_chr4_plus.bw",
                         package = "BRGenomics")
ps_m_file <- system.file("extdata", "PROseq_dm6_chr4_minus.bw",
                         package = "BRGenomics")

test_that("bigWig files (PROseq files) are found", {
    expect_is(ps_p_file, "character")
    expect_is(ps_m_file, "character")
    expect_true(nchar(ps_p_file) > 1)
    expect_true(nchar(ps_m_file) > 1)
})

paired_p_file <- system.file("extdata", "PROseq_dm6_chr4_plus.bedGraph",
                             package = "BRGenomics")
paired_m_file <- system.file("extdata", "PROseq_dm6_chr4_minus.bedGraph",
                             package = "BRGenomics")

test_that("bedGraph files (paired PROseq files) are found", {
    expect_is(paired_p_file, "character")
    expect_is(paired_m_file, "character")
    expect_true(nchar(paired_p_file) > 1)
    expect_true(nchar(paired_m_file) > 1)
})

ps_bam <- system.file("extdata", "PROseq_dm6_chr4.bam", package = "BRGenomics")

test_that("Bam file is found", {
    expect_is(ps_bam, "character")
    expect_true(nchar(ps_bam) > 1)
})

# Check tidyChromosome function -------------------------------------------

context("tidyChromosome")

bw <- GRanges(seqnames = c("chr2L", "chrM", "chrX", "chrY", "unassigned"),
              ranges = IRanges(1:5, 2:6))
genome(bw) <- "dm6"

test_that("sex chromosomes", {
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
})

test_that("trim non-standard, keep sex", {
    # default keeps sex, trims M and non-standard
    trim_default <- tidyChromosomes(bw, keep.X = TRUE, keep.Y = TRUE,
                                    keep.M = FALSE,
                                    keep.nonstandard = FALSE)
    expect_equal(length(trim_default), length(bw)-2)
    expect_equal(seqinfo(trim_default),
                 sortSeqlevels(dropSeqlevels(seqinfo(bw),
                                             c("chrM", "unassigned"))))
})

test_that("use pre-assigned genome", {
    expect_identical(tidyChromosomes(bw),
                     tidyChromosomes(bw, genome = "dm6"))
})

test_that("trim from list input", {
    trim_ls <- tidyChromosomes(list(bw, bw))
    expect_is(trim_ls, "list")
    expect_identical(trim_ls[[1]], tidyChromosomes(bw))
})

test_that("trim seqinfo object", {
    si <- seqinfo(bw)
    expect_identical(tidyChromosomes(si), seqinfo(tidyChromosomes(bw)))
})


# Stranded bigWig import --------------------------------------------------

context("bigWig/bedGraph import")

if (.Platform$OS.type == "unix") {

    ps <- import_bigWig(ps_p_file, ps_m_file, "dm6")

    test_that("PROseq files import", {
        expect_is(ps, "GRanges")
    })

    test_that("Imported PROseq formatted correctly", {
        expect_equal(c(genome(ps), use.names = FALSE), "dm6")
        expect_true(all(width(ps) == 1))
        expect_equal(length(ps), 150)
        expect_equal(sum(score(ps)), 190)
        expect_is(ps$score, "integer")
        expect_equal(as.character(unique(seqnames(ps))), "chr4")
        expect_true(all(as.character(strand(ps)) %in% c("+", "-")))
    })

    test_that("can import from lists of files", {
        psl <- import_bigWig(c(ps_p_file, ps_p_file),
                             c(ps_m_file, ps_m_file),
                             genome = "dm6", ncores = 1)
        expect_is(psl, "list")
        expect_equal(length(psl), 2L)
        expect_identical(psl[[1]], ps)

        # must be same length
        expect_error(import_bigWig(c(ps_p_file, ps_p_file),
                                   ps_m_file,
                                   genome = "dm6", ncores = 1))
    })

    test_that("can import single strands alone", {
        psp <- import_bigWig(c(ps_p_file, ps_p_file), genome = "dm6",
                             ncores = 1)
        expect_equal("+", levels(droplevels(strand(psp[[1]]))))

        psm <- import_bigWig(NULL, c(ps_m_file, ps_m_file), genome = "dm6",
                             ncores = 1)
        expect_equal("-", levels(droplevels(strand(psm[[1]]))))
    })

}

# Stranded bedGraph import ------------------------------------------------


ps_paired <- import_bedGraph(paired_p_file, paired_m_file, "dm6")

test_that("Paired PROseq files import", {
    expect_is(ps_paired, "GRanges")
})

test_that("Paired PROseq formatted correctly", {
    expect_equal(c(genome(ps_paired), use.names = FALSE), "dm6")
    expect_equal(length(ps_paired), 164)
    expect_equal(sum(score(ps_paired)), 190)
    expect_is(ps_paired$score, "integer")
    expect_equal(as.character(unique(seqnames(ps_paired))), "chr4")
    expect_true(all(as.character(strand(ps_paired)) %in% c("+", "-")))
})

test_that("can import from lists of files", {
    pairedl <- import_bedGraph(c(paired_p_file, paired_p_file),
                           c(paired_m_file, paired_m_file),
                         genome = "dm6", ncores = 1)
    expect_is(pairedl, "list")
    expect_equal(length(pairedl), 2L)
    expect_identical(pairedl[[1]], ps_paired)

    # must be same length
    expect_error(import_bedGraph(paired_p_file,
                                 c(paired_m_file, paired_m_file),
                                 genome = "dm6", ncores = 1))
})

test_that("can import single strands alone", {
    pairedp <- import_bedGraph(c(paired_p_file, paired_p_file),
                               genome = "dm6", ncores = 1)
    expect_equal("+", levels(droplevels(strand(pairedp[[1]]))))

    pairedm <- import_bedGraph(NULL, c(paired_m_file, paired_m_file),
                               genome = "dm6", ncores = 1)
    expect_equal("-", levels(droplevels(strand(pairedm[[1]]))))
})


# Bam file import ---------------------------------------------------------

context("Bam file import")

ps_reads <- import_bam(ps_bam, revcomp = TRUE, paired_end = FALSE)

test_that("Bam file imports", {
    expect_is(import_bam(ps_bam, paired_end = FALSE), "GRanges")
})

test_that("trim.to option", {
    expect_equal(unique(width(import_bam(ps_bam, trim.to = "5p",
                                         paired_end = FALSE))), 1)
    expect_equal(unique(width(import_bam(ps_bam, trim.to = "3p",
                                         paired_end = FALSE))), 1)
    expect_equal(unique(width(import_bam(ps_bam, trim.to = "center",
                                         paired_end = FALSE))), 1)
    # makes disjoint
    expect_true(isDisjoint(import_bam(ps_bam, trim.to = "5p",
                                      paired_end = FALSE)))
})

test_that("reads counted, regardless of trim.to", {
    expect_equal(sum(score(ps_reads)),
                 sum(score(import_bam(ps_bam, trim.to = "3p",
                                      paired_end = FALSE))))
})


test_that("revcomp option", {
    no_revcomp <- import_bam(ps_bam, paired_end = FALSE)
    expect_equal(sum(as.character(strand(ps_reads)) == "+"),
                 sum(as.character(strand(no_revcomp)) == "-"))
})

test_that("unstranded option", {
    ps_nostrand <- import_bam(ps_bam, ignore.strand = TRUE, paired_end = FALSE)
    expect_true(length(ps_nostrand) <= length(ps_reads))
    expect_equal(unique(as.character(strand(ps_nostrand))), "*")
    expect_equal(sum(score(ps_nostrand)), sum(score(ps_reads)))
})

test_that("shift option", {
    ps_paired <- import_bam(ps_bam, revcomp = TRUE, paired_end = FALSE,
                            shift = c(0, -1))
    expect_identical(resize(ps_reads, 1, "start"),
                     resize(ps_paired, 1, "start"))

    reads.3p <- start(resize(ps_reads, 1, "end"))
    paired.3p <- start(resize(ps_paired, 1, "end"))
    pranges <- as.logical(strand(ps_reads) == "+")
    expect_identical(reads.3p[pranges], paired.3p[pranges] + 1L)
    expect_identical(reads.3p[!pranges], paired.3p[!pranges] - 1L)

    expect_error(import_bam(ps_bam, shift = c(1, 2, 3)))
})

test_that("paired-end test", {
    expect_identical(import_bam(ps_bam, revcomp = TRUE),
                     ps_reads)
})

test_that("chunking", {
    expect_identical(import_bam(ps_bam, revcomp = TRUE, paired_end = FALSE,
                                yieldSize = 10),
                     ps_reads)
})

test_that("import from list", {
    ps_reads_ls <- import_bam(list(ps_bam, ps_bam), revcomp = TRUE,
                              paired_end = FALSE, ncores = 1)
    expect_is(ps_reads_ls, "list")
    expect_equal(length(ps_reads_ls), 2)
    expect_identical(ps_reads_ls[[1]], ps_reads)
})


test_that("bam helper fxns evaluate", {
    # args are already tested; check that fxns exists and pass args
    ps <- import_bam_PROseq(ps_bam)
    expect_true(isBRG(ps))

    ps5 <- import_bam_PROcap(ps_bam)
    expect_true(isBRG(ps5))

    ps_atac <- import_bam_ATACseq(ps_bam, paired_end = FALSE)
    expect_true(is(ps_atac, "GRanges"))
    expect_true(!isBRG(ps_atac))
})
