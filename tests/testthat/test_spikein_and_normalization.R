context("counting spike-in reads")
library(BRGenomics)

# make a list of datasets with spike-in chromosomes -----------------------

gr1_rep1 <- GRanges(seqnames = c("chr1", "chr2", "spikechr1", "spikechr2"),
                    ranges = IRanges(start = 1:4, width = 1),
                    strand = "+")
gr3_rep2 <- gr3_rep1 <- gr2_rep2 <- gr2_rep1 <- gr1_rep2 <- gr1_rep1


# exp + spike = total
score(gr1_rep1) <- 1 # 2 + 2 = 4
score(gr2_rep1) <- 2 # 4 + 4 = 8
score(gr3_rep1) <- c(2, 2, 1, 1) # 4 + 2 = 6

score(gr1_rep2) <- c(1, 1, 2, 1) # 2 + 3 = 5
score(gr2_rep2) <- c(2, 3, 2, 2) # 5 + 4 = 9
score(gr3_rep2) <- c(4, 4, 2, 2) # 8 + 4 = 12

grl <- list(gr1_rep1, gr2_rep1, gr3_rep1,
            gr1_rep2, gr2_rep2, gr3_rep2)
names(grl) <- c("gr1_rep1", "gr2_rep1", "gr3_rep1",
                "gr1_rep2", "gr2_rep2", "gr3_rep2")


# test getting spike-in counts --------------------------------------------

s_counts <- getSpikeInCounts(grl, si_names = c("spikechr1", "spikechr2"),
                             ncores = 2)
sr <- s_counts$spike_reads
er <- s_counts$exp_reads

test_that("can get read counts from given spike-in chromosomes", {
    expect_is(s_counts, "data.frame")
    expect_equal(s_counts$sample, names(grl))
    expect_equal(s_counts$total_reads, c(4, 8, 6, 5, 9, 12))
    expect_equal(er, c(2, 4, 4, 2, 5, 8))
    expect_equal(sr, c(2, 4, 2, 3, 4, 4))
})

test_that("can get read counts from pattern-matched spike-in chromosomes", {
    expect_equivalent(s_counts,
                      getSpikeInCounts(grl, si_pattern = "spike", ncores = 2))
})

test_that("can get spike-in counts for GRanges input", {
    expect_equivalent(s_counts[1, ],
                      getSpikeInCounts(gr1_rep1, si_pattern = "spike",
                                       ncores = 1))
})


test_that("can get spike-in counts for NULL field", {
    null_counts <- getSpikeInCounts(grl, si_pattern = "spike", field = NULL,
                                    ncores = 2)
    expect_equivalent(s_counts$sample, null_counts$sample)
    expect_true(all(null_counts$total_reads == 4))
    expect_true(all(null_counts$exp_reads == 2))
})

test_that("can get spike-in counts for multiple fields", {
    grl_mix <- grl
    grl_mix$gr1_rep1$newscore <- 2 # switching gr1_rep1 and gr2_rep1
    grl_mix$gr2_rep1$newscore <- 1
    mix_counts <- getSpikeInCounts(grl_mix, si_pattern = "spike",
                                   field = c("newscore", "newscore", "score",
                                             "score", "score", "score"),
                                   ncores = 2)
    expect_equivalent(s_counts[1:2, -1], mix_counts[2:1, -1])
})


# test subsetting by spike-in reads ---------------------------------------

context("subsetting GRanges by spike-in reads")

test_that("can filter out spike-in reads", {
    gr <- removeSpikeInReads(gr1_rep1, si_pattern = "spike", ncores = 1)
    expect_equivalent(gr, gr1_rep1[1:2])

    trimlist <- removeSpikeInReads(grl[1:2], si_pattern = "spike", ncores = 2)
    expect_is(trimlist, "list")
    expect_equivalent(trimlist[[1]], gr)
    expect_equivalent(trimlist[[2]], gr2_rep1[1:2])
})


test_that("can isolate spike-in reads", {
    gr <- getSpikeInReads(gr1_rep1, si_pattern = "spike", ncores = 1)
    expect_equivalent(gr, gr1_rep1[3:4])

    trimlist <- getSpikeInReads(grl[1:2], si_pattern = "spike", ncores = 2)
    expect_is(trimlist, "list")
    expect_equivalent(trimlist[[1]], gr)
    expect_equivalent(trimlist[[2]], gr2_rep1[3:4])
})


# test normalization factor calculations ----------------------------------

context("normalization factor calculations")

test_that("RPM normalization works", {
    nf_rpm <- getSpikeInNFs(grl, si_pattern = "spike", method = "RPM",
                            ncores = 2)
    expect_is(nf_rpm, "numeric")
    expect_true(all(nf_rpm * er == 1e6))
})

test_that("simple spike-in read normalization", {
    nf_snr <- getSpikeInNFs(grl, si_pattern = "spike", method = "SNR",
                            batch_norm = FALSE, ncores = 2)
    expect_is(nf_snr, "numeric")
    expect_true(all(nf_snr * sr == min(sr)))
})

test_that("batch-corrected spike-in read normalization", {
    nf_snrb <- getSpikeInNFs(grl, si_pattern = "spike", ctrl_pattern = "gr1",
                             method = "SNR", batch_norm = TRUE, ncores = 2)
    expect_is(nf_snrb, "numeric")

    # in each replicate, spike-in normalization
    expect_true(all(nf_snrb[1:3] * sr[1:3] == min(sr[1:3])))
    expect_true(all(nf_snrb[4:6] * sr[4:6] == min(sr[4:6])))
    # spike-in reads across replicates are independent
    expect_true(nf_snrb[1] * sr[1] != nf_snrb[4] * sr[4])

    # but negative control normalized readcounts are the same
    expect_equal(nf_snrb[1] * er[1], nf_snrb[4] * er[4])
})


si_ratio <- er / sr # spike-in normalized read counts

test_that("spike-in normalized RPM vs. control", {
    nf_srpmcb <- getSpikeInNFs(grl, si_pattern = "spike", ctrl_pattern = "gr1",
                                method = "SRPMC", batch_norm = TRUE,
                                ncores = 2)

    expect_true(all(er[c(1,4)] * nf_srpmcb[c(1,4)] == 1e6))

    expect_equivalent(si_ratio[1:3] / si_ratio[1],
                      nf_srpmcb[1:3] * er[1:3] / 1e6)
    expect_equivalent(si_ratio[4:6] / si_ratio[4],
                      nf_srpmcb[4:6] * er[4:6] / 1e6)
})


test_that("spike-in normalized RPM vs. control, no batch normalization", {
    # without batch normalization; takes average number of reads across controls
    nf_srpmc <- getSpikeInNFs(grl, si_pattern = "spike", ctrl_pattern = "gr1",
                               method = "SRPMC", batch_norm = FALSE,
                               ncores = 2)

    expect_true(all(er[c(1,4)] * nf_srpmc[c(1,4)] != 1e6))
    expect_equal(mean( er[c(1,4)] * nf_srpmc[c(1,4)] ), 1e6)

    expect_equivalent(si_ratio / mean(si_ratio[c(1, 4)]),
                      nf_srpmc * er / 1e6)
})


test_that("error if give names and regex for control samples", {
    expect_error(getSpikeInNFs(grl, ctrl_pattern = "gr1",
                               ctrl_names = c("gr1_rep1", "gr1_rep2")))
})






