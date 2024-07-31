context("Generating DESeqDataSets")
library(BRGenomics)
library(DESeq2)

data("PROseq")
data("txs_dm6_chr4")

# split into 6 toy datasets
ps_a_rep1 <- PROseq[seq(1, length(PROseq), 6)]
ps_b_rep1 <- PROseq[seq(2, length(PROseq), 6)]
ps_c_rep1 <- PROseq[seq(3, length(PROseq), 6)]
ps_a_rep2 <- PROseq[seq(4, length(PROseq), 6)]
ps_b_rep2 <- PROseq[seq(5, length(PROseq), 6)]
ps_c_rep2 <- PROseq[seq(6, length(PROseq), 6)]

ps_list <- list(A_rep1 = ps_a_rep1, A_rep2 = ps_a_rep2,
                B_rep1 = ps_b_rep1, B_rep2 = ps_b_rep2,
                C_rep1 = ps_c_rep1, C_rep2 = ps_c_rep2)

# Generate DESeqDataSets --------------------------------------------------

txs <- txs_dm6_chr4[1:50]

dds <- getDESeqDataSet(ps_list, txs, quiet = TRUE, ncores = 1)

test_that("can generate simple DESeqDataSet", {
    expect_is(dds, "DESeqDataSet")
    expect_equivalent(dim(dds), c( length(txs), length(ps_list) ))
    expect_equal(design(dds), ~condition)
    expect_equivalent(rowRanges(dds), txs)
})

test_that("metadata in simple DESeqDataSet correct", {
    expect_equivalent(names(colData(dds)), c("condition", "replicate"))
    expect_equivalent(rownames(colData(dds)), names(ps_list))
    # expect_equivalent(levels(colData(dds)$condition), c("A", "B", "C"))
    expect_equivalent(unique(as.character(colData(dds)$condition)),
                      c("A", "B", "C"))
    # expect_equivalent(levels(colData(dds)$replicate), c("rep1", "rep2"))
    expect_equivalent(unique(as.character(colData(dds)$replicate)),
                      c("rep1", "rep2"))
    expect_equal(names(assays(dds)), "counts")
    expect_null(names(dds))
})

test_that("can use multicore to make DESeqDataSet", {
    expect_is(getDESeqDataSet(ps_list, txs, quiet = TRUE, ncores = 1),
              "DESeqDataSet")
    expect_equivalent(getDESeqDataSet(ps_list, txs, quiet = TRUE, ncores = 1),
                      dds)
})

test_that("can add sizeFactors to DESeqDataSet", {
    expect_null(sizeFactors(dds))
    expect_equivalent(1:6, sizeFactors(getDESeqDataSet(ps_list,
                                                       txs,
                                                       sizeFactors = 1:6,
                                                       quiet = TRUE,
                                                       ncores = 1)))
})

test_that("gene names, discontinuous ranges supported; ordering maintained", {
    dds_dsc <- getDESeqDataSet(ps_list, txs, gene_names = txs$gene_id,
                               quiet = TRUE, ncores = 1)
    expect_is(dds_dsc, "DESeqDataSet")
    expect_equivalent(dim(dds_dsc), c(length(unique(txs$gene_id)),
                                      length(ps_list)))
    expect_equivalent(colSums(assay(dds)), colSums(assay(dds_dsc)))
    expect_equivalent(colData(dds), colData(dds_dsc))
    expect_equivalent(rownames(assay(dds_dsc)), rowRanges(dds_dsc)$gene_id)
    expect_equivalent(rownames(assay(dds_dsc)), unique(txs$gene_id))
})

multi_gr <- mergeGRangesData(ps_list, multiplex = TRUE, ncores = 1)

# cutting tests to try and avoid bioc timeouts...
# test_that("can get dds from multiplexed GRanges", {
#     ddsm <- getDESeqDataSet(multi_gr, txs, quiet = TRUE, ncores = 1)
#     expect_equivalent(dds, ddsm)
#     expect_identical(assay(dds), assay(ddsm))
#     expect_identical(colData(dds), colData(ddsm))
# 
#     # force fields to take
#     expect_equivalent(dds, getDESeqDataSet(multi_gr, txs, quiet = TRUE,
#                                            sample_names = names(ps_list),
#                                            ncores = 1))
# 
#     # don't take all fields
#     multi_AB <- getDESeqDataSet(multi_gr, txs, quiet = TRUE,
#                                 field = names(ps_list)[1:4],
#                                 ncores = 1)
#     expect_equivalent(colData(multi_AB), colData(dds)[1:4, ])
#     expect_equivalent(assay(multi_AB), assay(dds)[, 1:4])
# })

test_that("error if no sample_names found, or don't match input data", {
    expect_error(getDESeqDataSet(list(ps_a_rep1, ps_a_rep2,
                                      ps_b_rep1, ps_b_rep2),
                                 txs, ncores = 1))
    # not enough sample_names (or user should have used 'field' to subset)
    expect_error(getDESeqDataSet(multi_gr, txs, quiet = TRUE,
                                 sample_names = names(ps_list)[1:4],
                                 ncores = 1))
})

test_that("error if sample_names lack replicate specifiers", {
    expect_error(getDESeqDataSet(list(a_1 = ps_a_rep1, a_2 = ps_a_rep2,
                                      b_1 = ps_b_rep1, b_2 = ps_b_rep2),
                                 txs, ncores = 1))
})

test_that("error if gene_names not matched to regions.gr", {
    expect_error(getDESeqDataSet(ps_list, txs,
                                 gene_names = txs$gene_id[-1]))
})


# Getting DESeq2 Results (using reduced dispersion matrices) --------------

context("Get DESeq2 results with pairwise dispersion estimates")

names(dds) <- txs$tx_name
res <- getDESeqResults(dds, contrast.numer = "A", contrast.denom = "B",
                       quiet = TRUE)

test_that("can get single result from simple dds", {
    expect_is(res, "DESeqResults")
    expect_equal(dim(res)[1], length(txs))
    expect_true(all( rownames(res) == names(dds) ))
})

# res2 <- getDESeqResults(dds, comparisons = list(c("B", "A"), c("C", "A")),
#                         ncores = 1, quiet = TRUE)
res2 <- getDESeqResults(dds, comparisons = list(c("B", "A")),
                        ncores = 1, quiet = TRUE)

test_that("get proper results from list of comparisons", {
    expect_is(res2, "list")
    # expect_equivalent(names(res2), c("B_vs_A", "C_vs_A"))
    expect_equivalent(names(res2), c("B_vs_A"))
    expect_is(res2[[1]], "DESeqResults")
    expect_equal(unique(sapply(res2, nrow)), length(txs))
})

test_that("can use multicore to get results from list of comparisons", {
    # cl <- list(c("B", "A"), c("C", "A"))
    cl <- list(c("B", "A"))
    expect_equivalent(res2, getDESeqResults(dds, comparisons = cl, ncores = 1))
})

# test_that("messages when quiet = FALSE", {
#     expect_message(getDESeqResults(dds, contrast.numer = "A",
#                                    contrast.denom = "B", quiet = FALSE))
# })

test_that("error when invalid arguments", {
    expect_error(getDESeqResults(dds,
                                 contrast.numer = "B", contrast.denom = "A",
                                 comparisons = list(c("B", "A"))))
    expect_error(getDESeqResults(dds, comparisons = list("B", "A")))
    expect_error(getDESeqResults(dds, comparisons = c("B", "A")))
    expect_error(getDESeqResults(dds, contrast.numer = "B"))
    expect_error(getDESeqResults(dds, comparisons = list(c(2, 1), c(3, 1))))
    expect_error(getDESeqResults(dds, comparisons = list(c("B", "C", "A"),
                                                         c("C", "B", "A"))))
})

test_that("can apply size factors to whole dds, and comparison-only", {
    # res2sf <- getDESeqResults(dds, comparisons = list(c("B", "A"), c("C", "A")),
    #                           sizeFactors = rep(1, 6), quiet = TRUE, ncores = 1)
    res2sf <- getDESeqResults(dds, comparisons = list(c("B", "A")),
                              sizeFactors = rep(1, 6), quiet = TRUE, ncores = 1)
    expect_false(all(res2[[1]]$log2FoldChange == res2sf[[1]]$log2FoldChange))
    res2sf2 <- getDESeqResults(dds, comparisons = list(c("B", "A")),
                               sizeFactors = rep(1, 4), quiet = TRUE,
                               ncores = 1)
    expect_equivalent(res2sf[[1]]$log2FoldChange, res2sf2[[1]]$log2FoldChange)
})

# test_that("error if size factors the wrong length", {
#     expect_error(getDESeqResults(dds, "B", "A", sizeFactors = rep(1, 5)))
#     expect_error(getDESeqResults(dds, comparisons = list(c("B", "A"),
#                                                          c("C", "A")),
#                                  sizeFactors = rep(1, 5)))
# })

# quiet is only working for list of length 1; the rest get results output...

# test_that("warning if overwriting size factors", {
#     dds_sf <- dds
#     sizeFactors(dds_sf) <- rep(1, 6)
#     expect_warning(getDESeqResults(dds_sf, "B", "A", sizeFactors = rep(1, 6)),
#                    quiet = TRUE, ncores = 1)
#     expect_warning(getDESeqResults(dds_sf, "B", "A", sizeFactors = rep(1, 4)),
#                    quiet = TRUE, ncores = 1)
#     expect_warning(getDESeqResults(dds_sf, comparisons = list(c("B", "A")),
#                                    sizeFactors = rep(1, 6), ncores = 1))
# })

test_that("arguments can be passed to DESeq call", {
    res_alt <- getDESeqResults(dds, contrast.numer = "A", contrast.denom = "B",
                               args.DESeq = list(fitType = "mean"),
                               quiet = TRUE)
    expect_is(res_alt, "DESeqResults")
    expect_equivalent(names(res_alt), names(res))
    expect_equivalent(dim(res_alt), dim(res))
    expect_true(all(res$baseMean == res_alt$baseMean))
    expect_false(all(res$pvalue == res_alt$pvalue))
})

test_that("arguments can be passed to results call", {
    res_alt <- getDESeqResults(dds, contrast.numer = "A", contrast.denom = "B",
                               args.results = list(altHypothesis = "greater"),
                               quiet = TRUE)
    expect_is(res_alt, "DESeqResults")
    expect_equivalent(names(res_alt), names(res))
    expect_equivalent(dim(res_alt), dim(res))
    expect_true(all(res$baseMean == res_alt$baseMean))
    expect_false(all(res$pvalue == res_alt$pvalue))
})

test_that("error if passing arguments not correct", {
    expect_error(getDESeqResults(dds, contrast.numer = "A",
                                 contrast.denom = "B",
                                 args.DESeq = list("mean"), quiet = TRUE))
})
