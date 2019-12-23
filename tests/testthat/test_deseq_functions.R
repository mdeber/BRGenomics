context("Generating DESeqDataSets")
library(BRGenomics)
library(DESeq2)

data("PROseq")
data("txs_dm6_chr4")

# split into 4 toy datasets
ps_a_rep1 <- PROseq[seq(1, length(PROseq), 6)]
ps_b_rep1 <- PROseq[seq(2, length(PROseq), 6)]
ps_c_rep1 <- PROseq[seq(3, length(PROseq), 6)]
ps_a_rep2 <- PROseq[seq(4, length(PROseq), 6)]
ps_b_rep2 <- PROseq[seq(5, length(PROseq), 6)]
ps_c_rep2 <- PROseq[seq(6, length(PROseq), 6)]

ps_list <- list(A_rep1 = ps_a_rep1,
                A_rep2 = ps_a_rep2,
                B_rep1 = ps_b_rep1,
                B_rep2 = ps_b_rep2,
                C_rep1 = ps_c_rep1,
                C_rep2 = ps_c_rep2)


# Generate DESeqDataSets --------------------------------------------------

txs <- txs_dm6_chr4[1:50]

dds <- getDESeqDataSet(ps_list, txs, quiet = TRUE, ncores = 1)

test_that("can generate simple DESeqDataSet", {
    expect_is(dds, "DESeqDataSet")
    expect_equivalent(dim(dds), c( length(txs), length(ps_list) ))
    expect_equal(dds@design, ~condition)
    expect_equivalent(dds@rowRanges, txs)
})

test_that("metadata in simple DESeqDataSet correct", {
    expect_equivalent(names(dds@colData), c("condition", "replicate"))
    expect_equivalent(rownames(dds@colData), names(ps_list))
    expect_equivalent(levels(dds@colData$condition), c("A", "B", "C"))
    expect_equivalent(levels(dds@colData$replicate), c("rep1", "rep2"))
    expect_equal(names(dds@assays), "counts")
    expect_null(dds@NAMES)
})

test_that("can use multicore to make DESeqDataSet", {
    expect_is(getDESeqDataSet(ps_list, txs, quiet = T), "DESeqDataSet")
    expect_equivalent(getDESeqDataSet(ps_list, txs, quiet = T), dds)
})

test_that("can add sizeFactors to DESeqDataSet", {
    expect_null(sizeFactors(dds))
    expect_equivalent(1:6, sizeFactors(getDESeqDataSet(ps_list,
                                                       txs,
                                                       sizeFactors = 1:6,
                                                       quiet = T)))
})

test_that("gene names and discontinuous ranges supported", {
    dds_dsc <- getDESeqDataSet(ps_list, txs,
                               gene_names = txs$gene_id, quiet = T)
    expect_is(dds_dsc, "DESeqDataSet")
    expect_equivalent(dim(dds_dsc), c(length(unique(txs$gene_id)),
                                      length(ps_list)))
    expect_equivalent(colSums(assay(dds)), colSums(assay(dds_dsc)))
    expect_equivalent(dds@colData, dds_dsc@colData)
    expect_equivalent(rownames(assay(dds_dsc)),
                      unique(txs$gene_id))
})

test_that("error if no sample_names found", {
    expect_error(getDESeqDataSet(list(ps_a_rep1, ps_a_rep2,
                                      ps_b_rep1, ps_b_rep2),
                                 txs, ncores = 1))
})

test_that("error if sample_names lack replicate specifiers", {
    expect_error(getDESeqDataSet(list(a_1 = ps_a_rep1,
                                      a_2 = ps_a_rep2,
                                      b_1 = ps_b_rep1,
                                      b_2 = ps_b_rep2),
                                 txs, ncores = 1))
})

test_that("error if gene_names not matched to regions.gr", {
    expect_error(getDESeqDataSet(ps_list, txs,
                                 gene_names = txs$gene_id[-1]))
})


test_that("message if quiet = FALSE", {
    expect_message(getDESeqDataSet(ps_list, txs, quiet = FALSE))
})



# Getting DESeq2 Results (using reduced dispersion matrices) --------------


context("Get DESeq2 results with pairwise dispersion estimates")

names(dds) <- txs$tx_name
res <- getDESeqResults(dds,
                       contrast.numer = "A",
                       contrast.denom = "B",
                       quiet = TRUE)

test_that("can get single result from simple dds", {
    expect_is(res, "DESeqResults")
    expect_equal(dim(res)[1], length(txs))
    expect_true(all( rownames(res) == names(dds) ))
})

res_2 <- getDESeqResults(dds,
                         comparisons.list = list(c("B", "A"), c("C", "A")),
                         ncores = 1, quiet = TRUE)

test_that("get proper results from list of comparisons", {
    expect_is(res_2, "list")
    expect_equivalent(names(res_2), c("B_vs_A", "C_vs_A"))
    expect_is(res_2[[1]], "DESeqResults")
    expect_equal(unique(sapply(res_2, nrow)), length(txs))
})

test_that("can use multicore to get results from list of comparisons", {
    cl <- list(c("B", "A"), c("C", "A"))
    expect_equivalent(res_2, getDESeqResults(dds, comparisons.list = cl))
})

test_that("messages when quiet = FALSE", {
    expect_message(getDESeqResults(dds,
                                   contrast.numer = "A",
                                   contrast.denom = "B", quiet = FALSE))
})


test_that("arguments can be passed to DESeq call", {
    res_alt <- getDESeqResults(dds,
                               contrast.numer = "A",
                               contrast.denom = "B",
                               args_DESeq = list(fitType = "mean"),
                               quiet = TRUE)
    expect_is(res_alt, "DESeqResults")
    expect_equivalent(names(res_alt), names(res))
    expect_equivalent(dim(res_alt), dim(res))
    expect_true(all(res$baseMean == res_alt$baseMean))
    expect_false(all(res$pvalue == res_alt$pvalue))
})


test_that("arguments can be passed to results call", {
    res_alt <- getDESeqResults(dds,
                               contrast.numer = "A",
                               contrast.denom = "B",
                               args_results = list(altHypothesis = "greater"),
                               quiet = TRUE)
    expect_is(res_alt, "DESeqResults")
    expect_equivalent(names(res_alt), names(res))
    expect_equivalent(dim(res_alt), dim(res))
    expect_true(all(res$baseMean == res_alt$baseMean))
    expect_false(all(res$pvalue == res_alt$pvalue))
})


