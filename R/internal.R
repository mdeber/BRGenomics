# ========================================================================= #
# Interal helper functions
# ------------------------------------------------------------------------- #

## for output messages, removes newlines and leading spaces from strings

.nicemsg <- function(...) strwrap(..., prefix = " ", initial = "")


## superimposes evenly spaced bins over a vector, and performs FUN on them
## (output length = nbins)

.binVector <- function(x, binsize = NULL, nbins = NULL, FUN = sum) {
    if (missing(nbins))  nbins <- floor(length(x) / binsize)
    if (missing(binsize))  binsize <- floor(length(x) / nbins)

    # bin_breaks <- seq(1, binsize * nbins, binsize)
    # bin_idx <- findInterval(seq_len(binsize * nbins), bin_breaks)
    # aggregate(x[seq_len(binsize*nbins)], by = list(bin_idx), FUN = FUN)[, 2]

    mat <- matrix(x[seq_len(nbins*binsize)], nrow = binsize)
    apply(mat, 2, FUN)
}


## remove undesired chromosomes from GRanges

.tidy_chrom <- function(gr,
                        keep_X = TRUE,
                        keep_Y = TRUE,
                        keep_M = FALSE,
                        keep_nonstandard = FALSE) {

    chrom <- standardChromosomes(gr)

    if (keep_nonstandard) chrom <- seqlevels(gr)
    if (!keep_X)  chrom <- chrom[ chrom != "chrX" ]
    if (!keep_Y)  chrom <- chrom[ chrom != "chrY" ]
    if (!keep_M)  chrom <- chrom[ (chrom != "chrM") & (chrom != "chrMT") ]
    gr <- keepSeqlevels(gr, chrom, pruning.mode = "tidy")
    sortSeqlevels(gr)
}


### convenience function for getting TxDb objects

.getTxDb <- function(genome) {
    if (genome == "hg38") {
        db.txs <-
            TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    } else if (genome == "hg19") {
        db.txs <-
            TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    } else if (genome == "mm10") {
        db.txs <-
            TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
    } else if (genome == "mm9") {
        db.txs <-
            TxDb.Mmusculus.UCSC.mm9.knownGene::TxDb.Mmusculus.UCSC.mm9.knownGene
    } else if (genome == "dm6") {
        db.txs <-
            TxDb.Dmelanogaster.UCSC.dm6.ensGene::TxDb.Dmelanogaster.UCSC.dm6.ensGene
    } else if (genome == "dm3") {
        db.txs <-
            TxDb.Dmelanogaster.UCSC.dm3.ensGene::TxDb.Dmelanogaster.UCSC.dm3.ensGene
    }
}


