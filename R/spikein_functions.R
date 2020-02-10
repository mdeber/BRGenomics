
#' Filtering and counting spike-in reads
#'
#' @param dataset.gr A GRanges object or a list of GRanges objects.
#' @param si_pattern A regular expression that matches spike-in chromosomes. Can
#'   be used in addition to, or as an alternative to \code{si_names}.
#' @param si_names A character vector giving the names of the spike-in
#'   chromosomes. Can be used in addition to, or as an alternative to
#'   \code{si_pattern}.
#' @param field The metadata field in \code{dataset.gr} that contains
#'   readcounts. If each range is an individual read, set \code{field = NULL}.
#' @param sample_names An optional character vector used to rename the datasets
#'   in \code{dataset.gr}
#' @param ncores The number of cores to use for computations.
#'
#' @return A dataframe containing total readcounts, experimental (non-spike-in)
#'   readcounts, and spike-in readcounts.
#' @author Mike DeBerardine
#'
#' @import GenomicRanges
#' @export
#'
#' @examples
getSpikeInCounts <- function(dataset.gr, si_pattern = NULL, si_names = NULL,
                             field = "score", sample_names = NULL,
                             ncores = detectCores()) {

    if (!is.list(dataset.gr)) {
        name_in <- deparse(substitute(dataset.gr))
        dataset.gr <- list(dataset.gr)
        names(dataset.gr) <- name_in
    }

    if (!is.null(sample_names))  names(dataset.gr) <- sample_names

    spike_chrom <- .get_spike_chrom(dataset.gr, si_pattern, si_names, ncores)
    .get_spikecounts(dataset.gr, spike_chrom, field, ncores)
}

#' @importFrom parallel mclapply
#' @importFrom GenomicRanges seqinfo
.get_spike_chrom <- function(X, si_pattern, si_names, ncores = 1) {

    if (is.list(X)) {
        chrom <- mclapply(X, function(x) names(seqinfo(x)), mc.cores = ncores)
        chrom <- unique(unlist(chrom))
    } else {
        chrom <- names(seqinfo(X))
    }

    spike_chrom <- c()

    if (!is.null(si_pattern))
        spike_chrom <- grep(si_pattern, chrom, value = TRUE)

    if (!is.null(si_names))
        spike_chrom <- c(spike_chrom, si_names)

    unique(spike_chrom)
}

#' @importFrom parallel mclapply mcMap
#' @importFrom GenomicRanges seqnames mcols
.get_spikecounts <- function(dataset.list, spike_chrom, field, ncores) {

    snames <- names(dataset.list)

    if (is.null(field)) {
        cl <- mclapply(dataset.list, function(x) {
            si <- seqnames(x) %in% spike_chrom
            data.frame(total_reads = length(x),
                       exp_reads = sum(!si),
                       spike_reads = sum(si))
        }, mc.cores = ncores)

    } else {
        cl <- mcMap(function(x, field) {
            si <- seqnames(x) %in% spike_chrom
            data.frame(total_reads = sum( mcols(x)[[field]] ),
                       exp_reads = sum( mcols(x[!si])[[field]] ),
                       spike_reads = sum( mcols(x[si])[[field]] ))
        }, dataset.list, field, mc.cores = ncores)
    }

    names(cl) <- snames
    .dfList2df(cl)
}


#' @importFrom GenomeInfoDb dropSeqlevels
#' @importFrom parallel mclapply
#' @rdname getSpikeInCounts
#' @export
removeSpikeInReads <- function(dataset.gr, si_pattern = NULL, si_names = NULL,
                               field = "score", ncores = detectCores()) {
    spike_chrom <- .get_spike_chrom(dataset.gr, si_pattern, si_names, ncores)

    if (is.list(dataset.gr)) {
        mclapply(dataset.gr, dropSeqlevels, spike_chrom, pruning.mode = "tidy")
    } else {
        dropSeqlevels(dataset.gr, spike_chrom, pruning.mode = "tidy")
    }
}


#' @importFrom GenomeInfoDb keepSeqlevels
#' @importFrom parallel mclapply
#' @rdname getSpikeInCounts
#' @export
getSpikeInReads <- function(dataset.gr, si_pattern = NULL, si_names = NULL,
                            field = "score", ncores = detectCores()) {
    spike_chrom <- .get_spike_chrom(dataset.gr, si_pattern, si_names, ncores)

    if (is.list(dataset.gr)) {
        mclapply(dataset.gr, keepSeqlevels, spike_chrom, pruning.mode = "tidy")
    } else {
        keepSeqlevels(dataset.gr, spike_chrom, pruning.mode = "tidy")
    }
}



