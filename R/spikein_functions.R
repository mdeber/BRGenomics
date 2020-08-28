
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
#' @param expand_ranges Logical indicating if ranges in \code{dataset.gr} should
#'   be treated as descriptions of single molecules (\code{FALSE}), or if ranges
#'   should be treated as representing multiple adjacent positions with the same
#'   signal (\code{TRUE}). See \code{\link[BRGenomics:getCountsByRegions]{
#'   getCountsByRegions}}.
#' @param ncores The number of cores to use for computations.
#'
#' @return A dataframe containing total readcounts, experimental (non-spike-in)
#'   readcounts, and spike-in readcounts.
#' @author Mike DeBerardine
#'
#' @importFrom methods is
#' @import GenomicRanges
#' @export
#'
#' @examples
#' #--------------------------------------------------#
#' # Make list of dummy GRanges
#' #--------------------------------------------------#
#'
#' gr1_rep1 <- GRanges(seqnames = c("chr1", "chr2", "spikechr1", "spikechr2"),
#'                     ranges = IRanges(start = 1:4, width = 1),
#'                     strand = "+")
#' gr2_rep2 <- gr2_rep1 <- gr1_rep2 <- gr1_rep1
#'
#' # set readcounts
#' score(gr1_rep1) <- c(1, 1, 1, 1) # 2 exp + 2 spike = 4 total
#' score(gr2_rep1) <- c(2, 2, 1, 1) # 4 exp + 2 spike = 6 total
#' score(gr1_rep2) <- c(1, 1, 2, 1) # 2 exp + 3 spike = 5 total
#' score(gr2_rep2) <- c(4, 4, 2, 2) # 8 exp + 4 spike = 12 total
#'
#' grl <- list(gr1_rep1, gr2_rep1,
#'             gr1_rep2, gr2_rep2)
#'
#' names(grl) <- c("gr1_rep1", "gr2_rep1",
#'                 "gr1_rep2", "gr2_rep2")
#'
#' grl
#'
#' #--------------------------------------------------#
#' # Count spike-in reads
#' #--------------------------------------------------#
#'
#' # by giving names of all spike-in chromosomes
#' getSpikeInCounts(grl, si_names = c("spikechr1", "spikechr2"), ncores = 1)
#'
#' # or by matching the string/regular expression "spike" in chromosome names
#' getSpikeInCounts(grl, si_pattern = "spike", ncores = 1)
#'
#' #--------------------------------------------------#
#' # Filter out spike-in reads
#' #--------------------------------------------------#
#'
#' removeSpikeInReads(grl, si_pattern = "spike", ncores = 1)
#'
#' #--------------------------------------------------#
#' # Return spike-in reads
#' #--------------------------------------------------#
#'
#' getSpikeInReads(grl, si_pattern = "spike", ncores = 1)
getSpikeInCounts <- function(dataset.gr, si_pattern = NULL, si_names = NULL,
                             field = "score", sample_names = NULL,
                             expand_ranges = FALSE,
                             ncores = getOption("mc.cores", 2L)) {

    if (!is.list(dataset.gr) && !is(dataset.gr, "GRangesList")) {
        name_in <- deparse(substitute(dataset.gr))
        dataset.gr <- list(dataset.gr)
        names(dataset.gr) <- name_in
    }

    if (!is.null(sample_names))
        names(dataset.gr) <- sample_names

    spike_chrom <- .get_spike_chrom(dataset.gr, si_pattern, si_names, ncores)
    FUN <- if (expand_ranges) .get_spikecounts_expand else .get_spikecounts
    FUN(dataset.gr, spike_chrom, field, ncores)
}

#' @importFrom methods is
#' @importFrom parallel mclapply
#' @importFrom GenomicRanges seqinfo
.get_spike_chrom <- function(X, si_pattern, si_names, ncores = 1L) {

    if (is.list(X) || is(X, "GRangesList")) {
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

#' @importFrom parallel mclapply mcMap
#' @importFrom GenomicRanges seqnames mcols
.get_spikecounts_expand <- function(dataset.list, spike_chrom, field, ncores) {

    snames <- names(dataset.list)

    if (is.null(field))
        stop(.nicemsg("Cannot use field = NULL when expand_ranges = TRUE.
                      expand_ranges is for collapsed coverage objects, which
                      always have associated signal for each position"))

    cl <- mcMap(function(x, field) {
        si <- as.logical(seqnames(x) %in% spike_chrom)
        data.frame(total_reads = sum( mcols(x)[[field]] * width(x)),
                   exp_reads = sum( mcols(x)[[field]][!si] * width(x)[!si]),
                   spike_reads = sum( mcols(x)[[field]][si] * width(x)[si]))
    }, dataset.list, field, mc.cores = ncores)

    names(cl) <- snames
    .dfList2df(cl)
}

#' @importFrom methods is
#' @importFrom GenomeInfoDb dropSeqlevels
#' @importFrom parallel mclapply
#' @rdname getSpikeInCounts
#' @export
removeSpikeInReads <- function(dataset.gr, si_pattern = NULL, si_names = NULL,
                               field = "score",
                               ncores = getOption("mc.cores", 2L)) {
    spike_chrom <- .get_spike_chrom(dataset.gr, si_pattern, si_names, ncores)

    if (is.list(dataset.gr) || is(dataset.gr, "GRangesList")) {
        mclapply(dataset.gr, dropSeqlevels, spike_chrom, pruning.mode = "tidy")
    } else {
        dropSeqlevels(dataset.gr, spike_chrom, pruning.mode = "tidy")
    }
}

#' @importFrom methods is
#' @importFrom GenomeInfoDb keepSeqlevels
#' @importFrom parallel mclapply
#' @rdname getSpikeInCounts
#' @export
getSpikeInReads <- function(dataset.gr, si_pattern = NULL, si_names = NULL,
                            field = "score",
                            ncores = getOption("mc.cores", 2L)) {
    spike_chrom <- .get_spike_chrom(dataset.gr, si_pattern, si_names, ncores)

    if (is.list(dataset.gr) || is(dataset.gr, "GRangesList")) {
        mclapply(dataset.gr, keepSeqlevels, spike_chrom, pruning.mode = "tidy")
    } else {
        keepSeqlevels(dataset.gr, spike_chrom, pruning.mode = "tidy")
    }
}
