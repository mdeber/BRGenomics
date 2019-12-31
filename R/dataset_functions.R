# ========================================================================= #
# Format GRanges datasets for base-pair resolution analysis
# ------------------------------------------------------------------------- #


#' Make base-pair resolution GRanges object
#'
#' Splits up all ranges in \code{gr} to be each 1 basepair wide. For any range
#' that is split up, all metadata information belonging to that range is
#' inherited by its daughter ranges, and therefore the transformation is
#' non-destructive.
#'
#' @param dataset.gr A disjoint GRanges object
#'
#' @details Note that this function doesn't perform any transformation on the
#'   metadata in the input. This function assumes that for an input GRanges
#'   object, any metadata for each range is equally correct when inherited by
#'   each individual base in that range. In other words, the dataset's "signal"
#'   (usually readcounts) is derived from a single basepair position.
#'
#'   The motivating case for this function is a bigWig file (e.g. one imported
#'   by \code{rtracklayer}), as bigWig files typically use run-length
#'   compression on the data signal (the 'score' column), such that adjacent
#'   bases sharing the same signal are combined into a single range. The
#'   base-pair resolution GRanges objects produced by this function remove this
#'   compression, resulting in each index (each range) of the GRanges object
#'   addressing a single genomic position.
#'
#' @section Generating basepair-resolution GRanges from whole reads: If working
#'   with a GRanges object containing whole reads, one can obtain base-pair
#'   resolution information by using the strand-specific function
#'   \code{\link[GenomicRanges:resize]{GenomicRanges::resize}} to select a
#'   single base from each read: set \code{width = 1} and use the \code{fix}
#'   argument to choose the strand-specific 5' or 3' end. Then, strand-specific
#'   coverage can be calculated using
#'   \code{\link[BRGenomics:getStrandedCoverage]{getStrandedCoverage}}.
#'
#' @author Mike DeBerardine
#' @seealso \code{\link[BRGenomics:getStrandedCoverage]{getStrandedCoverage}},
#'   \code{\link[GenomicRanges:resize]{GenomicRanges::resize()}}
#' @export
#' @importFrom GenomicRanges GRanges GPos mcols mcols<- isDisjoint findOverlaps
#' @examples
#' #--------------------------------------------------#
#' # Make a bigWig file single width
#' #--------------------------------------------------#
#'
#' # get local address for an included bigWig file
#' bw_file <- system.file("extdata", "PROseq_dm6_chr4_plus.bw",
#'                        package = "BRGenomics")
#'
#' # BRGenomics::import_bigWig automatically applies makeGRangesBRG;
#' # therefore will import using rtracklayer
#' bw <- rtracklayer::import.bw(bw_file)
#' strand(bw) <- "+"
#'
#' range(width(bw))
#' length(bw)
#'
#' # make basepair-resolution (single-width)
#' gr <- makeGRangesBRG(bw)
#'
#' range(width(gr))
#' length(gr)
#' length(gr) == sum(width(bw))
#' sum(score(gr)) == sum(score(bw) * width(bw))
#'
#' #--------------------------------------------------#
#' # Reverse using getStrandedCoverage
#' #--------------------------------------------------#
#' # -> for more examples, see getStrandedCoverage
#'
#' undo <- getStrandedCoverage(gr)
#'
#' range(width(undo))
#' length(undo) == length(bw)
#' all(score(undo) == score(bw))
makeGRangesBRG <- function(dataset.gr) {

    if (!isDisjoint(dataset.gr)) {
        stop("Input dataset.gr is not disjoint. See documentation")
        return(geterrmessage())
    }

    # Make all widths = 1
    gr_bp <- GRanges(GPos(dataset.gr))

    # Add back all metadata
    hits <- findOverlaps(dataset.gr, gr_bp) # find corresponding indices
    mcols(gr_bp) <- mcols(dataset.gr)[hits@from, ]
    names(mcols(gr_bp)) <- names(mcols(dataset.gr))

    return(sort(gr_bp))
}



#' Get strand-specific coverage
#'
#' Computes strand-specific coverage signal, and returns a GRanges object with
#' signal in the "score" metadata column. Function also works for
#' non-strand-specific data. Note that output is not automatically converted
#' into a "basepair-resolution" GRanges object.
#'
#' @param dataset.gr A GRanges object either containing ranges for each read, or
#'   one in which readcounts for individual ranges are contained in metadata
#'   (typically in the "score" field).
#' @param field The name of the metadata field that contains readcounts. If no
#'   metadata field contains readcounts, and each range represents a single
#'   read, set to NULL.
#'
#' @author Mike DeBerardine
#' @seealso \code{\link[BRGenomics:makeGRangesBRG]{makeGRangesBRG}},
#'   \code{\link[GenomicRanges:coverage]{GenomicRanges::coverage}}
#' @export
#' @importFrom GenomicRanges mcols
#' @examples
#' #--------------------------------------------------#
#' # Using included full-read data
#' #--------------------------------------------------#
#' # -> whole-read coverage sacrifices meaningful readcount
#' #    information, but can be useful for visualization,
#' #    e.g. for looking at RNA-seq data in a genome browser
#'
#' data("PROseq_paired")
#'
#' PROseq_paired[1:6]
#'
#' getStrandedCoverage(PROseq_paired)[1:6]
#'
#' #--------------------------------------------------#
#' # Getting coverage from single bases of single reads
#' #--------------------------------------------------#
#'
#' # included PROseq data is already single-base coverage
#' data("PROseq")
#' range(width(PROseq))
#'
#' # undo coverage for the first 100 positions
#' ps <- PROseq[1:100]
#' ps_reads <- rep(ps, times = ps$score)
#' mcols(ps_reads) <- NULL
#'
#' ps_reads[1:6]
#'
#' # re-create coverage
#' getStrandedCoverage(ps_reads, field = NULL)[1:6]
#'
#' #--------------------------------------------------#
#' # Reversing makeGRangesBRG
#' #--------------------------------------------------#
#' # -> getStrandedCoverage doesn't return single-width
#' #    GRanges, which is useful because getting coverage
#' #    will merge adjacent bases with equivalent scores
#'
#' # included PROseq data is already single-width
#' range(width(PROseq))
#' isDisjoint(PROseq)
#'
#' ps_cov <- getStrandedCoverage(PROseq)
#'
#' range(width(ps_cov))
#' sum(score(PROseq)) == sum(score(ps_cov) * width(ps_cov))
#'
#' # -> Look specifically at ranges that could be combined
#' neighbors <- c(shift(PROseq, 1), shift(PROseq, -1))
#' hits <- findOverlaps(PROseq, neighbors)
#' idx <- unique(from(hits)) # indices for PROseq with neighbor
#'
#' PROseq[idx]
#'
#' getStrandedCoverage(PROseq[idx])
getStrandedCoverage <- function(dataset.gr, field = "score") {

    if (!is.null(field) && !(field %in% names(mcols(dataset.gr)))) {
        msg <- .nicemsg("The given value for 'field' is not found in
                        mcols(dataset.gr). If no field contains signal counts
                        for each range, set field = NULL")
        stop(message = msg)
        return(geterrmessage())
    }

    p_cov <- .get_stranded_cov(dataset.gr, "+", field)
    m_cov <- .get_stranded_cov(dataset.gr, "-", field)
    n_cov <- .get_stranded_cov(dataset.gr, "*", field)

    cov_gr <- c(p_cov, m_cov, n_cov)
    return(sort(cov_gr))
}

#' @importFrom GenomicRanges strand coverage mcols GRanges seqinfo strand<-
#'   score
.get_stranded_cov <- function(dataset.gr, strand_i, field) {
    gr <- subset(dataset.gr, strand == strand_i)
    if (is.null(field)) {
        cv <- coverage(gr)
    } else {
        cv <- coverage(gr, weight = mcols(gr)[[field]])
    }

    cv_gr <- GRanges(cv, seqinfo = seqinfo(dataset.gr))
    strand(cv_gr) <- strand_i
    subset(cv_gr, score != 0)
}


#' Randomly subsample reads from GRanges dataset
#'
#' Random subsampling is not performed on ranges, but on reads. Readcounts
#' should be given as a metadata field (usually "score"), and should normally be
#' integers. If normalized readcounts are given, an attempt will be made to
#' infer the normalization factor based on the least-common-multiple of the
#' signal found in the specified field. This function can also subsample ranges
#' directly if \code{field = NULL}, but the \code{sample} function can be used
#' in this scenario.
#'
#' @param dataset.gr A GRanges object in which signal (e.g. readcounts) are
#'   contained within metadata.
#' @param n Number of reads to subsample. Either \code{n} or \code{prop} can be
#'   given.
#' @param prop Proportion of total signal to subsample.
#' @param field The metadata field of \code{dataset.gr} that contains readcounts
#'   for reach position. If each range represents a single read, set \code{field
#'   = NULL}
#'
#' @author Mike DeBerardine
#' @export
#' @importFrom GenomicRanges mcols mcols<- countOverlaps
#' @examples
#' data("PROseq") # load included PROseq data
#'
#' #--------------------------------------------------#
#' # sample 10% of the reads of a GRanges with signal coverage
#' #--------------------------------------------------#
#'
#' ps_sample <- subsampleGRanges(PROseq, prop = 0.1)
#'
#' # cannot predict number of ranges (positions) that will be sampled
#' length(PROseq)
#' length(ps_sample)
#'
#' # 1/10th the score is sampled
#' sum(score(PROseq))
#' sum(score(ps_sample))
#'
#' #--------------------------------------------------#
#' # Sample 10% of ranges (e.g. if each range represents one read)
#' #--------------------------------------------------#
#'
#' ps_sample <- subsampleGRanges(PROseq, prop = 0.1, field = NULL)
#'
#' length(PROseq)
#' length(ps_sample)
#'
#' # Alternatively
#' ps_sample <- sample(PROseq, 0.1 * length(PROseq))
#' length(ps_sample)
subsampleGRanges <- function(dataset.gr,
                             n = NULL,
                             prop = NULL,
                             field = "score") {

    if (!xor(is.null(n), is.null(prop))) {
        msg <- "Must give either 'n' or 'prop' for subsampling, but not both."
        stop(message = msg)
        return(geterrmessage())
    }

    if (is.null(field)) {
        if (is.null(n)) n <- floor(prop * length(dataset.gr))
        return(sample(dataset.gr, n))
    }

    signal_counts <- mcols(dataset.gr)[[field]]
    if (is.null(n))  n <- round(prop * sum(signal_counts))

    if (all( round(signal_counts, 3) %% 1 == 0 )) {
        normed_signal <- FALSE
    } else {
        normed_signal <- TRUE
        lcm <- min(signal_counts)
        unnorm_signal <- signal_counts / lcm
        if (!all( round(unnorm_signal, 3) %% 1 == 0 )) {
            stop(message = .nicemsg("Signal given in 'field' are not whole
                                    numbers, and unable to infer a normalization
                                    factor."))
            return(geterrmessage())
        } else {
            warning(.nicemsg("Signal given in 'field' are not whole numbers. A
                             normalization factor was inferred based on the
                             least common multiple."))
            signal_counts <- round(unnorm_signal)
        }
    }

    # avoid expanding GRanges, and sample associated indices
    idx <- rep(seq_along(dataset.gr), times = signal_counts)
    gr_sample <- dataset.gr[sample(idx, n)]
    gr_out <- unique(gr_sample)
    mcols(gr_out)[field] <- countOverlaps(gr_out, gr_sample)

    if (normed_signal) mcols(gr_out)[field] <- mcols(gr_out)[[field]] * lcm
    return(sort(gr_out))
}



#' Merge base-pair resolution GRanges objects
#'
#' Merges 2 or more GRanges objects. For each object, the range widths must all
#' be 1, and the \code{score} metadata column contains coverage information at
#' each site. This function returns a single GRange object containing all sites
#' of the input objects, and the sum of all scores at all sites.
#'
#' @param ... Any number of GRanges objects in which signal (e.g. readcounts)
#'   are contained within metadata.
#' @param field One or more metadata fields to be combined, typically the
#'   "score" field. Fields typically contain coverage information.
#' @param ncores More than one core can be used to coerce non-single-width
#'   GRanges objects using \code{makeGRangesBRG}.
#' @author Mike DeBerardine
#' @seealso \code{\link[BRGenomics:makeGRangesBRG]{makeGRangesBRG}}
#' @export
#' @importFrom parallel detectCores mclapply
#' @importFrom GenomicRanges GRanges width findOverlaps mcols mcols<-
#' @examples
#' data("PROseq") # load included PROseq data
#'
#' #--------------------------------------------------#
#' # divide & recombine PROseq (no overlapping positions)
#' #--------------------------------------------------#
#'
#' thirds <- floor( (1:3)/3 * length(PROseq) )
#' ps_1 <- PROseq[1:thirds[1]]
#' ps_2 <- PROseq[(thirds[1]+1):thirds[2]]
#' ps_3 <- PROseq[(thirds[2]+1):thirds[3]]
#'
#' # re-merge
#' length(PROseq)
#' length(ps_1)
#' length(mergeGRangesData(ps_1, ps_2))
#' length(mergeGRangesData(ps_1, ps_2, ps_3))
#'
#' #--------------------------------------------------#
#' # combine PRO-seq with overlapping positions
#' #--------------------------------------------------#
#'
#' gr1 <- PROseq[10:13]
#' gr2 <- PROseq[12:15]
#'
#' PROseq[10:15]
#'
#' mergeGRangesData(gr1, gr2)
mergeGRangesData <- function(..., field = "score", ncores = detectCores()) {
    data_in <- list(...)

    if (length(data_in) < 2) {
        warning("Running mergeGRangesData on less than 2 GRanges objects")
        if (length(data_in) == 1)  return(makeGRangesBRG(data_in[[1]]))
        return(GRanges())
    }

    # check if input data are base-pair resolution coverage data
    width_check <- vapply(data_in, function(x) all(width(x) == 1),
                          FUN.VALUE = logical(1))
    if (any(width_check == FALSE)) {
        warning(.nicemsg("One or more inputs are not single-width GRanges
                         objects. Will coerce them using makeGRangesBRG()..."),
                immediate. = TRUE)
        data_in <- mclapply(data_in, makeGRangesBRG, mc.cores = ncores)
    }

    # merge input GRanges
    gr <- data_in[[1]]

    for (i in data_in[-1]) {
        # add signal from overlapping sites
        hits <- findOverlaps(i, gr)
        signal_i <- mcols(i)[[field]][hits@from]
        signal_gr <- mcols(gr)[[field]][hits@to]
        mcols(gr)[[field]][hits@to] <- signal_gr + signal_i

        # add ranges from non-overlapping sites
        gr <- c(gr, i[!(seq_along(i) %in% hits@from)])
        gr <- sort(gr)
    }
    return(gr)
}
