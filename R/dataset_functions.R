# ========================================================================= #
# Format GRanges datasets for base-pair resolution analysis
# ------------------------------------------------------------------------- #


#' Make base-pair resolution GRanges object
#'
#' Splits up all ranges in \code{gr} to be each 1 basepair wide. All information
#' is preserved, including all metadata. To wit, \code{length(output.gr) =
#' sum(width(dataset.gr))}.
#'
#' @param gr A disjoint GRanges object.
#'
#' Note that this function doesn't perform any transformation on the metadata
#' in the input; for any ranges of width > 1, the metadata is simply copied
#' to the daughters of that range (whose widths are all equal to 1).
#'
#' This function is intended to work on datasets at single-base resolution. Data
#' of this type is often formatted as a bigWig file, and any data imported from
#' a bigWig file by rtracklayer is suitable for processing. bigWig files will
#' typically use run-length compression on the data signal (the 'score'
#' column), such that when imported by rtracklayer, adjacent bases sharing the
#' same signal will combined into a single range. The base-pair resolution
#' GRanges objects produced by this function remove this compression, resulting
#' in each index (each range) of the GRanges object addressing a single genomic
#' position.
#'
#' To properly use base-pair resolution information, the user should be
#' selecting a single-base from each read, which can be accomplished using
#' GenomicRanges::resize(). Then, single-base coverage can be calculated using
#' getStrandedCoverage.
#'
#' @author Mike DeBerardine
#' @export
makeGRangesBPres <- function(dataset.gr) {

    if (!isDisjoint(dataset.gr)) {
        stop(.nicemsg("Input dataset.gr is not disjoint. See documentation"))
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
#' @param field The name of the field that contains readcounts. If no metadata
#'   field contains readcounts, and each range represents a single read, set to
#'   NULL.
#'
#' @author Mike DeBerardine
#'
#' @export
getStrandedCoverage <- function(gr, field = "score") {

    if (!(field %in% names(mcols(gr)))) {
        msg <- .nicemsg("The given value for 'field' is not found in mcols(gr).
                        If no field contains signal counts for each range,
                        set field = NULL")
        stop(message = msg)
        return(geterrmessage())
    }

    p_gr <- subset(gr, strand == "+")
    m_gr <- subset(gr, strand == "-")
    n_gr <- subset(gr, strand == "*")

    if (is.null(field)) {
        p_cov <- coverage(p_gr)
        m_cov <- coverage(m_gr)
        n_cov <- coverage(n_gr)
    } else {
        p_cov <- coverage(p_gr, weight = mcols(p_gr)[[field]])
        m_cov <- coverage(m_gr, weight = mcols(m_gr)[[field]])
        n_cov <- coverage(n_gr, weight = mcols(n_gr)[[field]])
    }

    p_cov <- GRanges(p_cov, seqinfo = seqinfo(gr))
    m_cov <- GRanges(m_cov, seqinfo = seqinfo(gr))
    n_cov <- GRanges(n_cov, seqinfo = seqinfo(gr))

    # Adding strand info in GRanges() causes error if GRanges object is empty
    strand(p_cov) <- "+"
    strand(m_cov) <- "-"
    strand(n_cov) <- "*"

    cov_gr <- subset(c(p_cov, m_cov, n_cov), score != 0)
    return(sort(cov_gr))
}



#' Randomly subsample reads from GRanges dataset
#'
#' Random subsampling is not performed on ranges, but on reads. Readcounts
#' should be given as a metadata field (usually "score"), and should normally be
#' integers. If normalized readcounts are given, an attempt will be made to
#' infer the normalization factor based on the least-common-multiple of the
#' signal found in the specified field.
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
#' @examples
subsampleGRanges <- function(dataset.gr,
                             n = NULL,
                             prop = NULL,
                             field = "score") {

    if (!xor(missing(n), missing(prop))) {
        msg <- "Must give either 'n' or 'prop' for subsampling, but not both."
        stop(message = msg)
        return(geterrmessage())
    }

    if (is.null(field)) {
        signal_counts <- rep(1, length(dataset.gr))
    } else {
        signal_counts <- mcols(dataset.gr)[[field]]
    }

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

    if (!missing(prop))  n <- floor(prop * sum(signal_counts))

    # avoid expanding GRanges, and sample associated indices
    idx <- rep(seq_along(dataset.gr), times = signal_counts)
    gr_sample <- dataset.gr[sample(idx, n)]

    if (is.null(field)) {
        gr_out <- gr_sample
    } else {
        gr_out <- unique(gr_sample)
        mcols(gr_out)[field] <- countOverlaps(gr_out, gr_sample)
    }

    if (normed_signal) {
        mcols(gr_out)[field] <- mcols(gr_out)[[field]] * lcm
    }

    return(sort(gr_out))
}



#' Merge base-pair resolution GRanges objects
#'
#' Merges 2 or more GRanges objects. For each object, the range widths must all
#' be 1, and the \code{score} metadata column contains coverage information at
#' each site.
#'
#' @param dataset.gr A GRanges object fitting the above description.
#' @param ... One or more additional GRanges objects also fitting the above
#'   description.
#'
#' @return A single GRange object containing all sites of the input objects, and
#'   the sum of all scores at all sites.
#'
#' @author Mike DeBerardine
#' @export
#'
#' @seealso \code{\link{md.import.bigWigs}}
mergeGRangesData <- function(..., field = "score", ncores = detectCores()) {
    data_in <- list(...)

    if (length(data_in) < 2) {
        warning("Running mergeGRangesData on less than 2 GRanges objects")
        if (length(data_in) == 1)
            return(makeGRangesBPres(data_in[[1]]))
        return(GRanges())
    }

    # check if input data are base-pair resolution coverage data
    width_check <- vapply(data_in,
                          function(x) all(width(x) == 1),
                          FUN.VALUE = logical(1))
    if (any(width_check == FALSE)) {
        warning(.nicemsg("One or more inputs are not single-width GRanges,
                         objects. Will coerce them using
                         makeGRangesBPres()..."),
                immediate. = TRUE)
        data_in <- mclapply(data_in,
                            makeGRangesBPres,
                            mc.cores = ncores)
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
