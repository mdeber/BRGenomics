# ========================================================================= #
# Format GRanges datasets for base-pair resolution analysis
# ------------------------------------------------------------------------- #


#' Make base-pair resolution GRanges object
#'
#' Splits up all ranges in \code{dataset.gr} to be each 1 basepair wide. For any
#' range that is split up, all metadata information belonging to that range is
#' inherited by its daughter ranges, and therefore the transformation is
#' non-destructive.
#'
#' @param dataset.gr A disjoint GRanges object
#'
#' @return A GRanges object for which \code{length(output) ==
#'   sum(width(dataset.gr))}, and for which \code{all(width(output) == 1)}.
#'
#' @details Note that this function doesn't perform any transformation on the
#'   metadata in the input. This function assumes that for an input GRanges
#'   object, any metadata for each range is equally correct when inherited by
#'   each individual base in that range. In other words, the dataset's "signal"
#'   (usually readcounts) fundamentally belongs to a single basepair position.
#'
#' @section Motivation: The motivating case for this function is a bigWig file
#'   (e.g. one imported by \code{rtracklayer}), as bigWig files typically use
#'   run-length compression on the data signal (the 'score' column), such that
#'   adjacent bases sharing the same signal are combined into a single range. As
#'   basepair-resolution genomic data is typically sparse, this compression has
#'   a minimal impact on memory usage, and removing it greatly enhances data
#'   handling as each index (each range) of the GRanges object corresponds to a
#'   single genomic position.
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
#' @section On the use of GRanges instead of GPos: The
#'   \code{\link[GenomicRanges:GPos]{GPos}} class is a more suitable container
#'   for data of this type, as the GPos class is specific to 1-bp-wide ranges.
#'   However, in early testing, we encountered some kind of compatibility
#'   limitations with the newer GPos class, and have not re-tested it since. If
#'   you have feedback on switching to this class, please contact the author.
#'   Users can readily coerce a basepair-resolution GRanges object to a GPos
#'   object via \code{gp <- GPos(gr, score = score(gr))}.
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
#' Computes strand-specific coverage signal, and returns a GRanges object.
#' Function also works for non-strand-specific data.
#'
#' @param dataset.gr A GRanges object either containing ranges for each read, or
#'   one in which readcounts for individual ranges are contained in metadata
#'   (typically in the "score" field).
#' @param field The name of the metadata field that contains readcounts. If no
#'   metadata field contains readcounts, and each range represents a single
#'   read, set to NULL.
#'
#' @return A GRanges object with signal in the "score" metadata column. Note
#'   that the output is \emph{not} automatically converted into a
#'   \code{\link[BRGenomics:makeGRangesBRG]{"basepair-resolution"}} GRanges
#'   object.
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
#' should be given as a metadata field (usually "score"). This function can also
#' subsample ranges directly if \code{field = NULL}, but the \code{sample}
#' function can be used in this scenario.
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
#' @return A GRanges object identical in format to \code{dataset.gr}, but
#'   containing a random subset of its data. If \code{field != NULL}, the length
#'   of the output cannot be known \emph{a priori}, but the sum of its score
#'   can.
#'
#' @section Use with normalized readcounts: If the metadata field contains
#'   normalized readcounts, an attempt will be made to infer the normalization
#'   factor based on the least-common-multiple of the signal found in the
#'   specified field.
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



#' Merge basepair-resolution GRanges objects
#'
#' Merges 2 or more basepair-resolution (single-width) GRanges objects by
#' combining all of their ranges and associated signal (e.g. readcounts). If
#' \code{multiplex = TRUE}, the input datasets are reversibly combined into a
#' multiplexed GRanges containing a field for each input dataset.
#'
#' @param ... Any number of GRanges objects in which signal (e.g. readcounts)
#'   are contained within metadata. GRanges not single-width will be coerced
#'   using \code{\link[BRGenomics:makeGRangesBRG]{makeGRangesBRG}}. Lists of
#'   GRanges can also be passed, but they must be named lists if
#'   \code{multiplex = TRUE}. Multiple lists can be passed, but if any inputs
#'   are lists, then all inputs must be lists.
#' @param field One or more \emph{input} metadata fields to be combined,
#'   typically the "score" field. Fields typically contain coverage information.
#'   If only a single field is given (i.e. all input GRanges use the same
#'   field), that same field will be used for the output. Otherwise, the
#'   \code{score} metadata field will be used by default. The output metadata
#'   fields are different if \code{multiplex} is enabled.
#' @param multiplex When set to \code{FALSE} (the default), input GRanges are
#'   merged irreversibly into a single new GRange, effectively combining the
#'   reads from different experiments. When \code{multiplex = TRUE}, the input
#'   GRanges data are reversibly combined into a multiplexed GRanges object,
#'   such that each input GRanges object has its own metadata field in the
#'   output.
#' @param ncores Number of cores to use for computations.
#'
#' @return A disjoint, basepair-resolution (single-width) GRanges object
#'   comprised of all ranges found in the input GRanges objects.
#'
#'   If \code{multiplex = FALSE}, single fields from each input are combined
#'   into a single field in the output, the total signal of which is the sum of
#'   all input GRanges.
#'
#'   If \code{multiplex = TRUE}, each field of the output corresponds to an
#'   input GRanges object.
#'
#' @section Subsetting a multiplexed GRanges object: If \code{multiplex = TRUE},
#'   the datasets are only combined into a single object, but the data
#'   themselves are not combined. To subset \code{field_i}, corresponding to
#'   input \code{dataset_i}:
#'
#'   \code{multi.gr <- mergeGRangesData(gr1, gr2, multiplex = TRUE)}
#'   \code{subset(multi.gr, gr1 != 0, select = gr1)} # select gr1
#'
#' @author Mike DeBerardine
#' @seealso \code{\link[BRGenomics:makeGRangesBRG]{makeGRangesBRG}}
#' @export
#' @importFrom parallel detectCores mclapply mcmapply
#' @importFrom GenomicRanges mcols mcols<- sort
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
#' length(mergeGRangesData(ps_1, ps_2, ncores = 2))
#' length(mergeGRangesData(ps_1, ps_2, ps_3, ncores = 2))
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
#' mergeGRangesData(gr1, gr2, ncores = 2)
#'
#' #--------------------------------------------------#
#' # multiplex separate PRO-seq experiments
#' #--------------------------------------------------#
#'
#' multi.gr <- mergeGRangesData(gr1, gr2, multiplex = TRUE, ncores = 2)
#' multi.gr
#'
#' #--------------------------------------------------#
#' # subset a multiplexed GRanges object
#' #--------------------------------------------------#
#'
#' subset(multi.gr, gr1 > 0)
#'
#' subset(multi.gr, gr1 > 0, select = gr1)
mergeGRangesData <- function(..., field = "score", multiplex = FALSE,
                             ncores = detectCores()) {
    data_in <- list(...)
    if (any(vapply(data_in, is.list, logical(1))))  data_in <- unlist(data_in)

    if (is.null(names(data_in))) {
        exclude <- c("field", "multiplex", "ncores")
        in.names <- as.list(match.call())[-1]
        names(data_in) <- in.names[!names(in.names) %in% exclude]
    }

    # check length of fields vs. data_in
    field <- .check_merge_fields(data_in, field)

    # data must be *sorted*, base-pair resolution coverage data
    data_in <- .check_sorted_BRG(data_in, ncores)

    # merge ranges
    gr <- do.call(c, c(data_in, use.names = FALSE))
    mcols(gr) <- NULL
    gr <- unique(sort(gr))

    # Fastest to keep these steps separated (esp. for large datasets)
    idx <- mclapply(data_in, function(x) which(gr %in% x), mc.cores = ncores)

    counts <- mcmapply(function(dat, idx, field) {
        out <- rep(0L, length(gr))
        out[idx] <- mcols(dat)[[field]]
        out
    }, data_in, idx, field, mc.cores = ncores, SIMPLIFY = TRUE)

    if (multiplex) {
        mcols(gr)[names(data_in)] <- counts
    } else {
        field <- ifelse(length(unique(field)) == 1, field[1], "score")
        # use of apply will maintain integers; rowSums would make numeric
        mcols(gr)[field] <- apply(counts, 1, sum)
    }

    return(gr)
}

.check_merge_fields <- function(data_in, field) {

    if (length(field) > 1) {
        if (length(field) != length(data_in)) {
            stop("given fields not equal to number of datasets")
            return(geterrmessage())
        }
    } else {
        field <- rep(field, length(data_in))
    }

    return(field)
}

#' @importFrom GenomicRanges width sort
.check_sorted_BRG <- function(data_in, ncores) {

    widths_ok <- vapply(data_in, function(x) all(width(x) == 1),
                        FUN.VALUE = logical(1))

    if (all(widths_ok)) {
        # ensure all data is sorted
        data_in <- mclapply(data_in, sort, mc.cores = ncores)

    } else {
        warning(.nicemsg("One or more inputs are not single-width GRanges
                         objects. Will coerce them using makeGRangesBRG()..."),
                immediate. = TRUE)
        # (makeGRangesBRG sorts its output)
        data_in <- mclapply(data_in, makeGRangesBRG, mc.cores = ncores)
    }

    return(data_in)
}
