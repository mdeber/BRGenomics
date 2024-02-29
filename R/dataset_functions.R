# ========================================================================= #
# Format GRanges datasets for base-pair resolution analysis
# ------------------------------------------------------------------------- #

#' Constructing and checking for base-pair resolution GRanges objects
#'
#' \code{makeGRangesBRG} splits up all ranges in \code{dataset.gr} to be each 1
#' basepair wide. For any range that is split up, all metadata information
#' belonging to that range is inherited by its daughter ranges, and therefore
#' the transformation is non-destructive. \code{isBRG} checks whether an object
#' is a basepair resolution GRanges object.
#'
#' @param dataset.gr A disjoint GRanges object, or a list of such objects.
#' @param ncores If \code{dataset.gr} is a list, the number of cores to use for
#'   computations.
#'
#' @return \code{makeGRangesBRG} returns a GRanges object for which
#'   \code{length(output) == sum(width(dataset.gr))}, and for which
#'   \code{all(width(output) == 1)}.
#'
#'   \code{isBRG(x)} returns \code{TRUE} if \code{x} is a GRanges object with
#'   the above characteristics.
#'
#' @details Note that \code{makeGRangesBRG} doesn't perform any transformation
#'   on the metadata in the input. This function assumes that for an input
#'   GRanges object, any metadata for each range is equally correct when
#'   inherited by each individual base in that range. In other words, the
#'   dataset's "signal" (usually readcounts) fundamentally belongs to a single
#'   basepair position.
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
#'   \code{\link[GenomicRanges:intra-range-methods]{GenomicRanges::resize}} to
#'   select a single base from each read: set \code{width = 1} and use the
#'   \code{fix} argument to choose the strand-specific 5' or 3' end. Then,
#'   strand-specific coverage can be calculated using
#'   \code{\link[BRGenomics:getStrandedCoverage]{getStrandedCoverage}}.
#'
#' @section On the use of GRanges instead of GPos: The
#'   \code{\link[GenomicRanges:GPos-class]{GPos}} class is a more suitable
#'   container for data of this type, as the GPos class is specific to 1-bp-wide
#'   ranges. However, in early testing, we encountered some kind of
#'   compatibility limitations with the newer GPos class, and have not re-tested
#'   it since. If you have feedback on switching to this class, please contact
#'   the author. Users can readily coerce a basepair-resolution GRanges object
#'   to a GPos object via \code{gp <- GPos(gr, score = score(gr))}.
#'
#' @author Mike DeBerardine
#' @seealso \code{\link[BRGenomics:getStrandedCoverage]{getStrandedCoverage}},
#'   \code{\link[GenomicRanges:intra-range-methods]{GenomicRanges::resize()}}
#' @export
#' @importFrom methods is
#' @importFrom GenomicRanges width GPos mcols mcols<- isDisjoint findOverlaps
#'   GRanges sort
#' @importFrom S4Vectors from
#' @examples
#' if (.Platform$OS.type == "unix") {
#'
#'     #--------------------------------------------------#
#'     # Make a bigWig file single width
#'     #--------------------------------------------------#
#'
#'     # get local address for an included bigWig file
#'     bw_file <- system.file("extdata", "PROseq_dm6_chr4_plus.bw",
#'                            package = "BRGenomics")
#'
#'     # BRGenomics::import_bigWig automatically applies makeGRangesBRG;
#'     # therefore will import using rtracklayer
#'     bw <- rtracklayer::import.bw(bw_file)
#'     strand(bw) <- "+"
#'
#'     range(width(bw))
#'     length(bw)
#'
#'     # make basepair-resolution (single-width)
#'     gr <- makeGRangesBRG(bw)
#'
#'     isBRG(gr)
#'     range(width(gr))
#'     length(gr)
#'     length(gr) == sum(width(bw))
#'     sum(score(gr)) == sum(score(bw) * width(bw))
#'
#'     #--------------------------------------------------#
#'     # Reverse using getStrandedCoverage
#'     #--------------------------------------------------#
#'     # -> for more examples, see getStrandedCoverage
#'
#'     undo <- getStrandedCoverage(gr, ncores = 1)
#'
#'     isBRG(undo)
#'     range(width(undo))
#'     length(undo) == length(bw)
#'     all(score(undo) == score(bw))
#'
#' }
makeGRangesBRG <- function(dataset.gr, ncores = getOption("mc.cores", 2L)) {

    if (is.list(dataset.gr) || is(dataset.gr, "GRangesList"))
        return(mclapply(dataset.gr, makeGRangesBRG, mc.cores = ncores))

    if (!isDisjoint(dataset.gr))
        stop("Input dataset.gr is not disjoint. See documentation")

    # separate wide granges
    is_wide <- width(dataset.gr) > 1L
    gr_wide <- dataset.gr[is_wide]

    # make width 1, reverse map, and add metadata
    gp <- GPos(gr_wide)
    hits <- findOverlaps(gr_wide, gp)
    mcols(gp) <- mcols(gr_wide)[from(hits), , drop = FALSE]

    # combine and sort
    sort(c(dataset.gr[!is_wide], GRanges(gp)))
}

#' @rdname makeGRangesBRG
#' @param x Object to be tested.
#' @importFrom methods is
#' @importFrom GenomicRanges width isDisjoint
#' @export
isBRG <- function(x) {
    is(x, "GRanges") && all(width(x) == 1L) && isDisjoint(x)
}


#' Get strand-specific coverage
#'
#' Computes strand-specific coverage signal, and returns a GRanges object.
#' Function also works for non-strand-specific data.
#'
#' @param dataset.gr A GRanges object either containing ranges for each read, or
#'   one in which readcounts for individual ranges are contained in metadata
#'   (typically in the "score" field). \code{dataset.gr} can also be a list of
#'   such GRanges objects.
#' @param field The name of the metadata field that contains readcounts. If no
#'   metadata field contains readcounts, and each range represents a single
#'   read, set to NULL.
#' @param ncores Number of cores to use for calculating coverage. For a single
#'   dataset, the max that will be used is 3, one for each possible strand
#'   (plus, minus, and unstranded). More cores can be used if \code{dataset.gr}
#'   is a list.
#'
#' @return A GRanges object with signal in the "score" metadata column. Note
#'   that the output is \emph{not} automatically converted into a
#'   \code{\link[BRGenomics:makeGRangesBRG]{"basepair-resolution"}} GRanges
#'   object.
#'
#' @author Mike DeBerardine
#' @seealso \code{\link[BRGenomics:makeGRangesBRG]{makeGRangesBRG}},
#'   \code{\link[IRanges:coverage-methods]{GenomicRanges::coverage}}
#' @export
#' @importFrom methods is
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
#' getStrandedCoverage(PROseq_paired, ncores = 1)[1:6]
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
#' getStrandedCoverage(ps_reads, field = NULL, ncores = 1)[1:6]
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
#' ps_cov <- getStrandedCoverage(PROseq, ncores = 1)
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
#' getStrandedCoverage(PROseq[idx], ncores = 1)
getStrandedCoverage <- function(dataset.gr, field = "score",
                                ncores = getOption("mc.cores", 2L)) {

    if (is.list(dataset.gr) || is(dataset.gr, "GRangesList")) {
        if (is.null(field))
            field <- list(NULL)
        return(mcMap(getStrandedCoverage, dataset.gr, field, ncores = 1L,
                     mc.cores = ncores))
    }

    if (!is.null(field) && !(field %in% names(mcols(dataset.gr))))
        stop(.nicemsg("The given value for 'field' is not found in
                      mcols(dataset.gr). If no field contains signal counts for
                      each range, set field = NULL"))

    cvg_ls <- mclapply(c("+", "-", "*"), .get_stranded_cvg, dataset.gr, field,
                       mc.cores = min(ncores, 3L))
    sort(do.call(c, cvg_ls))
}

#' @import GenomicRanges
.get_stranded_cvg <- function(strand.i, dataset.gr, field) {
    gr <- dataset.gr[strand(dataset.gr) == strand.i]
    if (is.null(field)) {
        cvg <- coverage(gr)
    } else {
        cvg <- coverage(gr, weight = mcols(gr)[[field]])
    }
    cvg_gr <- GRanges(cvg)
    seqinfo(cvg_gr) <- seqinfo(dataset.gr) # setting in GRanges fails sometimes
    strand(cvg_gr) <- strand.i # (setting within GRanges() fails if length=0)
    cvg_gr[abs(score(cvg_gr)) > 1e-12] # precision issues with coverage()
}


#' Randomly subsample reads from GRanges dataset
#'
#' Random subsampling is not performed on ranges, but on reads. Readcounts
#' should be given as a metadata field (usually "score"). This function can also
#' subsample ranges directly if \code{field = NULL}, but the \code{sample}
#' function can be used in this scenario.
#'
#' @param dataset.gr A GRanges object in which signal (e.g. readcounts) are
#'   contained within metadata, or a list of such GRanges objects.
#' @param n,prop Either the number of reads to subsample (\code{n}), or the
#'   proportion of total \emph{signal} to subsample (\code{prop}). Either
#'   \code{n} or \code{prop} can be given, but not both. If \code{dataset.gr} is
#'   a list, or if \code{length(field) > 1}, users can supply a vector or list
#'   of \code{n} or \code{prop} values to match the individual datasets, but
#'   care should be taken to ensure that a value is given for each and every
#'   dataset.
#' @param field The metadata field of \code{dataset.gr} that contains readcounts
#'   for reach position. If each range represents a single read, set \code{field
#'   = NULL}. If multiple fields are given, and \code{dataset.gr} is not a list,
#'   then \code{dataset.gr} will be treated as a multiplexed GRanges, and each
#'   field will be treated as an indpendent dataset. See
#'   \code{\link[BRGenomics:mergeGRangesData]{mergeGRangesData}}.
#' @param expand_ranges Logical indicating if ranges in \code{dataset.gr} should
#'   be treated as descriptions of single molecules (\code{FALSE}), or if ranges
#'   should be treated as representing multiple adjacent positions with the same
#'   signal (\code{TRUE}). See \code{\link[BRGenomics:getCountsByRegions]{
#'   getCountsByRegions}}.
#' @param ncores Number of cores to use for computations. Multicore only used
#'   when \code{dataset.gr} is a list, or if \code{length(field) > 1}.
#'
#' @return A GRanges object identical in format to \code{dataset.gr}, but
#'   containing a random subset of its data. If \code{field != NULL}, the length
#'   of the output cannot be known \emph{a priori}, but the sum of its score
#'   can.
#'
#' @section Use with normalized readcounts: If the metadata field contains
#'   normalized readcounts, an attempt will be made to infer the normalization
#'   factor based on the lowest signal value found in the specified field.
#'
#' @author Mike DeBerardine
#' @export
#' @importFrom methods is
#' @importFrom GenomicRanges mcols mcols<- countOverlaps
#' @importFrom parallel mcMap
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
subsampleGRanges <- function(dataset.gr, n = NULL, prop = NULL,
                             field = "score", expand_ranges = FALSE,
                             ncores = getOption("mc.cores", 2L)) {

    .check_xor_args(n, prop)

    if (is.list(dataset.gr) || is(dataset.gr, "GRangesList")) {
        if (is.null(field))  field <- list(NULL)
        if (is.null(n))  n <- list(NULL)
        if (is.null(prop))  prop <- list(NULL)

        mcMap(.subsample_gr, dataset.gr, n, prop, field, expand_ranges,
              ncores = 1L, mc.cores = ncores)

    } else if (length(field) > 1L) {
        if (is.null(n))  n <- list(NULL)
        if (is.null(prop))  prop <- list(NULL)

        grl <- mcMap(.subsample_gr, list(dataset.gr), n, prop, field,
                     expand_ranges, ncores = 1L, mc.cores = ncores)
        names(grl) <- field
        mergeGRangesData(grl, field = field, multiplex = TRUE, ncores = ncores)

    } else {
        .subsample_gr(dataset.gr, n, prop, field, expand_ranges, ncores)

    }
}


.subsample_gr <- function(dataset.gr, n, prop, field, expand_ranges, ncores) {
    if (is.null(field)) {
        if (is.null(n)) n <- floor(prop * length(dataset.gr))
        return(sample(dataset.gr, n))
    }
    signal_counts <- mcols(dataset.gr)[[field]]

    if (expand_ranges)
        signal_counts <- signal_counts * width(dataset.gr)
    if (normed_signal <- !.close_int(signal_counts)) {
        lcm <- min(signal_counts[signal_counts != 0L])
        signal_counts <- .try_unnorm_signal(signal_counts, lcm)
    }
    nreads <- sum(signal_counts)
    if (is.null(n)) n <- round(prop * nreads)

    if (expand_ranges) {
        # 1. Sample read numbers (integers 1 to nreads)
        samp_reads <- sort(sample.int(nreads, n))
        # 2. For each range, get the associated read numbers
        csumreads <- cumsum(signal_counts)
        # 3. For each sampled read:
        #  - Get the associated range index -> idx.range
        #  - Get its number within each range (e.g. "2nd read in that range")
        #  - Get its position within each range (offset from start, in bp)
        idx.range <- findInterval(samp_reads + 1L, csumreads)
        read_in_range <- samp_reads - ( c(0L, csumreads)[idx.range] )
        pos_in_range <- ceiling( read_in_range / signal_counts[idx.range] )

        gr_sample <- shift(resize(dataset.gr[idx.range], 1L), pos_in_range)
        gr_out <- getStrandedCoverage(gr_sample, field = NULL, ncores = ncores)
        names(mcols(gr_out)) <- field
    } else {
        # avoid modifying GRanges, and sample associated indices (w/o replace)
        idx <- rep.int(seq_along(dataset.gr), times = signal_counts)
        gr_sample <- dataset.gr[sample(idx, n)] # critical to not use 'prob' arg
        gr_out <- unique(gr_sample)
        mcols(gr_out)[field] <- countOverlaps(gr_out, gr_sample)
    }
    if (normed_signal)
        mcols(gr_out)[field] <- mcols(gr_out)[[field]] * lcm
    sort(gr_out)
}

.try_unnorm_signal <- function(signal_counts, lcm) {
    unnorm_signal <- signal_counts / lcm
    if (!.close_int(unnorm_signal))
        stop(.nicemsg("Signal given in 'field' are not whole numbers, and
                      unable to infer a normalization factor."))

    warning(.nicemsg("Signal given in 'field' are not whole numbers. A
                     normalization factor was inferred based on the lowest
                     signal value."))
    as.integer(unnorm_signal)
}


#' Merge GRanges objects
#'
#' Merges 2 or more GRanges objects by combining all of their ranges and
#' associated signal (e.g. readcounts). If \code{multiplex = TRUE}, the input
#' datasets are reversibly combined into a multiplexed GRanges containing a
#' field for each input dataset.
#'
#' @param ... Any number of GRanges objects in which signal (e.g. readcounts)
#'   are contained within metadata. Lists of GRanges can also be passed, but
#'   they must be named lists if \code{multiplex = TRUE}. Multiple lists can be
#'   passed, but if any inputs are lists, then all inputs must be lists.
#' @param field One or more \emph{input} metadata fields to be combined,
#'   typically the "score" field. Fields typically contain coverage information.
#'   If only a single field is given (i.e. all input GRanges use the same
#'   field), that same field will be used for the output. Otherwise, the
#'   \code{"score"} metadata field will be used by default. The output metadata
#'   fields are different if \code{multiplex} is enabled.
#' @param multiplex When set to \code{FALSE} (the default), input GRanges are
#'   merged irreversibly into a single new GRange, effectively combining the
#'   reads from different experiments. When \code{multiplex = TRUE}, the input
#'   GRanges data are reversibly combined into a multiplexed GRanges object,
#'   such that each input GRanges object has its own metadata field in the
#'   output.
#' @param makeBRG If \code{TRUE} (the default), the output GRanges will be made
#'   "basepair-resolution" via \code{\link[BRGenomics:makeGRangesBRG]{
#'   makeGRangesBRG}}. This is always set to code{FALSE} if
#'   \code{exact_overlaps = TRUE}.
#' @param exact_overlaps By default (\code{FALSE}), any ranges that cover more
#'   than a single base are treated as "coverage" signal (see
#'   \code{\link[BRGenomics:getCountsByRegions]{getCountsByRegions}}). If
#'   \code{exact_overlaps = TRUE}, all input ranges are conserved exactly as
#'   they are; ranges are only combined if they overlap exactly, and the signal
#'   of any combined range is the sum of the input ranges that were merged.
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
#' @importFrom methods is as
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
#' length(mergeGRangesData(ps_1, ps_2, ncores = 1))
#' length(mergeGRangesData(ps_1, ps_2, ps_3, ncores = 1))
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
#' mergeGRangesData(gr1, gr2, ncores = 1)
#'
#' #--------------------------------------------------#
#' # multiplex separate PRO-seq experiments
#' #--------------------------------------------------#
#'
#' multi.gr <- mergeGRangesData(gr1, gr2, multiplex = TRUE, ncores = 1)
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
                             makeBRG = TRUE, exact_overlaps = FALSE,
                             ncores = getOption("mc.cores", 2L)) {
    data_in <- list(...)
    is_list <- function(x) is.list(x) || is(x, "GRangesList")
    if (any(vapply(data_in, is_list, logical(1L))))
        data_in <- unlist(data_in)

    if (is(data_in, "GRangesList")) # no performance difference w/ unlist, etc.
        data_in <- as(data_in, "list")

    # check field names match inputs
    if (length(field) == 1L) {
        if (is.null(field))
            field <- "score"
        field <- rep.int(field, length(data_in))
    } else if (length(field) != length(data_in)) {
        stop("Given fields not equal to number of datasets")
    }

    if (multiplex) {
        if (is.null(names(data_in))) {
            exclude <- c("field", "multiplex", "ncores")
            in.names <- as.list(match.call())[-1L]
            names(data_in) <- in.names[!names(in.names) %in% exclude]
        }
        .multiplex_gr(data_in, field, ncores)
    } else if (exact_overlaps) {
        .merge_gr_exact(data_in, field, ncores)
    } else {
        .merge_gr(data_in, field, makeBRG, ncores)
    }
}

#' @importFrom parallel mcMap
#' @importFrom GenomicRanges mcols mcols<-
.merge_gr <- function(data_in, field, makeBRG, ncores) {

    if (any(field != "score"))
        data_in <- mcMap(.rename_field_score, data_in, field, mc.cores = ncores)

    gr <- getStrandedCoverage(do.call(c, c(data_in, use.names = FALSE)),
                              field = "score", ncores = ncores)

    # keep field name if only one given
    if (length(unique(field)) == 1L)
        names(mcols(gr)) <- field[1L]

    if (makeBRG) return(makeGRangesBRG(gr, ncores = ncores))
    gr
}

#' @importFrom GenomicRanges mcols mcols<-
.rename_field_score <- function(gr, field) {
    mcols(gr) <- mcols(gr)[field]
    names(mcols(gr)) <- "score"
    gr
}

#' @importFrom parallel mcMap
#' @importFrom GenomicRanges sort duplicated findOverlaps score score<-
#' @importFrom S4Vectors from to
.merge_gr_exact <- function(data_in, field, ncores) {

    if (any(field != "score"))
        mcMap(.rename_field_score, data_in, field, mc.cores = ncores)

    gr <- do.call(c, c(data_in, use.names = FALSE))
    dupes <- duplicated(gr)
    gr_out <- sort(gr[!dupes])
    gr_dp <- gr[dupes]

    # get exact hits to redundant ranges, and aggregate their scores
    hits <- findOverlaps(gr_out, gr_dp, type = "equal")
    dpsc <- aggregate(score(gr_dp)[to(hits)], by = list(from(hits)), FUN = sum)
    names(dpsc) <- c("idx", "score")

    # add the scores from the redundant ranges
    score(gr_out)[dpsc$idx] <- score(gr_out)[dpsc$idx] + dpsc$score
    gr_out
}


#' @importFrom parallel mclapply mcmapply
#' @importFrom GenomicRanges width sort mcols mcols<-
.multiplex_gr <- function(data_in, field, ncores) {

    # data must be *sorted*, base-pair resolution coverage data
    if (all( vapply(data_in, function(x) all(width(x) == 1L), logical(1L)) )) {
        data_in <- mclapply(data_in, sort, mc.cores = ncores)
    } else {
        warning(.nicemsg("One or more inputs are not 'basepair resolution
                         GRanges' objects. Coercing them using
                         makeGRangesBRG()..."), immediate. = TRUE)
        data_in <- mclapply(data_in, makeGRangesBRG, mc.cores = ncores)
    }

    # merge ranges
    gr <- do.call(c, c(data_in, use.names = FALSE))
    mcols(gr) <- NULL
    gr <- unique(sort(gr))

    # (Fastest to keep these next steps separated, esp. for large datasets)

    # get dataframe of signal counts for each dataset in the output GRanges
    idx <- mclapply(data_in, function(x) which(gr %in% x), mc.cores = ncores)
    counts <- mcmapply(function(dat, idx, field) {
        out <- rep.int(0L, length(gr))
        out[idx] <- mcols(dat)[[field]]
        out
    }, data_in, idx, field, mc.cores = ncores, SIMPLIFY = TRUE)

    mcols(gr)[names(data_in)] <- counts
    gr
}



#' Merge replicates of basepair-resolution GRanges objects
#'
#' This simple convenience function uses
#' \code{\link[BRGenomics:mergeGRangesData]{ mergeGRangesData}} to combine
#' replicates (e.g. biological replicates) of basepair-resolution GRanges
#' objects.
#'
#' @param ... Either a list of GRanges objects, or any number of GRanges objects
#'   (see \code{\link[BRGenomics:mergeGRangesData]{mergeGRangesData}}). However,
#'   the names of the datasets must end in \code{"_rep#"}, where "#" is one or
#'   more characters indicating the replicate.
#' @param field The metadata field that contains count information for each
#'   range. \code{length(field)} should either be 1, or equal to the number of
#'   datasets.
#' @param sample_names Optional character vector with which to rename the
#'   datasets. This is useful if the sample names do not conform to the "_rep"
#'   naming scheme.
#' @param makeBRG,exact_overlaps See \code{\link[BRGenomics:mergeGRangesData]{
#'   mergeGRangesData}}.
#' @param ncores The number of cores to use. This function will try to maximize
#'   the use of the \code{ncores} given, but care should be taken as
#'   \code{mergeGRangesData} can be memory intensive. Excessive memory usage can
#'   cause dramatic reductions in performance.
#'
#' @return A list of GRanges objects.
#' @export
#'
#' @examples
#' data("PROseq")
#' ps_list <- list(a_rep1 = PROseq[seq(1, length(PROseq), 4)],
#'                 b_rep1 = PROseq[seq(2, length(PROseq), 4)],
#'                 a_rep2 = PROseq[seq(3, length(PROseq), 4)],
#'                 b_rep2 = PROseq[seq(4, length(PROseq), 4)])
#' mergeReplicates(ps_list, ncores = 1)
mergeReplicates <- function(..., field = "score", sample_names = NULL,
                            makeBRG = TRUE, exact_overlaps = FALSE,
                            ncores = getOption("mc.cores", 2L)) {
    data_in <- list(...)
    if (any(vapply(data_in, is.list, logical(1))))
        data_in <- unlist(data_in)

    if (!is.null(sample_names))
        names(data_in) <- sample_names

    if (is.null(names(data_in))) {
        in.names <- as.list(match.call())[-1L]
        names(data_in) <- setdiff(in.names, c("field", "multiplex", "ncores"))
    }

    if (!all(grepl("_rep.+", names(data_in), perl = TRUE)))
        stop(.nicemsg("All sample names must end with \"_rep\" followed by one
                      or more characters that indicate the replicate"))

    if (length(field) == 1L)
        field <- rep.int(field, length(data_in))

    basenames <- unique(sub("_rep.*", "", names(data_in)))

    nreps <- length(data_in) / length(basenames)
    if (length(basenames) > nreps) {
        ncores_out <- ncores
        ncores_in <- 1L
    } else {
        ncores_out <- 1L
        ncores_in <- ncores
    }

    mrg <- function(basename) {
        idx <- grep(basename, names(data_in))
        mergeGRangesData(data_in[idx], field = field[idx], makeBRG = makeBRG,
                         exact_overlaps = exact_overlaps, ncores = ncores_in)
    }
    data_out <- mclapply(basenames, mrg, mc.cores = ncores_out)
    names(data_out) <- basenames
    data_out
}

