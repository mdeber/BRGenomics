# ========================================================================= #
# Functions for modifying regions of interest (i.e. annotations)
# ------------------------------------------------------------------------- #

#' Extract Genebodies
#'
#' This function returns ranges that are defined relative to the strand-specific
#' start and end sites of regions of interest (usually genes).
#'
#' @param genelist A GRanges object containing genes of interest.
#' @param start Depending on \code{fix.start}, the distance from either the
#'   strand-specific start or end site to begin the returned ranges. If
#'   positive, the returned range will begin downstream of the reference
#'   position; negative numbers are used to return sites upstream of the
#'   reference. Set \code{start = 0} to return the reference position.
#' @param end Identical to the \code{start} argument, but defines the
#'   strand-specific end position of returned ranges. \code{end} must be
#'   downstream of \code{start}.
#' @param fix.start The reference point to use for defining the strand-specific
#'   start positions of returned ranges, either \code{"start"} or \code{"end"}.
#' @param fix.end The reference point to use for defining the strand-specific
#'   end positions of returned ranges, either \code{"start"} or \code{"end"}.
#'   Cannot be set to \code{"start"} if \code{fix.start = "end"}.
#' @param min.window When \code{fix.start = "start"} and \code{fix.end = "end"},
#'   \code{min.window} defines the minimum size (width) of a returned range.
#'   However, when \code{fix.end = fix.start}, all returned ranges have the same
#'   width, and \code{min.window} simply size-filters the input ranges.
#'
#' @return A GRanges object that may be shorter than \code{genelist} due to
#'   filtering of short ranges. For example, using the default arguments, genes
#'   shorter than 600 bp would be removed.
#'
#' @details Unlike
#'   \code{\link[GenomicRanges:promoters]{GenomicRanges::promoters}}, distances
#'   can be defined to be upstream or downstream by changing the sign of the
#'   argument, and both the start and end of the returned regions can be defined
#'   in terms of the strand-specific start or end site of the input ranges. For
#'   example, \code{genebodies(txs, -50, 150, fix.end = "start")} is equivalent
#'   to \code{promoters(txs, 50, 151)} (the downstream edge is off by 1 because
#'   \code{promoters} keeps the downstream interval closed). The default
#'   arguments return ranges that begin 300 bases downstream of the original
#'   start positions, and end 300 bases upstream of the original end positions.
#'
#' @author Mike DeBerardine
#' @seealso \code{\link[GenomicRanges:intra-range-methods]{intra-range-methods}}
#' @export
#' @importFrom GenomicRanges strand width resize ranges<-
#' @importFrom IRanges IRanges
#'
#' @examples
#' data("txs_dm6_chr4") # load included transcript data
#' txs <- txs_dm6_chr4[c(1, 2, 167, 168)]
#'
#' txs
#'
#' #--------------------------------------------------#
#' # genebody regions from 300 bp after the TSS to
#' # 300 bp before the polyA site
#' #--------------------------------------------------#
#'
#' genebodies(txs, 300, -300)
#'
#' #--------------------------------------------------#
#' # promoter-proximal region from 50 bp upstream of
#' # the TSS to 100 bp downstream of the TSS
#' #--------------------------------------------------#
#'
#' promoters(txs, 50, 101)
#'
#' genebodies(txs, -50, 100, fix.end = "start")
#'
#' #--------------------------------------------------#
#' # region from 500 to 1000 bp after the polyA site
#' #--------------------------------------------------#
#'
#' genebodies(txs, 500, 1000, fix.start = "end")
genebodies <- function(genelist, start = 300, end = -300,
                       fix.start = "start", fix.end = "end",
                       min.window = 0) {

    .check_args_genebodies(start, end, fix.start, fix.end)

    if (any(as.character(strand(genelist)) == "*")) {
        warning("Unstranded ranges were found and removed from genelist")
        genelist <- subset(genelist, strand != "*")
    }

    # Filter genelist based on min.window
    min_width <- .find_min_width(start, end, fix.start, fix.end, min.window)
    genelist <- subset(genelist, width >= min_width)

    # starts are at strand-specific beginnings of the genebodies
    sense_starts <- start(resize(genelist, 1, fix = fix.start))
    sense_ends <- start(resize(genelist, 1, fix = fix.end))

    # shift with strand-specificity
    is_plus <- as.character(strand(genelist)) == "+"
    idx.p <- which(is_plus)
    idx.m <- which(!is_plus)
    ranges(genelist)[idx.p] <- IRanges(start = sense_starts[idx.p] + start,
                                       end = sense_ends[idx.p] + end)
    ranges(genelist)[idx.m] <- IRanges(start = sense_ends[idx.m] - end,
                                       end = sense_starts[idx.m] - start)

    return(genelist)
}


.check_args_genebodies <- function(start, end, fix.start, fix.end) {
    if (!all(c(fix.start, fix.end) %in% c("start", "end"))) {
        stop("fix.start and fix.end must be 'start' or 'end'")
        return(geterrmessage())
    }

    if (fix.start == "end" & fix.end == "start") {
        stop("cannot have fix.end = start when fix.start = end")
        return(geterrmessage())
    }

    if ((fix.start == fix.end) & (end <= start)) {
        stop("If fix.end = fix.start, end must be greater than start")
        return(geterrmessage())
    }
}


.find_min_width <- function(start, end, fix.start, fix.end, min.window) {
    if (fix.start == "start" & fix.end == "end") {
        return(start - end + min.window)
    } else {
        return(min.window)
    }
}

# .find_min_width <- function(fix.start, fix.end, min.window) {
#     if (fix.start == "start") {
#         if (fix.end == "end") {
#             return(start - end + min.window)
#         } else {
#             # filter by distance from end of new interval to annotated end
#             return(end + min.window)
#         }
#     } else {
#         # filter by distance from annotated start to beginning of interval
#         min.window - start
#     }
# }



#' Find sites with max signal in regions of interest
#'
#' For each signal-containing region of interest, find the single site with the
#' most signal. Sites can be found at base-pair resolution, or defined for
#' larger bins.
#'
#' @param dataset.gr A GRanges object in which signal is contained in metadata
#'   (typically in the "score" field).
#' @param regions.gr A GRanges object containing regions of interest.
#' @param binsize The size of bin in which to calculate signal scores.
#' @param bin.centers Logical indicating if the centers of bins are returned, as
#'   opposed to the entire bin. By default, entire bins are returned.
#' @param field The metadata field of \code{dataset.gr} to be counted.
#' @param keep.signal Logical indicating if the signal value at the max site
#'   should be reported. If set to \code{TRUE}, the values are kept as a new
#'   \code{MaxSiteSignal} metadata column in the output GRanges.
#'
#' @return Output is a GRanges object with regions.gr metadata, but each range
#'   only contains the site within each \code{regions.gr} range that had the
#'   most signal. If \code{binsize > 1}, the entire bin is returned, unless
#'   \code{bin.centers = TRUE}, in which case a single-base site is returned.
#'   The site is set to the center of the bin, and if the binsize is even, the
#'   site is rounded to be closer to the beginning of the range.
#'
#'   The output may not be the same length as \code{regions.gr}, as regions
#'   without signal are not returned. If no regions have signal (e.g. as could
#'   happen if running this function on single regions), the function will
#'   return an empty GRanges object with intact metadata columns.
#'
#'   If \code{keep.signal = TRUE}, the output will also contain metadata for the
#'   signal at the max site, named \code{MaxSiteSignal}.
#'
#' @author Mike DeBerardine
#' @seealso \code{\link[BRGenomics:getCountsByPositions]{getCountsByPositions}}
#' @export
#' @importFrom GenomicRanges promoters resize
#' @importFrom IRanges subsetByOverlaps
#'
#' @examples
#' data("PROseq") # load included PROseq data
#' data("txs_dm6_chr4") # load included transcripts
#'
#' #--------------------------------------------------#
#' # first 50 bases of transcripts
#' #--------------------------------------------------#
#'
#' pr <- promoters(txs_dm6_chr4, 0, 50)
#' pr[1:3]
#'
#' #--------------------------------------------------#
#' # max sites
#' #--------------------------------------------------#
#'
#' getMaxPositionsBySignal(PROseq, pr[1:3], keep.signal = TRUE)
#'
#' #--------------------------------------------------#
#' # max sites in 5 bp bins
#' #--------------------------------------------------#
#'
#' getMaxPositionsBySignal(PROseq, pr[1:3], binsize = 5, keep.signal = TRUE)
getMaxPositionsBySignal <- function(dataset.gr, regions.gr, binsize = 1,
                                    bin.centers = FALSE, field = "score",
                                    keep.signal = FALSE) {
    # keep only ranges with signal
    regions.gr <- subsetByOverlaps(regions.gr, dataset.gr)

    # if no regions in regions.gr have signal, return an empty GRanges object
    if (length(regions.gr) == 0) {
        if (keep.signal)  regions.gr$MaxSiteSignal <- integer(0)
        return(regions.gr)
    }

    # Get list with 2 vectors: max bin for each gene, and score in max bin
    if (length(unique(width(regions.gr))) == 1) {
        maxsites <- .get_maxsite(dataset.gr, regions.gr, binsize, field)
    } else { # widths vary
        maxsites <- .get_maxsite_mw(dataset.gr, regions.gr, binsize, field)
    }

    # Make new GRanges object with only max site for each gene
    regions.max.gr <- regions.gr
    if (binsize == 1) {
        bin.centers <- maxsites$pos
        size <- 1
    } else {
        if (bin.centers) {
            bin.centers <- ceiling(binsize / 2) + (binsize * (maxsites$pos - 1))
            size <- 1
        } else {
            bin.centers <- binsize * maxsites$pos # end of bin
            size <- binsize
        }
    }
    regions.max.gr <- promoters(regions.gr, 0, bin.centers)
    regions.max.gr <- resize(regions.max.gr, size, fix = "end")

    if (keep.signal)  regions.max.gr$MaxSiteSignal <- maxsites$score
    return(regions.max.gr)
}


.get_maxsite <- function(dataset.gr, regions.gr, binsize, field) {

    mat <- getCountsByPositions(dataset.gr = dataset.gr,
                                regions.gr = regions.gr,
                                binsize = binsize,
                                field = field)

    max_pos <- apply(mat, 1, which.max)
    max_scores <- apply(mat, 1, max)

    list(pos = max_pos, score = max_scores)
}


#' @importFrom GenomicRanges width resize
.get_maxsite_mw <- function(dataset.gr, regions.gr, binsize, field) {

    # faster to simply expand all regions for initial counting
    widths <- width(regions.gr) # store widths
    suppressWarnings( regions.gr <- resize(regions.gr, max(width(regions.gr))) )

    mat <- getCountsByPositions(dataset.gr = dataset.gr,
                                regions.gr = regions.gr,
                                binsize = binsize,
                                field = field)
    # get number of bins for each region; & trim if remainder to widths/binsize
    bins_i <- floor(widths / binsize)
    countslist <- lapply(seq_len(nrow(mat)),
                         function(i) mat[ i, seq_len(bins_i[i]) ])

    max_pos <- vapply(countslist, which.max, FUN.VALUE = integer(1))
    max_scores <- vapply(countslist, max, FUN.VALUE = numeric(1))

    list(pos = max_pos, score = max_scores)
}



#' Subset regions of interest by quantiles of overlapping signal
#'
#' @description A convenience function to subset regions of interest by the
#'   amount of signal they contain, according to their quantile (i.e. their
#'   signal ranks).
#'
#' @param regions.gr A GRanges object containing regions of interest.
#' @param dataset.gr A GRanges object in which signal is contained in metadata
#'   (typically in the "score" field).
#' @param quantiles A value pair giving the lower quantile and upper quantile of
#'   regions to keep. Regions with signal quantiles below the lower quantile are
#'   removed, and likewise for regions with signal quantiles above the upper
#'   quantile. Quantiles must be in range \code{(0, 1)}. An empty GRanges object
#'   is returned if the lower quantile is set to \code{1} or if the upper
#'   quantile is set to \code{0}.
#' @param field The metadata field of \code{dataset.gr} to be counted, typically
#'   "score".
#' @param order.by.rank If \code{TRUE}, the output regions are sorted based on
#'   the amount of overlapping signal (in decreasing order). If \code{FALSE}
#'   (the default), genes are sorted by their positions.
#' @param density A logical indicating whether signal counts should be
#'   normalized to the width (chromosomal length) of ranges in
#'   \code{regions.gr}. By default, no length normalization is performed.
#' @param keep.signal Logical indicating if signal counts should be kept. If set
#'   to \code{TRUE}, the signal for each range (length-normalized if
#'   \code{density = TRUE}) are kept as a new \code{Signal} metadata column in
#'   the output GRanges object.
#'
#' @return A GRanges object of length \code{length(regions.gr) * (upper_quantile
#'   - lower_quantile)}.
#' @author Mike DeBerardine
#' @seealso \code{\link[BRGenomics:getCountsByRegions]{getCountsByRegions}}
#' @export
#' @importFrom stats quantile
#' @importFrom GenomicRanges width
#'
#' @examples
#' data("PROseq") # load included PROseq data
#' data("txs_dm6_chr4") # load included transcripts
#'
#' txs_dm6_chr4
#'
#' #--------------------------------------------------#
#' # get the top 50% of transcripts by signal
#' #--------------------------------------------------#
#'
#' subsetRegionsBySignal(txs_dm6_chr4, PROseq)
#'
#' #--------------------------------------------------#
#' # get the middle 50% of transcripts by signal
#' #--------------------------------------------------#
#'
#' subsetRegionsBySignal(txs_dm6_chr4, PROseq, quantiles = c(0.25, 0.75))
#'
#' #--------------------------------------------------#
#' # get the top 10% of transcripts by signal, and sort them by highest signal
#' #--------------------------------------------------#
#'
#' subsetRegionsBySignal(txs_dm6_chr4, PROseq, quantiles = c(0.9, 1),
#'                       order.by.rank = TRUE)
#'
#' #--------------------------------------------------#
#' # remove the most extreme 10% of regions, and keep scores
#' #--------------------------------------------------#
#'
#' subsetRegionsBySignal(txs_dm6_chr4, PROseq, quantiles = c(0.05, 0.95),
#'                       keep.signal = TRUE)
subsetRegionsBySignal <- function(regions.gr, dataset.gr, quantiles = c(0.5, 1),
                                  field = "score", order.by.rank = FALSE,
                                  density = FALSE, keep.signal = FALSE) {

    if (quantiles[1] == 1 | quantiles[2] == 0)  return(regions.gr[0])

    signal_counts <- getCountsByRegions(dataset.gr = dataset.gr,
                                        regions.gr = regions.gr,
                                        field = field)

    if (density)  signal_counts <- signal_counts / width(regions.gr)
    if (keep.signal)  regions.gr$Signal <- signal_counts

    idx_rank <- order(signal_counts) # increasing
    bounds <- quantile(seq_along(regions.gr), quantiles)
    idx <- window(idx_rank, bounds[1], bounds[2])

    if (order.by.rank) {
        return(rev(regions.gr[idx]))
    } else {
        return(sort(regions.gr[idx]))
    }
}


