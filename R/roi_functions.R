# ========================================================================= #
# Functions for modifying regions of interest (i.e. annotations)
# ------------------------------------------------------------------------- #

#' Extract Genebodies
#'
#' This function returns ranges that are defined relative to the strand-specific
#' start and end sites of regions of interest (usually genes). Unlike
#' \code{\link[GenomicRanges:promoters]{GenomicRanges::promoters}}, distances
#' can be upstream or downstream based on the sign, and both the start and end
#' of the returned regions can be defined in terms of either the start of end
#' site of the input ranges. For example, \code{genebodies(txs, -50, 150,
#' fix.end = "start")} is equivalent to \code{promoters(txs, 50, 150)}. The
#' default arguments return ranges that begin 300 bases downstream of the
#' original start positions, and end 300 bases upstream of the original end
#' positions.
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
#' @return A GRanges object that may be shorter than \code{genelist} due to loss
#'   of short ranges.
#' @author Mike DeBerardine
#' @seealso \code{\link[GenomicRanges:intra-range-methods]{intra-range-methods}}
#' @export
genebodies <- function(genelist,
                       start = 300,
                       end = -300,
                       fix.start = "start",
                       fix.end = "end",
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
    sense_starts <- ifelse(is_plus, sense_starts + start, sense_starts - start)
    sense_ends <- ifelse(is_plus, sense_ends + end, sense_ends - end)

    # actual GRanges starts/ends must be flipped for minus-strand genes
    starts <- ifelse(is_plus, sense_starts, sense_ends)
    ends <- ifelse(is_plus, sense_ends, sense_starts)

    ranges(genelist) <- IRanges(start = starts, end = ends)
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
#' @param regions.gr A GRanges object containing regions of interest.
#' @param dataset.gr A GRanges object in which signal is contained in metadata
#'   (typically in the "score" field).
#' @param binsize The size of bin in which to calculate signal scores.
#' @param bin.centers Logical indicating if the centers of bins are returned, as
#' opposed to the entire bin. If \code{TRUE},
#' @param field The metadata field of \code{dataset.gr} to be counted.
#' @param keep.score Logical indicating if the signal value at the max site
#'   should be reported. If set to \code{TRUE}, the values are kept as a new
#'   metadata column in \code{regions.gr}.
#'
#' @return Output is a GRanges object with regions.gr metadata, but each range
#'   only contains the site within each \code{regions.gr} range that had the
#'   most signal. If \code{binsize > 1}, the entire bin is returned, unless
#'   \code{bin.centers = TRUE}, in which case a single-base site is returned.
#'   The site is set to the center of the bin, and if the binsize is even, the
#'   site is rounded to be closer to the beginning of the range.
#'
#'   If \code{keep.score = TRUE}, the output will also contain metadata for the
#'   signal at the max site. The output is \emph{not} necessarily same length as
#'   \code{regions.gr}, as regions without signal are not returned. If \emph{no
#'   regions} have signal (e.g. as could happen if running this function on a
#'   single region), the function will return an empty GRanges object with
#'   intact metadata columns.
#'
#' @author Mike DeBerardine
#' @seealso \code{\link[BRGenomics:getCountsByPositions]{getCountsByPositions}}
#' @export
getMaxPositionsBySignal <- function(regions.gr,
                                    dataset.gr,
                                    binsize = 1,
                                    bin.centers = FALSE,
                                    field = "score",
                                    keep.score = FALSE) {

    # keep only ranges with signal
    regions.gr <- subsetByOverlaps(regions.gr, dataset.gr)

    # if no regions in regions.gr have signal, return an empty GRanges object
    if (length(regions.gr) == 0) {
        if (keep.score)  regions.gr$MaxSiteScore <- integer(0)
        return(regions.gr)
    }

    multi_width <- FALSE
    if (length(unique(width(regions.gr))) > 1) {
        multi_width <- TRUE
        # currently faster to simply expand all regions for initial counting
        widths <- width(regions.gr) # store widths
        suppressWarnings(
            regions.gr <- resize(regions.gr, max(width(regions.gr)))
        )
    }

    mat <- getCountsByPositions(dataset.gr = dataset.gr,
                                regions.gr = regions.gr,
                                binsize = binsize,
                                field = field)

    # Get vector with max bin for each gene, and another for the scores therein
    if (multi_width) {
        bins_i <- floor(widths / binsize) # number of bins within each region
        # remove last bins (if widths/binsize gives remainder)
        countslist <- lapply(1:nrow(mat),
                             function(i) mat[ i, seq_len(bins_i[i]) ])
        max_pos <- vapply(countslist, which.max, FUN.VALUE = integer(1))
        max_scores <- vapply(countslist, max, FUN.VALUE = numeric(1))
    } else {
        max_pos <- apply(mat, 1, which.max)
        max_scores <- apply(mat, 1, max)
    }

    # Make new GRanges object with only max site for each gene
    regions.max.gr <- regions.gr

    if (binsize == 1) {
        bin.centers <- max_pos
        size <- 1
    } else {
        if (bin.centers) {
            bin.centers <- floor(binsize / 2) + ( binsize * (max_pos - 1) )
            size <- 1
        } else {
            bin.centers <- binsize * max_pos # end of bin
            size <- binsize
        }
    }

    regions.max.gr <- GenomicRanges::promoters(regions.gr, 0, bin.centers)
    regions.max.gr <- GenomicRanges::resize(regions.max.gr, size, fix = "end")

    if (keep.score)  regions.max.gr$MaxSiteScore <- max_scores

    return(regions.max.gr)
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
#'   regions to keep. Regions with signal quantiles below than the lower
#'   quantile are removed, while regions with signal quantiles above the upper
#'   quantile are removed. Quantiles must be in range \code{(0, 1)}. An empty
#'   GRanges object is returned if \code{lower quantile = 1} or \code{upper
#'   quantile = 0}.
#' @param field The metadata field of \code{dataset.gr} to be counted.
#' @param order.by.rank If \code{TRUE}, the output regions are sorted based on
#'   the amount of signal contained (in decreasing order). If \code{FALSE} (the
#'   default), genes are sorted by their positions.
#' @param density A logical indicating whether signal counts should be
#'   normalized to the width of ranges in \code{regions.gr}. By default, the
#'   function only considers the total signal in each range.
#'
#' @details Typical uses may include removing the 5% of genes with the lowest
#'   signal (\code{lower_quantile = 0.05}) and the 5% with the highest signal
#'   (\code{upper_quantile = 0.95}), or returning the middle 50% of genes by
#'   signal (\code{lower_quantile = 0.25}, \code{upper_quantile = 0.75}). If
#'   \code{lower_quantile = 0} and \code{upper_quantile = 1}, all regions are
#'   returned, but the returned regions will be sorted by position, or by score
#'   if \code{order.by.rank = TRUE}.
#'
#' @return A GRanges object of length \code{length(regions.gr) * (upper_quantile
#'   - lower_quantile)}.
#' @author Mike DeBerardine
#' @seealso \code{\link[BRGenomics:getCountsByRegion]{getCountsByRegion}}
#' @export
subsetRegionsBySignal <- function(regions.gr,
                                  dataset.gr,
                                  quantiles = c(0.5, 1),
                                  field = "score",
                                  order.by.rank = FALSE,
                                  density = FALSE) {

    if (quantiles[1] == 1 | quantiles[2] == 0)  return(regions.gr[0])

    signal_counts <- getCountsByRegions(dataset.gr = dataset.gr,
                                        regions.gr = regions.gr,
                                        field = field)

    if (density == TRUE)  signal_counts <- signal_counts / width(regions.gr)

    idx_rank <- order(signal_counts) # increasing
    bounds <- quantile(seq_along(regions.gr), quantiles)
    idx <- window(idx_rank, bounds[1], bounds[2])

    if (order.by.rank) {
        return(rev(regions.gr[idx]))
    } else {
        return(sort(regions.gr[idx]))
    }
}


