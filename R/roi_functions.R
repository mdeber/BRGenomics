
# for updating documentation:
#   for when fix_start != fix_end, filters intuitively for minimum window size
#   if fix_start == fix_end, the size of all the returned ranges will be the same;
      # for fix_start = "end", must be min_window between annotated
      #     start and start of interval;
      # for fix_end = "start", must be min_window between interval end and annotated
      #     end
#       but could make it so that...
#       min_window filters by the length within the annotated region that's
#       outside the returned range (i.e. for fix_end = "start", the distance
#       between the interval end to the annotated end; for fix_start = "end",
#       the distance between the annotated start to the beginning of the interval)
#


#' Extract Genebodies
#'
#' Analagous to \code{\link[GenomicRanges:promoters]{GenomicRanges::promoters}},
#' this function returns ranges that start and end downstream of the TSS. When
#' \code{fix = "start"}, the function behaves differently depending on the sign
#' of the \code{end} parameter. Currently, there's no way to set the ranges to
#' continue a fixed distance past the original end sites when \code{fix =
#' "start"}.
#'
#' @param genelist A GRanges object containing genes of interest.
#' @param start When \code{fix = "start"}, the distance downstream of the TSS
#'   where the new ranges should begin.
#' @param end Where the ranges should end. When \code{end = 0}, the returned
#'   ranges keep the original end sites. When \code{end < 0}, the new ranges
#'   will end \code{abs(end)} number of bases before the original end site. When
#'   \code{end > 0} and \code{fix = "start"}, the new ends are fixed to
#'   \code{end} number of bases from the \emph{original} start site.
#' @param min.window The minimum size of a returned genebody length, after
#'   accounting for \code{start} and \code{end} parameters.
#' @param fix If set to "end", function will return ranges centered around the
#'   end of the ranges. Negative values for \code{start} will begin the output
#'   ranges a fixed distance before the end of the input ranges; positive values
#'   will begin the ranges after the ends of the input ranges. The behavior of
#'   \code{end} is the same as for \code{start}, and errors will be returned if
#'   \code{end < start}.
#'
#' @return A GRanges object that may be shorter than \code{genelist} due to loss
#'   of short ranges.
#' @author Mike DeBerardine
#' @export
genebodies <- function(genelist,
                       start = 300,
                       end = -300,
                       fix_start = "start",
                       fix_end = "end",
                       min_window = 500) {

    if (!all(c(fix_start, fix_end) %in% c("start", "end"))) {
        stop("fix_start and fix_end must be 'start' or 'end'")
        return(geterrmessage())
    }

    if (fix_start == "end" & fix_end == "start") {
        stop("cannot have fix_end = start when fix_start = end")
        return(geterrmessage())
    }

    if ((fix_start == fix_end) & (end <= start)) {
        stop("If fix_end = fix_start, end must be greater than start")
        return(geterrmessage())
    }

    if (any(as.character(strand(genelist)) == "*")) {
        warning("Unstranded ranges were found and removed from genelist")
        genelist <- subset(genelist, strand != "*")
    }

    # Filter genelist based on min_window

    if (fix_start == "start") {
        if (fix_end == "end") {
            min_width <- start - end + min_window
        } else {
            # filter by distance from end of new interval to annotated end
            min_width <- end + min_window
        }
    } else {
        # filter by distance from annotated start to beginning of interval
        min_width <- min_window - start
    }

    genelist <- subset(genelist, width >= min_width)

    # Get genebodies

    # starts are at strand-specific beginnings of the genebodies
    sense_starts <- start(resize(genelist, 1, fix = fix_start))
    sense_ends <- start(resize(genelist, 1, fix = fix_end))

    # shift with strand-specificity
    is_plus <- as.character(strand(genelist)) == "+"
    sense_starts <- ifelse(is_plus, sense_starts + start, sense_starts - start)
    sense_ends <- ifelse(is_plus, sense_ends + end, sense_ends - end)

    # final starts/ends must be flipped for minus-strand genes
    starts <- ifelse(is_plus, sense_starts, sense_ends)
    ends <- ifelse(is_plus, sense_ends, sense_starts)

    ranges(genelist) <- IRanges(start = starts, end = ends)
    return(genelist)
}





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
#' @param field The metadata field of \code{dataset.gr} to be counted.
#' @param keep.score Logical indicating if the signal value at the max site
#'   should be reported. If set to \code{TRUE}, the values are kept as a new
#'   metadata column in \code{regions.gr}.
#'
#' @return Output is a GRanges object with regions.gr metadata, but each range
#'   only contains the site within each \code{regions.gr} range that had the
#'   most signal. If \code{binsize != 1}, a single site is still returned, but
#'   its position is set to the center of the bin. If the binsize is even, the
#'   site is rounded to be closer to the beginning of the range. If
#'   \code{keep.score = TRUE}, then the output will also have metadata for score
#'   at the max site. The output is \emph{not} necessarily same length as
#'   \code{regions.gr}, as regions without signal are not returned. If \emph{no
#'   regions} have signal (e.g. as could happen if running this function on a
#'   single region), the function will return an empty GRanges object with
#'   intact metadata columns.
#' @author Mike DeBerardine
#' @export
getMaxPositionsBySignal <- function(regions.gr,
                                    dataset.gr,
                                    binsize = 1,
                                    field = "score",
                                    keep.score = F) {

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
                                field = field,
                                remove_empty = FALSE)

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
        bin_centers <- max_pos
    } else {
        bin_centers <- floor(binsize / 2) + ( binsize * (max_pos - 1) )
    }

    regions.max.gr <- GenomicRanges::promoters(regions.gr, 0, bin_centers)
    regions.max.gr <- GenomicRanges::resize(regions.max.gr, 1, fix = "end")

    if (keep.score) {
        regions.max.gr$MaxSiteScore <- max_scores
    }

    return(regions.max.gr)
}



#' Subset regions of interest by highest signal
#'
#' Subsets regions based on signal in a dataset, taking only the top quantile of
#' regions.
#'
#' @param regions.gr A GRanges object containing regions of interest.
#' @param dataset.gr A GRanges object in which signal is contained in metadata
#'   (typically in the "score" field).
#' @param regions_quantile The proportion of regions.gr to return, e.g. if
#'   \code{regions_quantile = 0.2}, the top 20\% of regions by signal are
#'   returned.
#' @param field The metadata field of \code{dataset.gr} to be counted.
#' @param order_by_rank Logical indicating if genes should be returned in order
#'   of their expression. If \code{FALSE} (the default), genes are sorted by
#'   their positions.
#' @param density A logical indicating whether signal counts should be
#'   normalized to the width of ranges in \code{regions.gr}. By default, the
#'   function only considers the total signal in each range.
#'
#' @return A GRanges object of length \code{length(regions.gr) *
#'   regions_quantile}.
#' @author Mike DeBerardine
#' @export
subsetRegionsBySignal <- function(regions.gr,
                                  dataset.gr,
                                  regions_quantile,
                                  field = "score",
                                  order_by_rank = FALSE,
                                  density = FALSE) {
    # number of genes to return
    num.genes <- floor(regions_quantile * length(regions.gr))

    # find signal in genes
    hits <- findOverlaps(regions.gr, dataset.gr)
    signal.df <- aggregate(mcols(dataset.gr)[[field]][hits@to],
                           by = list(hits@from),
                           FUN = sum)
    names(signal.df) <- c("gene", "signal")

    if (density == TRUE) {
        signal.df$signal <- signal.df$signal / width(regions.gr[signal.df$gene])
    }

    # sort by signal, and return only the top num.genes
    signal.df <- signal.df[ order(signal.df$signal, decreasing = T), ]
    top.regions.gr <- regions.gr[ signal.df$gene[seq_len(num.genes)] ]

    if (order_by_rank) {
        return(top.regions.gr)
    } else {
        return(sort(top.regions.gr))
    }
}




