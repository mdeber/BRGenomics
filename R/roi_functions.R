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
#'   \code{\link[GenomicRanges:intra-range-methods]{GenomicRanges::promoters}},
#'   distances can be defined to be upstream or downstream by changing the sign
#'   of the argument, and both the start and end of the returned regions can be
#'   defined in terms of the strand-specific start or end site of the input
#'   ranges. For example, \code{genebodies(txs, -50, 150, fix.end = "start")} is
#'   equivalent to \code{promoters(txs, 50, 151)} (the downstream edge is off by
#'   1 because \code{promoters} keeps the downstream interval closed). The
#'   default arguments return ranges that begin 300 bases downstream of the
#'   original start positions, and end 300 bases upstream of the original end
#'   positions.
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

    fix.start <- match.arg(fix.start, c("start", "end", "center"))
    fix.end <- match.arg(fix.end, c("end", "start", "center"))

    if (fix.start == "end" & fix.end == "start")
        stop("cannot have fix.end = start when fix.start = end")

    if ((fix.start == fix.end) & (end <= start))
        stop("If fix.end = fix.start, end must be greater than start")

    if ("*" %in% levels(droplevels(strand(genelist)))) {
        warning("Unstranded ranges were found and removed from genelist")
        genelist <- subset(genelist, strand != "*")
    }

    # Filter genelist based on min.window
    if (fix.start == "start" & fix.end == "end")
        min.window <- start - end + min.window
    genelist <- genelist[width(genelist) >= min.window]

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
    genelist
}



#' Intersect or reduce ranges according to gene names
#'
#' These functions divide up regions of interest according to associated names,
#' and perform an inter-range operation on them. \code{intersectByGene} returns
#' the "consensus" segment that is common to all input ranges, and returns no
#' more than one range per gene. \code{reduceByGene} collapses the input ranges
#' into one or more non-overlapping ranges that encompass all segments from the
#' input ranges.
#'
#' @param regions.gr A GRanges object containing regions of interest. If
#'   \code{regions.gr} has the class \code{list}, \code{GRangesList}, or
#'   \code{CompressedGRangesList}, it will be treated as if each list element is
#'   a gene, and the GRanges within are the ranges associated with that gene.
#' @param gene_names A character vector with the same length as
#'   \code{regions.gr}.
#' @param disjoin Logical. If \code{disjoin = TRUE}, the output GRanges is
#'   disjoint, and each output range can overlap only a single gene. If
#'   \code{FALSE}, segments from different genes can overlap.
#'
#' @return A GRanges object whose individual ranges are named for the associated
#'   gene.
#'
#' @details These functions modify regions of interest that have associated
#'   names, such that several ranges share the same name, e.g. transcripts with
#'   associated gene names. Both functions "combine" the ranges on a
#'   gene-by-gene basis.
#'
#'   \strong{intersectByGene}
#'
#'   \emph{For each unique gene}, the segment that overlaps \emph{all} input
#'   ranges is returned. If no single range can be constructed that overlaps all
#'   input ranges, no range is returned for that gene (i.e. the gene is
#'   effectively filtered).
#'
#'   In other words, for all the ranges associated with a gene, the
#'   most-downstream start site is selected, and the most upstream end site is
#'   selected.
#'
#'   \strong{reduceByGene}
#'
#'   \emph{For each unique gene}, the associated ranges are
#'   \code{\link[IRanges:inter-range-methods]{reduced}} to produce one or
#'   more non-overlapping ranges. The output range(s) are effectively a
#'   \code{union} of the input ranges, and cover every input base.
#'
#'   With \code{disjoin = FALSE}, no genomic segment is overlapped by more than
#'   one range \emph{of the same gene}, but ranges from different genes can
#'   overlap. With \code{disjoin = TRUE}, the output ranges are disjoint, and no
#'   genomic position is overlapped more than once. Any segment that overlaps
#'   more than one gene is removed, but any segment (i.e. any section of an
#'   input range) that overlaps only one gene is still maintained.
#'
#' @section Typical Uses:
#'
#'   A typical use for \code{intersectByGene} is to avoid transcript isoform
#'   selection, as the returned range is found in every isoform.
#'
#'   \code{reduceByGene} can be used to count any and all reads that overlap any
#'   part of a gene's annotation, but without double-counting any of them. With
#'   \code{disjoin = FALSE}, no reads will be double-counted for the same gene,
#'   but the same read can be counted for multiple genes. With \code{disjoin =
#'   TRUE}, no read can be double-counted.
#'
#' @author Mike DeBerardine
#'
#' @importFrom GenomicRanges split start end mcols<-
#' @importFrom IRanges IRanges
#' @importFrom methods is as
#' @export
#'
#' @examples
#' # Make example data:
#' #  Ranges 1 and 2 overlap,
#' #  Ranges 3 and 4 are adjacent
#' gr <- GRanges(seqnames = "chr1",
#'               ranges = IRanges(start = c(1, 3, 7, 10),
#'                                end = c(4, 5, 9, 11)))
#' gr
#'
#' #--------------------------------------------------#
#' # intersectByGene
#' #--------------------------------------------------#
#'
#' intersectByGene(gr, c("A", "A", "B", "B"))
#'
#' intersectByGene(gr, c("A", "A", "B", "C"))
#'
#' gr2 <- gr
#' end(gr2)[1] <- 10
#' gr2
#'
#' intersectByGene(gr2, c("A", "A", "B", "C"))
#'
#' intersectByGene(gr2, c("A", "A", "A", "C"))
#'
#' #--------------------------------------------------#
#' # reduceByGene
#' #--------------------------------------------------#
#'
#' # For a given gene, overlapping/adjacent ranges are combined;
#' #  gaps result in multiple ranges for that gene
#' gr
#'
#' reduceByGene(gr, c("A", "A", "A", "A"))
#'
#' # With disjoin = FALSE, ranges from different genes can overlap
#' gnames <- c("A", "B", "B", "B")
#' reduceByGene(gr, gnames)
#'
#' # With disjoin = TRUE, segments overlapping >1 gene are removed as well
#' reduceByGene(gr, gnames, disjoin = TRUE)
#'
#' # Will use one more example to demonstrate how all
#' #  unambiguous segments are identified and returned
#' gr2
#'
#' gnames
#' reduceByGene(gr2, gnames, disjoin = TRUE)
#'
#' #--------------------------------------------------#
#' # reduceByGene, then aggregate counts by gene
#' #--------------------------------------------------#
#'
#' # Consider if you did getCountsByRegions on the last output,
#' #  you can aggregate those counts according to the genes
#' gr2_redux <- reduceByGene(gr2, gnames, disjoin = TRUE)
#' counts <- c(5, 2, 3) # if these were the counts-by-regions
#' aggregate(counts ~ names(gr2_redux), FUN = sum)
#'
#' # even more convenient if using a melted dataframe
#' df <- data.frame(gene = names(gr2_redux),
#'                  reads = counts)
#' aggregate(reads ~ gene, df, FUN = sum)
#'
#' # can be extended to multiple samples
#' df <- rbind(df, df)
#' df$sample <- rep(c("s1", "s2"), each = 3)
#' df$reads[4:6] <- c(3, 1, 2)
#' df
#'
#' aggregate(reads ~ sample*gene, df, FUN = sum)
intersectByGene <- function(regions.gr, gene_names) {

    if (is.list(regions.gr))
        regions.gr <- as(regions.gr, "GRangesList")

    if (is(regions.gr, "GRangesList")) {
        names(regions.gr) <- gene_names
        regions.gr <- unlist(regions.gr)
        gene_names <- names(regions.gr)
    }

    grl <- split(regions.gr, gene_names)

    max.starts <- vapply(start(grl), max, numeric(1))
    min.ends <- vapply(end(grl), min, numeric(1))

    # initialize output GRanges; same order as grl
    idx <- which(!duplicated(gene_names))
    gr <- unlist(split(regions.gr[idx], gene_names[idx]))
    mcols(gr) <- NULL

    # remove genes for which no consensus is possible
    idx.drop <- which(max.starts >= min.ends)
    if (length(idx.drop) > 0) {
        gr <- gr[-idx.drop]
        grl <- grl[-idx.drop]
        max.starts <- max.starts[-idx.drop]
        min.ends <- min.ends[-idx.drop]
    }

    ranges(gr) <- IRanges(max.starts, min.ends)
    names(gr) <- names(grl)
    sort(gr)
}


#' @rdname intersectByGene
#' @param disjoin Logical. If \code{disjoin = TRUE}, the output GRanges is
#'   disjoint, and each output range will match a single gene name. If
#'   \code{FALSE}, segments from different genes can overlap.
#'
#' @importFrom GenomicRanges split reduce disjoin
#' @importFrom methods is as
#' @export
reduceByGene <- function(regions.gr, gene_names, disjoin = FALSE) {

    if (is.list(regions.gr))
        regions.gr <- as(regions.gr, "GRangesList")

    if (is(regions.gr, "GRangesList")) {
        names(regions.gr) <- gene_names
        regions.gr <- unlist(reduce(regions.gr))
    } else {
        regions.gr <- unlist(reduce(split(regions.gr, gene_names)))
    }

    if (disjoin) {
        gene_names <- names(regions.gr) # reserve
        regions.gr <- disjoin(regions.gr, with.revmap = TRUE)

        # drop ranges with multiple mappings, unless from the same gene
        gn <- lapply(regions.gr$revmap, function(i) unique(gene_names[i]))
        idx <- lengths(gn) == 1L # ranges to keep
        regions.gr <- regions.gr[idx]
        regions.gr$revmap <- NULL
        names(regions.gr) <- unlist(gn[idx])
    }
    regions.gr
}





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
#' @param expand_ranges Logical indicating if ranges in \code{dataset.gr} should
#'   be treated as descriptions of single molecules (\code{FALSE}), or if ranges
#'   should be treated as representing multiple adjacent positions with the same
#'   signal (\code{TRUE}). See \code{\link[BRGenomics:getCountsByRegions]{
#'   getCountsByRegions}}.
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
                                    keep.signal = FALSE,
                                    expand_ranges = FALSE) {
    # keep only ranges with signal
    regions.gr <- subsetByOverlaps(regions.gr, dataset.gr)

    # if no regions in regions.gr have signal, return an empty GRanges object
    if (length(regions.gr) == 0) {
        if (keep.signal)  regions.gr$MaxSiteSignal <- integer(0)
        return(regions.gr)
    }

    # Get list with 2 vectors: max bin for each gene, and score in max bin
    mw <- length(unique(width(regions.gr))) > 1L
    FUN <- if (mw) .get_maxsite_mw else .get_maxsite
    maxsites <- FUN(dataset.gr, regions.gr, binsize, field, expand_ranges)

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
    regions.max.gr
}


.get_maxsite <- function(dataset.gr, regions.gr, binsize, field,
                         expand_ranges) {

    mat <- getCountsByPositions(dataset.gr, regions.gr, binsize = binsize,
                                field = field, expand_ranges = expand_ranges)
    max_pos <- apply(mat, 1, which.max)
    max_scores <- apply(mat, 1, max)
    list(pos = max_pos, score = max_scores)
}


#' @importFrom GenomicRanges width resize
.get_maxsite_mw <- function(dataset.gr, regions.gr, binsize, field,
                            expand_ranges) {

    countslist <- getCountsByPositions(
        dataset.gr = dataset.gr, regions.gr = regions.gr, binsize = binsize,
        simplify.multi.widths = "list", field = field,
        expand_ranges = expand_ranges
    )
    max_pos <- vapply(countslist, which.max, integer(1))
    max_scores <- vapply(countslist, max, numeric(1))
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
#' @param expand_ranges Logical indicating if ranges in \code{dataset.gr} should
#'   be treated as descriptions of single molecules (\code{FALSE}), or if ranges
#'   should be treated as representing multiple adjacent positions with the same
#'   signal (\code{TRUE}). See \code{\link[BRGenomics:getCountsByRegions]{
#'   getCountsByRegions}}.
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
                                  density = FALSE, keep.signal = FALSE,
                                  expand_ranges = FALSE) {

    if (quantiles[1] == 1 || quantiles[2] == 0)
        return(regions.gr[0])

    signal_counts <- getCountsByRegions(dataset.gr, regions.gr, field = field,
                                        expand_ranges = expand_ranges)

    if (density)
        signal_counts <- signal_counts / width(regions.gr)
    if (keep.signal)
        regions.gr$Signal <- signal_counts

    idx_rank <- order(signal_counts) # increasing
    bounds <- quantile(seq_along(regions.gr), quantiles)
    idx <- window(idx_rank, bounds[1], bounds[2])

    if (order.by.rank)
        return(rev(regions.gr[idx]))
    sort(regions.gr[idx])
}
