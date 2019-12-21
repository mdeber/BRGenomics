

#' Signal counts in regions of interest
#'
#' Count signal (e.g. read coverage) of data in each region of interest.
#'
#' @param dataset.gr A GRanges object in which signal is contained in metadata
#'   (typically in the "score" field).
#' @param regions.gr A GRanges object containing all the regions of interest.
#' @param field The metadata field of \code{dataset.gr} to be counted. If
#'   \code{length(field) > 1}, a dataframe is returned containing the counts for
#'   each region in each field.
#'
#' @return Returns a vector the same length as \code{regions.gr} containing
#'   signal found in each range.
#' @author Mike DeBerardine
#' @export
getCountsByRegions <- function(dataset.gr, regions.gr, field = "score") {
    if (length(field) == 1) {
        hits <- findOverlaps(regions.gr, dataset.gr)
        counts <- aggregate(mcols(dataset.gr)[[field]][hits@to],
                            by = list(hits@from),
                            FUN = sum)
        names(counts) <- c("gene.idx", "signal")
        counts.all <- rep(0, length(regions.gr)) # include regions without hits
        counts.all[counts$gene.idx] <- counts$signal
        return(counts.all)
    } else {
        counts <- lapply(field, function(i) getCountsByRegions(dataset.gr,
                                                               regions.gr,
                                                               field = i))
        names(counts) <- field
        return(as.data.frame(counts))
    }
}



#' Get signal count matrix within regions of interest
#'
#' @param dataset.gr A GRanges object in which signal is contained in metadata
#'   (typically in the "score" field).
#' @param regions.gr A GRanges object containing all the regions of interest.
#'   All ranges must have the same width!
#' @param binsize Size of bins (in bp) to use for counting within each range of
#'   \code{regions.gr}. Note that counts will \emph{not} be length-normalized.
#' @param field The metadata field of \code{dataset.gr} to be counted.
#' @param remove_empty Logical indicating whether regions without signal should
#'   be removed.
#'
#' @return A matrix containing a row for each range in regions.gr, and a column
#'   for each bin.
#' @author Mike DeBerardine
#' @export
getCountsByPositions <- function(dataset.gr,
                                 regions.gr,
                                 binsize = 1,
                                 field = "score",
                                 remove_empty = FALSE) {
    if (length(unique(width(regions.gr))) > 1) {
        stop(message = "Not all ranges in regions.gr are the same width")
        return(geterrmessage())
    }

    # function requires widths of ranges in dataset.gr are also equal to 1
    if (length(unique(width(dataset.gr))) > 1)
        dataset.gr <- makeGRangesBPres(dataset.gr)

    if (remove_empty)  regions.gr <- subsetByOverlaps(regions.gr, dataset.gr)

    # initialize signal matrix of dim = (region, position within region)
    mat <- matrix(0, length(regions.gr), unique(width(regions.gr)))

    hits <- findOverlaps(regions.gr, dataset.gr)

    # find (x, y) = (region, position)
    x <- hits@from
    y1 <- start(dataset.gr[hits@to]) # site of signal
    y2 <- start(resize(regions.gr[hits@from], 1)) # beginning of window
    y <- abs(y1 - y2) + 1 # position of signal within region
    z <- mcols(dataset.gr)[[field]][hits@to]

    mat[cbind(x, y)] <- z

    if (binsize > 1) {
        mat <- apply(mat, 1, .binVector, binsize = binsize)
        mat <- t(mat) # apply will cbind rather than rbind
    }

    return(mat)
}



#' Calculate pausing indices from user-supplied promoters & genebodies
#'
#' Pausing index (PI) is calculated for each gene (within matched
#' \code{promoters.gr} and \code{genebodies.gr}) as promoter signal counts
#' divided by genebody signal counts. If \code{length.normalize = TRUE}
#' (recommended), the signal counts within each range in \code{promoters.gr} and
#' \code{genebodies.gr} are divided by their respective range widths (region
#' lengths) before pausing indices are calculated.
#'
#' @param dataset.gr A GRanges object in which signal is contained in metadata
#'   (typically in the "score" field).
#' @param promoters.gr A GRanges object containing all the regions of interest.
#'   The sum of all signal counts within is the pause index numerator.
#' @param genebodies.gr A GRanges object containing all the regions of interest.
#'   The sum of all signal counts within is the pause index denominator.
#' @param field The metadata field of \code{dataset.gr} to be counted, i.e. that
#'   contains the readcounts of interest.
#' @param length_normalize A logical indicating if signal counts within regions
#'   of interest should be length normalized. The default is TRUE, which is
#'   recommended, especially if input regions don't all have the same width.
#' @param remove_empty A logical indicating if genes without any signal should
#'   be removed. The default is FALSE.
#'
#' @return A vector of length given by the length of the genelist (or possibly
#'   shorter if \code{remove_empty = TRUE}).
#' @author Mike DeBerardine
#' @export
getPausingIndices <- function(dataset.gr,
                              promoters.gr,
                              genebodies.gr,
                              field = "score",
                              length_normalize = TRUE,
                              remove_empty = FALSE) {

    if (length(promoters.gr) != length(genebodies.gr)) {
        stop(message = .nicemsg("Number of ranges in promoters.gr != number of
                                ranges in genebodies.gr"))
        return(geterrmessage())
    }

    counts_pr <- getCountsByRegions(dataset.gr = dataset.gr,
                                    regions.gr = promoters.gr,
                                    field = field)
    counts_gb <- getCountsByRegions(dataset.gr = dataset.gr,
                                    regions.gr = genebodies.gr,
                                    field = field)

    if (remove_empty) {
        idx_pr <- which(counts_pr != 0)
        idx_gb <- which(counts_gb != 0)
        idx <- intersect(idx_pr, idx_gb)
        counts_pr <- counts_pr[idx]
        counts_gb <- counts_gb[idx]
    }

    if (length_normalize) {
        counts_pr <- counts_pr / width(promoters.gr)
        counts_gb <- counts_gb / width(genebodies.gr)
    }

    return(counts_pr / counts_gb)
}





