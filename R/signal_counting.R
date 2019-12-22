

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
#' @param ncores Multiple cores can only be used if \code{length(field) > 1}.
#'
#' @return Returns a vector the same length as \code{regions.gr} containing
#'   signal found in each range.
#' @author Mike DeBerardine
#' @export
getCountsByRegions <- function(dataset.gr,
                               regions.gr,
                               field = "score",
                               ncores = detectCores()) {
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
        # recursive call
        counts <- mclapply(field,
                           function(i) getCountsByRegions(dataset.gr,
                                                          regions.gr,
                                                          field = i),
                           mc.cores = ncores)
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
#' @param bin_FUN If \code{binsize > 1}, the function used to aggregate the
#'   signal at each base within each bin. The default is to sum the signal in
#'   each bin, but any function operating on a numeric vector can be used.
#' @param simplify_multi_widths A string indicating the output format if the
#'   ranges in \code{regions.gr} have variable widths. Default = \code{"list"}.
#'   See details below.
#' @param field The metadata field of \code{dataset.gr} to be counted. If
#'   \code{length(field) > 1}, the output is a list whose elements contain the
#'   output for generated each field.
#' @param ncores Multiple cores can only be used if \code{length(field) > 1}.
#'
#' @details If the widths of all ranges in \code{regions.gr} are equal, a matrix
#'   is returned containing a row for each range in \code{regions.gr}, and a
#'   column for each bin. For input \code{regions.gr} with varying widths,
#'   setting \code{simplify_multi_widths = "list"} will output a list of
#'   variable-length vectors, with each vector corresponding to an input region.
#'   If \code{simplify_multi_widths = "pad zero"} or \code{simplify_multi_widths
#'   = "pad NA"}, the output is a matrix containing a row for each range in
#'   \code{regions.gr}, and a column for each position in each range. The number
#'   of columns is determined by the largest range in \code{regions.gr}, and
#'   columns corresponding to positions outside of each range are either set to
#'   \code{0} or \code{NA}, depending on the argument.
#'
#' @author Mike DeBerardine
#' @export
getCountsByPositions <- function(dataset.gr,
                                 regions.gr,
                                 binsize = 1,
                                 bin_FUN = sum,
                                 simplify_multi_widths = c("list",
                                                           "pad zero",
                                                           "pad NA"),
                                 field = "score",
                                 ncores = detectCores()) {

    # function makes practical use of single-width dataset.gr, but the output
    # is always valid regardless of the input data type
    if (any(width(dataset.gr) != 1)) dataset.gr <- makeGRangesBPres(dataset.gr)

    multi_width <- length(unique(width(regions.gr))) > 1 # (logical)

    if (length(field) > 1) {
        # recursive call
        reslist <- mclapply(field, function(i) {
            getCountsByPositions(dataset.gr = dataset.gr,
                                 regions.gr = regions.gr,
                                 binsize = binsize,
                                 bin_FUN = bin_FUN,
                                 simplify_multi_widths = simplify_multi_widths,
                                 field = i)
        }, mc.cores = ncores)
        names(reslist) <- field
        return(reslist)
    }

    if (multi_width) {

        # check 'simplify_multi_widths' argument
        if (missing(simplify_multi_widths))  simplify_multi_widths <- "list"
        if (!simplify_multi_widths %in% c("list", "pad zero", "pad NA")) {
            stop(message = .nicemsg("regions.gr has multiple widths, but an
                                    invalid argument for simplify_multi_widths
                                    was given. See documentation"))
            return(geterrmessage())
        }

        # expand all regions to be the same width
        widths <- width(regions.gr) # save widths
        suppressWarnings(
            regions.gr <- resize(regions.gr, max(widths))
        )
    }

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
        mat <- apply(mat, 1, function(x) .binVector(x, binsize = binsize,
                                                    FUN = bin_FUN))
        mat <- t(mat) # apply will cbind rather than rbind
    }

    if (multi_width) {
        nbins_i <- floor(widths / binsize) # number of bins within each region

        if (simplify_multi_widths == "list") {

            return(mapply(function(i, nbin) mat[i, 1:nbin],
                          1:nrow(mat), nbins_i))

        } else {
            arridx_pad <- vapply(nbins_i,
                                 function(nbin) seq_len(ncol(mat)) > nbin,
                                 logical(ncol(mat)))
            arridx_pad <- t(arridx_pad) # sapply/vapply cbinds the rows
            arridx_pad <- which(arridx_pad, arr.ind = TRUE)

            if (simplify_multi_widths == "pad zero") {
                mat[arridx_pad] <- 0
            } else {
                mat[arridx_pad] <- NA
            }
        }
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

    if (length_normalize) {
        counts_pr <- counts_pr / width(promoters.gr)
        counts_gb <- counts_gb / width(genebodies.gr)
    }

    if (remove_empty) {
        idx_pr <- which(counts_pr != 0)
        idx_gb <- which(counts_gb != 0)
        idx <- intersect(idx_pr, idx_gb)
        counts_pr <- counts_pr[idx]
        counts_gb <- counts_gb[idx]
    }

    return(counts_pr / counts_gb)
}





