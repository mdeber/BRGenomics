

#' Get signal counts in regions of interest
#'
#' Get the sum of the signal in \code{dataset.gr} that overlaps each range in
#' \code{regions.gr}.
#'
#' @param dataset.gr A GRanges object in which signal is contained in metadata
#'   (typically in the "score" field).
#' @param regions.gr A GRanges object containing regions of interest.
#' @param field The metadata field of \code{dataset.gr} to be counted. If
#'   \code{length(field) > 1}, a dataframe is returned containing the counts for
#'   each region in each field.
#' @param ncores Multiple cores can only be used if \code{length(field) > 1}.
#'
#' @return An atomic vector the same length as \code{regions.gr} containing
#' the sum of the signal overlapping each range of \code{regions.gr}.
#'
#' @author Mike DeBerardine
#' @seealso \code{\link[BRGenomics:getCountsByPositions]{getCountsByPositions}}
#' @export
#' @importFrom stats aggregate
#' @importFrom GenomicRanges findOverlaps mcols
#' @importFrom parallel detectCores mclapply
#'
#' @examples
#' data("PROseq") # load included PROseq data
#' data("txs_dm6_chr4") # load included transcripts
#'
#' counts <- getCountsByRegions(PROseq, txs_dm6_chr4)
#'
#' length(txs_dm6_chr4)
#' length(counts)
#' head(counts)
#'
#' # Assign as metadata to the transcript GRanges
#' txs_dm6_chr4$PROseq <- counts
#'
#' txs_dm6_chr4[1:6]
getCountsByRegions <- function(dataset.gr, regions.gr, field = "score",
                               ncores = detectCores()) {
    if (length(field) == 1) {
        hits <- findOverlaps(regions.gr, dataset.gr)
        counts <- aggregate(mcols(dataset.gr)[[field]][hits@to],
                            by = list(hits@from), FUN = sum)
        names(counts) <- c("gene.idx", "signal")
        counts.all <- rep(0, length(regions.gr)) # include regions without hits
        counts.all[counts$gene.idx] <- counts$signal
        return(counts.all)

    } else {
        # recursive call
        call_each <- function(x) getCountsByRegions(dataset.gr, regions.gr, x)
        counts <- mclapply(field, call_each, mc.cores = ncores)
        names(counts) <- field
        return(as.data.frame(counts))
    }
}



#' Get signal counts at each position within regions of interest
#'
#' Get the sum of the signal in \code{dataset.gr} that overlaps each position
#' within each range in \code{regions.gr}. If binning is used (i.e. positions
#' are wider than 1 bp), any function can be used to summarize the signal
#' overlapping each bin.
#'
#' @param dataset.gr A GRanges object in which signal is contained in metadata
#'   (typically in the "score" field).
#' @param regions.gr A GRanges object containing regions of interest.
#' @param binsize Size of bins (in bp) to use for counting within each range of
#'   \code{regions.gr}. Note that counts will \emph{not} be length-normalized.
#' @param FUN If \code{binsize > 1}, the function used to aggregate the signal
#'   within each bin. By default, the signal is summed, but any function
#'   operating on a numeric vector can be used.
#' @param simplify.multi.widths A string indicating the output format if the
#'   ranges in \code{regions.gr} have variable widths. Default is \code{"list"}.
#'   See details below.
#' @param field The metadata field of \code{dataset.gr} to be counted. If
#'   \code{length(field) > 1}, the output is a list whose elements contain the
#'   output for generated each field.
#' @param ncores Multiple cores can only be used if \code{length(field) > 1}.
#'
#' @return If the widths of all ranges in \code{regions.gr} are equal, a matrix
#'   is returned that contains a row for each region of interest, and a column
#'   for each position (each base if \code{binsize = 1}) within each region.
#'
#' @section Use of multi-width regions of interest: If the input
#'   \code{regions.gr} contains ranges of varying widths, setting
#'   \code{simplify.multi.widths = "list"} will output a list of variable-length
#'   vectors, with each vector corresponding to an individual input region. If
#'   \code{simplify.multi.widths = "pad 0"} or \code{"pad NA"}, the output is a
#'   matrix containing a row for each range in \code{regions.gr}, but the number
#'   of columns is determined by the largest range in \code{regions.gr}. For
#'   each region of interest, columns that correspond to positions outside of
#'   the input range are set, depending on the argument, to \code{0} or
#'   \code{NA}.
#'
#' @author Mike DeBerardine
#' @seealso \code{\link[BRGenomics:getCountsByRegions]{getCountsByRegions}}
#' @export
#' @importFrom parallel detectCores mclapply
#' @importFrom GenomicRanges width
#'
#' @examples
#' data("PROseq") # load included PROseq data
#' data("txs_dm6_chr4") # load included transcripts
#'
#' #--------------------------------------------------#
#' # counts from 0 to 50 bp after the TSS
#' #--------------------------------------------------#
#'
#' txs_pr <- promoters(txs_dm6_chr4, 0, 50) # first 50 bases
#' countsmat <- getCountsByPositions(PROseq, txs_pr)
#' countsmat[10:15, 41:50] # show only 41-50 bp after TSS
#'
#' #--------------------------------------------------#
#' # redo with 10 bp bins from 0 to 100
#' #--------------------------------------------------#
#'
#' # column 5 is sums of rows shown above
#'
#' txs_pr <- promoters(txs_dm6_chr4, 0, 100)
#' countsmat <- getCountsByPositions(PROseq, txs_pr, binsize = 10)
#' countsmat[10:15, ]
#'
#' #--------------------------------------------------#
#' # same as the above, but with the average signal in each bin
#' #--------------------------------------------------#
#'
#' countsmat <- getCountsByPositions(PROseq, txs_pr, binsize = 10, FUN = mean)
#' countsmat[10:15, ]
#'
#' #--------------------------------------------------#
#' # standard deviation of signal in each bin
#' #--------------------------------------------------#
#'
#' countsmat <- getCountsByPositions(PROseq, txs_pr, binsize = 10, FUN = sd)
#' round(countsmat[10:15, ], 2)
getCountsByPositions <- function(dataset.gr, regions.gr, binsize = 1, FUN = sum,
                                 simplify.multi.widths = c("list",
                                                           "pad 0",
                                                           "pad NA"),
                                 field = "score", ncores = detectCores()) {

    # function makes practical use of single-width dataset.gr, but the output
    # is always valid regardless of the input data type
    if (any(width(dataset.gr) != 1)) dataset.gr <- makeGRangesBRG(dataset.gr)

    hits <- findOverlaps(regions.gr, dataset.gr) # get hits early (once)

    if (length(field) > 1) {

        reslist <- mclapply(field, function(field.i) {
            .get_cbp(hits, dataset.gr, regions.gr, binsize = binsize,
                     FUN = FUN, simplify.multi.widths, field = field.i)
        }, mc.cores = ncores)

        names(reslist) <- field
        return(reslist)

    } else {
        .get_cbp(hits, dataset.gr, regions.gr, binsize = binsize,
                 FUN = FUN, simplify.multi.widths, field = field)
    }
}


.get_cbp <- function(hits, dataset.gr, regions.gr, binsize,
                     FUN, simplify.multi.widths, field) {

    multi_width <- length(unique(width(regions.gr))) > 1

    if (multi_width) {
        .get_cbp_mw(hits, dataset.gr, regions.gr, binsize, FUN,
                    simplify.multi.widths, field)
    } else {
        .get_signal_mat(hits, dataset.gr, regions.gr, binsize, FUN, field)
    }

}


#' @importFrom GenomicRanges width resize
.get_cbp_mw <- function(hits, dataset.gr, regions.gr, binsize, FUN,
                        simplify.multi.widths, field) {

    simplify.multi.widths <- match.arg(simplify.multi.widths,
                                       c("list", "pad 0", "pad NA"))

    # expand all regions to be the same width, then get counts
    widths <- width(regions.gr) # save widths
    suppressWarnings( regions.gr <- resize(regions.gr, max(widths)) )
    mat <- .get_signal_mat(hits, dataset.gr, regions.gr, binsize, FUN, field)

    nbins <- floor(widths / binsize) # number of bins to keep for each region

    if (simplify.multi.widths == "list") {
        # list of vectors whose lengths determined by width of range i
        clist <- mapply(function(row.i, nbins.i) mat[row.i, seq_len(nbins.i)],
                        seq_len(nrow(mat)), nbins)
        return(clist)

    } else {
        arridx_pad <- vapply(nbins,
                             function(nbins.i) seq_len(ncol(mat)) > nbins.i,
                             FUN.VALUE = logical(ncol(mat)))
        arridx_pad <- t(arridx_pad) # sapply/vapply cbinds the rows
        arridx_pad <- which(arridx_pad, arr.ind = TRUE)
        if (simplify.multi.widths == "pad 0") {
            mat[arridx_pad] <- 0
        } else {
            mat[arridx_pad] <- NA
        }
        return(mat)
    }
}


#' @importFrom GenomicRanges start mcols
.get_signal_mat <- function(hits, dataset.gr, regions.gr, binsize, FUN, field) {

    # initialize signal matrix of dim = (region, position within region)
    mat <- matrix(0, length(regions.gr), unique(width(regions.gr)))

    # find (x, y) = (region, position)
    x <- hits@from
    y1 <- start(dataset.gr[hits@to]) # site of signal
    y2 <- start(resize(regions.gr[hits@from], 1)) # beginning of window
    y <- abs(y1 - y2) + 1 # position of signal within region
    z <- mcols(dataset.gr)[[field]][hits@to]
    mat[cbind(x, y)] <- z

    if (binsize > 1) {
        mat <- apply(mat, 1, function(x) .binVector(x, binsize = binsize,
                                                    FUN = FUN))
        mat <- t(mat) # apply will cbind rather than rbind
    }
    return(mat)
}


#' Calculate pausing indices from user-supplied promoters & genebodies
#'
#' Pausing index (PI) is calculated for each gene (within matched
#' \code{promoters.gr} and \code{genebodies.gr}) as promoter-proximal (or pause
#' region) signal counts divided by genebody signal counts. If
#' \code{length.normalize = TRUE} (recommended), the signal counts within each
#' range in \code{promoters.gr} and \code{genebodies.gr} are divided by their
#' respective range widths (region lengths) before pausing indices are
#' calculated.
#'
#' @param dataset.gr A GRanges object in which signal is contained in metadata
#'   (typically in the "score" field).
#' @param promoters.gr A GRanges object containing promoter-proximal regions of
#'   interest.
#' @param genebodies.gr A GRanges object containing genebody regions of
#'   interest.
#' @param field The metadata field of \code{dataset.gr} to be counted. If
#'   \code{length(field) > 1}, a dataframe is returned containing the pausing
#'   indices for each region in each field.
#' @param length.normalize A logical indicating if signal counts within regions
#'   of interest should be length normalized. The default is \code{TRUE}, which
#'   is recommended, especially if input regions don't all have the same width.
#' @param remove.empty A logical indicating if genes without any signal in
#'   \code{promoters.gr} should be removed. No genes are filtered by default.
#' @param ncores Multiple cores can only be used if \code{length(field) > 1}.
#'
#' @return A vector parallel to the input genelist, unless \code{remove.empty =
#'   TRUE}, in which case the vector may be shorter. If \code{length(field) >
#'   1}, a dataframe is returned, containing a column for each field.
#'
#' @author Mike DeBerardine
#' @seealso \code{\link[BRGenomics:getCountsByRegions]{getCountsByRegions}}
#' @export
#' @importFrom parallel detectCores
#'
#' @examples
#' data("PROseq") # load included PROseq data
#' data("txs_dm6_chr4") # load included transcripts
#'
#' #--------------------------------------------------#
#' # Get promoter-proximal and genebody regions
#' #--------------------------------------------------#
#'
#' # genebodies from +300 to 300 bp before the poly-A site
#' gb <- genebodies(txs_dm6_chr4, 300, -300, min.window = 400)
#'
#' # get the transcripts that are large enough (>1kb in size)
#' txs <- subset(txs_dm6_chr4, tx_name %in% gb$tx_name)
#'
#' # for the same transcripts, promoter-proximal region from 0 to +100
#' pr <- promoters(txs, 0, 100)
#'
#' #--------------------------------------------------#
#' # Calculate pausing indices
#' #--------------------------------------------------#
#'
#' pidx <- getPausingIndices(PROseq, pr, gb)
#'
#' length(txs)
#' length(pidx)
#' head(pidx)
#'
#' #--------------------------------------------------#
#' # Without length normalization
#' #--------------------------------------------------#
#'
#' head( getPausingIndices(PROseq, pr, gb, length.normalize = FALSE) )
#'
#' #--------------------------------------------------#
#' # Removing empty means the values no longer match the genelist
#' #--------------------------------------------------#
#'
#' pidx_signal <- getPausingIndices(PROseq, pr, gb, remove.empty = TRUE)
#'
#' length(pidx_signal)
getPausingIndices <- function(dataset.gr, promoters.gr, genebodies.gr,
                              field = "score", length.normalize = TRUE,
                              remove.empty = FALSE, ncores = detectCores()) {

    if (length(promoters.gr) != length(genebodies.gr)) {
        stop(message = .nicemsg("Number of ranges in promoters.gr != number of
                                ranges in genebodies.gr"))
        return(geterrmessage())
    }

    counts_pr <- getCountsByRegions(dataset.gr = dataset.gr,
                                    regions.gr = promoters.gr,
                                    field = field, ncores = ncores)
    counts_gb <- getCountsByRegions(dataset.gr = dataset.gr,
                                    regions.gr = genebodies.gr,
                                    field = field, ncores = ncores)

    if (length.normalize) {
        if (length(field) > 1) {
            counts_pr <- .lnorm_multifields(counts_pr, promoters.gr, field)
            counts_gb <- .lnorm_multifields(counts_gb, genebodies.gr, field)
        } else {
            counts_pr <- counts_pr / GenomicRanges::width(promoters.gr)
            counts_gb <- counts_gb / GenomicRanges::width(genebodies.gr)
        }
    }

    if (remove.empty) {
        if (length(field) > 1) {
            idx <- lapply(counts_pr, function(x) which(x != 0))
            idx <- Reduce(union, idx)
            counts_pr <- counts_pr[idx, ]
            counts_gb <- counts_gb[idx, ]
        } else {
            idx <- which(counts_pr != 0)
            counts_pr <- counts_pr[idx]
            counts_gb <- counts_gb[idx]
        }
    }

    return(counts_pr / counts_gb)
}

.lnorm_multifields <- function(counts, regions, field) {
    counts <- lapply(counts, "/", GenomicRanges::width(regions))
    counts <- as.data.frame(counts)
    names(counts) <- field
    return(counts)
}



