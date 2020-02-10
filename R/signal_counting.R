

#' Get signal counts in regions of interest
#'
#' Get the sum of the signal in \code{dataset.gr} that overlaps each range in
#' \code{regions.gr}.
#'
#' @param dataset.gr A GRanges object in which signal is contained in metadata
#'   (typically in the "score" field), or a named list of such GRanges objects.
#'   If a list is given, a dataframe is returned containing the counts in each
#'   region for each dataset.
#' @param regions.gr A GRanges object containing regions of interest.
#' @param field The metadata field of \code{dataset.gr} to be counted. If
#'   \code{length(field) > 1}, a dataframe is returned containing the counts for
#'   each region in each field. If \code{field} not found in
#'   \code{names(mcols(dataset.gr))}, will default to using all fields found in
#'   \code{dataset.gr}.
#' @param NF An optional normalization factor by which to multiply the counts.
#'   If given, \code{length(NF)} must be equal to \code{length(field)}.
#' @param blacklist An optional GRanges object containing regions that should be
#'   excluded from signal counting.
#' @param ncores Multiple cores will only be used if \code{dataset.gr} is a list
#'   of multiple datasets, or if \code{length(field) > 1}.
#'
#' @return An atomic vector the same length as \code{regions.gr} containing the
#'   sum of the signal overlapping each range of \code{regions.gr}. If
#'   \code{dataset.gr} is a list of multiple GRanges, or if \code{length(field)
#'   > 1}, a dataframe is returned.
#'
#' @author Mike DeBerardine
#' @seealso \code{\link[BRGenomics:getCountsByPositions]{getCountsByPositions}}
#' @export
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
                               NF = NULL, blacklist = NULL,
                               ncores = detectCores()) {



    if (!is.null(blacklist))
        dataset.gr <- .blacklist(dataset.gr, blacklist, ncores)

    if (is.list(dataset.gr)) {
        NF <- .check_nfs(dataset.gr, NF, field)
        cl <- mcMap(.get_cbr, dataset.gr, list(regions.gr), field, NF,
                    mc.cores = ncores)
        return(as.data.frame(cl))
    }

    # convenience: if field not in dataset.gr, use all fields that are
    field <- .check_fields(dataset.gr, field)
    NF <- .check_nfs(dataset.gr, NF, field)

    if (length(field) > 1) {
        cl <- mcMap(.get_cbr, list(dataset.gr), list(regions.gr), field, NF)
        names(cl) <- field
        as.data.frame(cl)
    } else {
        .get_cbr(dataset.gr, regions.gr, field, NF)
    }
}


#' @importFrom parallel mclapply
.blacklist <- function(dataset.gr, blacklist, ncores) {
    if (is.list(dataset.gr)) {
        mclapply(dataset.gr, subsetByOverlaps, blacklist, invert = TRUE,
                 mc.cores = ncores)
    } else {
        subsetByOverlaps(dataset.gr, blacklist, invert = TRUE)
    }
}


.check_nfs <- function(dataset.gr, NF, field) {
    if (is.list(dataset.gr)) {
        n <- length(dataset.gr)
    } else {
        n <- length(field)
    }
    if (is.null(NF))  NF <- rep(1L, n)

    if (length(NF) != n) {
        if (is.list(dataset.gr)) {
            msg <- .nicemsg("If dataset.gr is a list and an NF is given, then
                            length(NF) must equal length(dataset.gr)")
        } else {
            msg <- .nicemsg("If NF is given and length(field) > 1, then
                            length(NF) must equal length(field)")
        }
        stop(msg)
        return(geterrmessage())
    }
    return(NF)
}


#' @importFrom GenomicRanges findOverlaps mcols
.get_cbr <- function(dataset.gr, regions.gr, field, NF) {
    hits <- findOverlaps(regions.gr, dataset.gr)
    counts <- aggregate(mcols(dataset.gr)[[field]][hits@to],
                        by = list(hits@from), FUN = sum)
    names(counts) <- c("idx", "signal")
    counts.all <- rep(0L, length(regions.gr)) # include regions without hits
    counts.all[counts$idx] <- counts$signal
    counts.all * NF
}


#' @importFrom GenomicRanges mcols
.check_fields <- function(dataset.gr, field) {
    fnames <- names(mcols(dataset.gr))
    if (!all(field %in% fnames)) {
        message(.nicemsg("field not found in dataset.gr; will default to using
                         all fields in dataset.gr"))
        field <- fnames
    }
    return(field)
}


#' Get signal counts at each position within regions of interest
#'
#' Get the sum of the signal in \code{dataset.gr} that overlaps each position
#' within each range in \code{regions.gr}. If binning is used (i.e. positions
#' are wider than 1 bp), any function can be used to summarize the signal
#' overlapping each bin.
#'
#' @param dataset.gr A GRanges object in which signal is contained in metadata
#'   (typically in the "score" field), or a named list of such GRanges objects.
#' @param regions.gr A GRanges object containing regions of interest.
#' @param binsize Size of bins (in bp) to use for counting within each range of
#'   \code{regions.gr}. Note that counts will \emph{not} be length-normalized.
#' @param FUN If \code{binsize > 1}, the function used to aggregate the signal
#'   within each bin. By default, the signal is summed, but any function
#'   operating on a numeric vector can be used.
#' @param simplify.multi.widths A string indicating the output format if the
#'   ranges in \code{regions.gr} have variable widths. By default, an error is
#'   returned. See details below.
#' @param field The metadata field of \code{dataset.gr} to be counted. If
#'   \code{length(field) > 1}, the output is a list whose elements contain the
#'   output for generated each field. If \code{field} not found in
#'   \code{names(mcols(dataset.gr))}, will default to using all fields found in
#'   \code{dataset.gr}.
#' @param NF An optional normalization factor by which to multiply the counts.
#'   If given, \code{length(NF)} must be equal to \code{length(field)}.
#' @param melt A logical indicating if the count matrices should be melted. If
#'   set to \code{TRUE}, a dataframe is returned in containing columns for
#'   "region", "position", and "signal". If \code{dataset.gr} is a list of
#'   multiple GRanges, or if \code{length(field) > 1}, a single dataframe is
#'   returned, which contains an additional column "sample", which contains
#'   individual sample names. If used with multi-width \code{regions.gr}, the
#'   resulting dataframe will only contain positions that are found within each
#'   respective region.
#' @param blacklist An optional GRanges object containing regions that should be
#'   excluded from signal counting.
#' @param NA_blacklisted A logical indicating if NA values should be returned
#'   for blacklisted regions. By default, signal in the blacklisted sites is
#'   ignored, i.e. the reads are excluded. If \code{NA_blacklisted = TRUE},
#'   those positions are set to \code{NA} in the final output.
#' @param ncores Multiple cores will only be used if \code{dataset.gr} is a list
#'   of multiple datasets, or if \code{length(field) > 1}.
#'
#' @return If the widths of all ranges in \code{regions.gr} are equal, a matrix
#'   is returned that contains a row for each region of interest, and a column
#'   for each position (each base if \code{binsize = 1}) within each region. If
#'   \code{dataset.gr} is a list, a parallel list is returned containing a
#'   matrix for each input dataset.
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
#' @importFrom parallel detectCores mcMap
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
#' round(countsmat[10:15, ], 1)
getCountsByPositions <- function(dataset.gr, regions.gr, binsize = 1, FUN = sum,
                                 simplify.multi.widths = c("error", "list",
                                                           "pad 0", "pad NA"),
                                 field = "score", NF = NULL, melt = FALSE,
                                 blacklist = NULL, NA_blacklisted = FALSE,
                                 ncores = detectCores()) {

    smw <- match.arg(simplify.multi.widths, c("error", "list",
                                              "pad 0", "pad NA"))

    # apply blacklisting; xy positions (matrix row, column) of blacklist hits;
    # delay application if multiwidth -> send down arg as xy.bl
    if (smw == "error") {
        xy.bl <- NULL
    } else {
        xy.bl <- NA_blacklisted
    }

    if (!is.null(blacklist)) {
        dataset.gr <- .blacklist(dataset.gr, blacklist, ncores)
        if (NA_blacklisted & smw == "error") {
            blacklist <- GPos(reduce(blacklist))
            hits.bl <- findOverlaps(regions.gr, blacklist)
            xy.bl <- .get_positions_in_regions(hits.bl, blacklist, regions.gr)
        }
    }

    # function dispatch (for lists)
    if (is.list(dataset.gr)) {
        NF <- .check_nfs(dataset.gr, NF, field)
        cl <- mcMap(.get_cbp, dataset.gr, list(regions.gr), binsize, list(FUN),
                    smw, field, NF, melt, list(blacklist), list(xy.bl),
                    ncores = 1, mc.cores = ncores)
        if (melt)  cl <- .dfList2df(cl, prepend = FALSE)
        return(cl)

    } else {
        # convenience: if field not in dataset.gr, use all fields that are
        field <- .check_fields(dataset.gr, field)
        NF <- .check_nfs(dataset.gr, NF, field)
        .get_cbp(dataset.gr, regions.gr, binsize, FUN, smw, field, NF, melt,
                 blacklist, xy.bl, ncores)
    }

}


.get_cbp <- function(dataset.gr, regions.gr, binsize, FUN, smw, field, NF, melt,
                     blacklist, xy.bl, ncores) {

    # function makes practical use of single-width dataset.gr, but the output
    # is always valid regardless of the input data type
    if (any(width(dataset.gr) != 1)) dataset.gr <- makeGRangesBRG(dataset.gr)

    # get hits early (so happens once for if multiple fields)
    hits <- findOverlaps(regions.gr, dataset.gr)

    # function dispatch
    if (length(field) > 1) {
        call_multifield <- function(field.i, NF.i) {
            .split_cbp(hits, dataset.gr, regions.gr, binsize = binsize,
                       FUN = FUN, smw = smw, field = field.i, NF = NF.i,
                       melt = melt, blacklist = blacklist, xy.bl = xy.bl)
        }
        cl <- mcMap(call_multifield, field, NF, mc.cores = ncores)
        names(cl) <- field
        if (melt)  cl <- .dfList2df(cl, prepend = FALSE)
        return(cl)

    } else {
        .split_cbp(hits, dataset.gr, regions.gr, binsize = binsize, FUN = FUN,
                   smw, field = field, NF = NF, melt = melt,
                   blacklist = blacklist, xy.bl = xy.bl)
    }
}

.split_cbp <- function(hits, dataset.gr, regions.gr, binsize, FUN, smw, field,
                       NF, melt, blacklist, xy.bl) {

    multi_width <- length(unique(width(regions.gr))) > 1

    if (multi_width) {
        .get_cbp_mw(hits, dataset.gr, regions.gr, binsize, FUN, smw, field, NF,
                    melt, blacklist, xy.bl)
    } else {
        .check_mw_arg(smw, melt, mw = FALSE)
        .get_signal_mat(hits, dataset.gr, regions.gr, binsize, FUN, field, NF,
                        melt, blacklist, xy.bl)
    }

}


.check_mw_arg <- function(smw, melt, mw) {
    if (mw) { # if dispatched from .get_cbp_bw
        if (smw == "error") {
            stop(message = .nicemsg("regions.gr contains ranges with multiple
                                    widths, but simplify.multi.widths is set to
                                    'error'. Did you mean to call
                                    getCountsByRegions instead?"))
            return(geterrmessage())
        }
        if (melt)  smw <- "list"
        return(smw)

    } else {
        # check is necessary because NA_blacklisting would fail otherwise, but
        # will give error in every case to not confuse users if they discover other
        # conditions when the error isn't returned
        if (smw != "error") {
            stop(message = .nicemsg("simplify.multi.widths changed from default,
                                    but regions.gr is not multiwidth"))
            return(geterrmessage())
        }
    }

}


#' @importFrom GenomicRanges width resize
.get_cbp_mw <- function(hits, dataset.gr, regions.gr, binsize, FUN, smw, field,
                        NF, melt, blacklist, xy.bl) {

    smw <- .check_mw_arg(smw, melt, mw = TRUE)

    # expand all regions to be the same width, then get counts
    widths <- width(regions.gr) # save widths
    suppressWarnings( regions.gr <- resize(regions.gr, max(widths)) )

    # check if NA_blacklisted (which for multiwidth is xy.bl)
    if (xy.bl) {
        blacklist <- GPos(reduce(blacklist))
        hits.bl <- findOverlaps(regions.gr, blacklist)
        xy.bl <- .get_positions_in_regions(hits.bl, blacklist, regions.gr)
    } else {
        xy.bl <- NULL
    }

    # get counts matrix
    mat <- .get_signal_mat(hits, dataset.gr, regions.gr, binsize, FUN, field,
                           NF, melt = FALSE, blacklist, xy.bl)
    nbins <- floor(widths / binsize) # number of bins to keep for each region

    if (smw == "list") { # includes if melt == TRUE
        # list of vectors whose lengths determined by width of range i
        cl <- Map(function(row.i, nbins.i) mat[row.i, seq_len(nbins.i)],
                  seq_len(nrow(mat)), nbins)
        if (melt)  cl <- .meltmw(cl)
        return(cl)

    } else {
        # get array indices for out-of-range bins
        arridx_pad <- vapply(nbins, function(n) seq_len(ncol(mat)) > n,
                             FUN.VALUE = logical(ncol(mat)))
        # (transpose as sapply/vapply cbinds the rows)
        arridx_pad <- which( t(arridx_pad), arr.ind = TRUE )
        mat[arridx_pad] <- ifelse(smw == "pad 0", 0, NA)
        return(mat)
    }
}


.meltmw <- function(cl) {
    lens <- lengths(cl)
    data.frame(region = rep(seq_along(cl), lens),
               position = unlist(lapply(lens, seq_len)),
               signal = unlist(cl))
}


#' @importFrom GenomicRanges start mcols
.get_signal_mat <- function(hits, dataset.gr, regions.gr, binsize, FUN, field,
                            NF, melt, blacklist, xy.bl) {

    # initialize signal matrix of dim = (region, position within region)
    rwidth <- width(regions.gr[1])
    if (length(rwidth) == 0) {
        stop(.nicemsg("Cannot make counts matrix because regions.gr is empty"))
        return(geterrmessage())
    }
    mat <- matrix(0L, length(regions.gr), rwidth)

    # find (x, y) = (region, position)
    xy <- .get_positions_in_regions(hits, dataset.gr, regions.gr)
    mat[xy] <- mcols(dataset.gr)[[field]][hits@to]

    if (!is.null(xy.bl))  mat[xy.bl] <- NA # apply NA_blacklisting

    if (binsize > 1) {
        mat <- apply(mat, 1, function(x) .binVector(x, binsize = binsize,
                                                    FUN = match.fun(FUN)))
        mat <- t(mat) # apply will cbind rather than rbind
    }
    mat <- mat * NF # apply normalization
    if (melt)  mat <- .meltmat(mat)
    return(mat)
}

.get_positions_in_regions <- function(hits, dataset.gr, regions.gr) {
    x <- hits@from
    y1 <- start(dataset.gr[hits@to]) # site of signal
    y2 <- start(resize(regions.gr[hits@from], 1)) # beginning of window
    y <- abs(y1 - y2) + 1 # position of signal within region
    cbind(x, y)
}


.meltmat <- function(mat) {
    mat <- t(mat) # (to sort by region, rather than position)
    df <- expand.grid( seq_len(nrow(mat)), seq_len(ncol(mat)) )
    df <- cbind(df, data.frame(as.vector(mat)))
    names(df) <- c("position", "region", "signal")
    df[, c(2, 1, 3)]
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
#'   (typically in the "score" field), or a named list of such GRanges objects.
#' @param promoters.gr A GRanges object containing promoter-proximal regions of
#'   interest.
#' @param genebodies.gr A GRanges object containing genebody regions of
#'   interest.
#' @param field The metadata field of \code{dataset.gr} to be counted. If
#'   \code{length(field) > 1}, a dataframe is returned containing the pausing
#'   indices for each region in each field. If \code{field} not found in
#'   \code{names(mcols(dataset.gr))}, will default to using all fields found in
#'   \code{dataset.gr}. If \code{dataset.gr} is a list, a single \code{field}
#'   should be given, or \code{length(field)} should be the equal to the number
#'   of datasets in \code{dataset.gr}.
#' @param length.normalize A logical indicating if signal counts within regions
#'   of interest should be length normalized. The default is \code{TRUE}, which
#'   is recommended, especially if input regions don't all have the same width.
#' @param remove.empty A logical indicating if genes without any signal in
#'   \code{promoters.gr} should be removed. No genes are filtered by default. If
#'   \code{dataset.gr} is a list of datasets, or if \code{length(field) > 1},
#'   regions are filtered unless they have promoter signal in all datasets.
#' @param blacklist An optional GRanges object containing regions that should be
#'   excluded from signal counting.
#' @param ncores Multiple cores will only be used if \code{dataset.gr} is a list
#'   of multiple datasets, or if \code{length(field) > 1}.
#'
#' @return A vector parallel to the input genelist, unless \code{remove.empty =
#'   TRUE}, in which case the vector may be shorter. If \code{length(field) >
#'   1}, a dataframe is returned, containing a column for each field.
#'
#' @author Mike DeBerardine
#' @seealso \code{\link[BRGenomics:getCountsByRegions]{getCountsByRegions}}
#' @export
#' @importFrom parallel detectCores mcMap
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
                              remove.empty = FALSE, blacklist = NULL,
                              ncores = detectCores()) {

    if (length(promoters.gr) != length(genebodies.gr)) {
        stop(message = .nicemsg("Number of ranges in promoters.gr not equal to
                                number of ranges in genebodies.gr"))
        return(geterrmessage())
    }

    if (!is.null(blacklist))
        dataset.gr <- .blacklist(dataset.gr, blacklist, ncores)

    # vectors for single dataset; dataframes for lists or multiple fields
    counts_pr <- getCountsByRegions(dataset.gr, promoters.gr, field = field,
                                    ncores = ncores)
    counts_gb <- getCountsByRegions(dataset.gr, genebodies.gr, field = field,
                                    ncores = ncores)

    if (is.list(dataset.gr)) {
        .pidx_multi(counts_pr, counts_gb, promoters.gr, genebodies.gr,
                    names(dataset.gr), length.normalize, remove.empty)
    } else if (length(field) > 1) {
        .pidx_multi(counts_pr, counts_gb, promoters.gr, genebodies.gr,
                    field, length.normalize, remove.empty)
    } else {
        .pidx_single(counts_pr, counts_gb, promoters.gr, genebodies.gr,
                     length.normalize, remove.empty)
    }
}


.pidx_single <- function(counts_pr, counts_gb, promoters.gr, genebodies.gr,
                         length.normalize, remove.empty) {

    if (length.normalize) {
        counts_pr <- counts_pr / GenomicRanges::width(promoters.gr)
        counts_gb <- counts_gb / GenomicRanges::width(genebodies.gr)
    }

    if (remove.empty) {
        idx <- which(counts_pr != 0)
        counts_pr <- counts_pr[idx]
        counts_gb <- counts_gb[idx]
    }

    return(counts_pr / counts_gb)
}


.pidx_multi <- function(counts_pr, counts_gb, promoters.gr, genebodies.gr,
                        dnames, length.normalize, remove.empty) {

    if (length.normalize) {
        counts_pr <- .lnorm_multi(counts_pr, promoters.gr, dnames)
        counts_gb <- .lnorm_multi(counts_gb, genebodies.gr, dnames)
    }

    if (remove.empty) {
        # idx to drop
        idx <- lapply(counts_pr, function(x) which(x == 0))
        idx <- unique(unlist(idx))
        counts_pr <- counts_pr[-idx, ]
        counts_gb <- counts_gb[-idx, ]
    }

    return(counts_pr / counts_gb)
}


.lnorm_multi <- function(counts, regions, dnames) {
    counts <- lapply(counts, "/", GenomicRanges::width(regions))
    counts <- as.data.frame(counts)
    names(counts) <- dnames
    return(counts)
}



