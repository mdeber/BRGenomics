

#' Get signal counts in regions of interest
#'
#' Get the sum of the signal in \code{dataset.gr} that overlaps each range in
#' \code{regions.gr}. This function is written to calculate \emph{readcounts}
#' overlapping each region, and \emph{not} "coverage signal" (see details
#' below).
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
#' @param melt If \code{melt = TRUE}, a dataframe is returned containing a
#'   column for regions and another column for signal. If multiple datasets are
#'   given (if \code{dataset.gr} is a list or if \code{length(field) > 1}), the
#'   output dataframe is melted to contain a third column indicating the sample
#'   names. (See section on return values below).
#' @param region_names If \code{melt = TRUE}, an optional vector of names for
#'   the regions in \code{regions.gr}. If left as \code{NULL}, indices of
#'   \code{regions.gr} are used instead.
#' @param ncores Multiple cores will only be used if \code{dataset.gr} is a list
#'   of multiple datasets, or if \code{length(field) > 1}.
#'
#' @return An atomic vector the same length as \code{regions.gr} containing the
#'   sum of the signal overlapping each range of \code{regions.gr}. If
#'   \code{dataset.gr} is a list of multiple GRanges, or if \code{length(field)
#'   > 1}, a dataframe is returned. If \code{melt = FALSE} (the default),
#'   dataframes have a column for each dataset and a row for each region. If
#'   \code{melt = TRUE}, dataframes contain one column to indicate regions
#'   (either by their indices, or by \code{region_names}, if given), another
#'   column to indicate signal, and a third column containing the sample name
#'   (unless \code{dataset.gr} is a single GRanges object).
#'
#' @details This function is designed to work with data in which each range
#'   represents one type of molecule, whether it's a single base (e.g. the 5'
#'   ends, 3' ends, or centers of reads) or entire reads (i.e. paired 5' and 3'
#'   ends of reads).
#'
#'   This is in contrast to standard run-length compressed GRanges object, as
#'   imported using \code{\link[rtracklayer:import.bw]{rtracklayer::import.bw}},
#'   in which a single range can represent multiple contiguous positions that
#'   share the same signal information. This function does \emph{not} expand
#'   run-length compressed coverage objects. As an example, a range of length 10
#'   with a score of 2 is treated as 2 reads (each spanning the same 10 bases),
#'   not 20 reads.
#'
#' @author Mike DeBerardine
#' @seealso \code{\link[BRGenomics:getCountsByPositions]{getCountsByPositions}}
#' @export
#' @importFrom parallel detectCores mcMap mclapply
#' @importFrom methods is
#' @importFrom GenomicRanges findOverlaps mcols
#' @importFrom IRanges subsetByOverlaps
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
                               NF = NULL, blacklist = NULL, melt = FALSE,
                               region_names = NULL, ncores = detectCores()) {

    if (!is.null(blacklist))
        dataset.gr <- .blacklist(dataset.gr, blacklist, ncores)

    if (is.list(dataset.gr) || is(dataset.gr, "GRangesList")) {
        NF <- .check_nfs(dataset.gr, NF, field)
        # subset to increase performance for large datasets
        dataset.gr <- mclapply(dataset.gr, subsetByOverlaps, regions.gr,
                               mc.cores = ncores)
        cl <- mcMap(.get_cbr, dataset.gr, list(regions.gr), field, NF,
                    mc.cores = ncores)
        cl <- as.data.frame(cl)
        if (melt) return(.melt_counts(cl, colnames(cl), region_names))
        return(cl)
    }

    # convenience: if field not in dataset.gr, use all fields that are
    field <- .check_fields(dataset.gr, field)
    NF <- .check_nfs(dataset.gr, NF, field)

    # subset to increase performance for large datasets
    dataset.gr <- subsetByOverlaps(dataset.gr, regions.gr)

    if (length(field) > 1) {
        cl <- mcMap(.get_cbr, list(dataset.gr), list(regions.gr), field, NF,
                    mc.cores = ncores)
        names(cl) <- field
        cl <- as.data.frame(cl)
        if (melt) return(.melt_counts(cl, field, region_names))
        return(cl)
    } else {
        counts <- .get_cbr(dataset.gr, regions.gr, field, NF)
        if (melt) return(.melt_counts(counts, snames = NULL, region_names))
        return(counts)
    }
}

#' @importFrom methods is
#' @importFrom parallel mclapply
#' @importFrom IRanges subsetByOverlaps
.blacklist <- function(dataset.gr, blacklist, ncores) {
    if (is.list(dataset.gr) || is(dataset.gr, "GRangesList")) {
        mclapply(dataset.gr, subsetByOverlaps, blacklist, invert = TRUE,
                 mc.cores = ncores)
    } else {
        subsetByOverlaps(dataset.gr, blacklist, invert = TRUE)
    }
}

#' @importFrom methods is
.check_nfs <- function(dataset.gr, NF, field) {
    if (is.list(dataset.gr) || is(dataset.gr, "GRangesList")) {
        n <- length(dataset.gr)
    } else {
        n <- length(field)
    }

    if (is.null(NF))  NF <- rep(1L, n)

    if (length(NF) != n) {
        if (is.list(dataset.gr))
            stop(.nicemsg("If dataset.gr is a list and an NF is given, then
                          length(NF) must equal length(dataset.gr)"))
        stop(.nicemsg("If NF is given and length(field) > 1, then length(NF)
                      must equal length(field)"))
    }
    NF
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
    field
}


.melt_counts <- function(df, snames, region_names) {
    if (is.vector(df))
        df <- data.frame(df)
    nr <- nrow(df) # number of regions
    ns <- ncol(df) # number of samples

    if (is.null(region_names))
        region_names <- seq_len(nr)

    if (length(snames) > 1) {
        df <- data.frame(region = rep(region_names, ns),
                         signal = unlist(df),
                         sample = rep(snames, each = nr))
    } else {
        df <- data.frame(region = region_names,
                         signal = unlist(df))
    }
    rownames(df) <- NULL
    df
}


#' Get signal counts at each position within regions of interest
#'
#' Get the sum of the signal in \code{dataset.gr} that overlaps each position
#' within each range in \code{regions.gr}. If binning is used (i.e. positions
#' are wider than 1 bp), any function can be used to summarize the signal
#' overlapping each bin. Just as in \code{\link[BRGenomics:getCountsByRegions]{
#' getCountsByRegions}}, this function calculates \emph{readcounts} overlapping
#' positions, and not "coverage signal".
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
#' @param blacklist An optional GRanges object containing regions that should be
#'   excluded from signal counting.
#' @param NA_blacklisted A logical indicating if NA values should be returned
#'   for blacklisted regions. By default, signal in the blacklisted sites is
#'   ignored, i.e. the reads are excluded. If \code{NA_blacklisted = TRUE},
#'   those positions are set to \code{NA} in the final output.
#' @param melt A logical indicating if the count matrices should be melted. If
#'   set to \code{TRUE}, a dataframe is returned in containing columns for
#'   "region", "position", and "signal". If \code{dataset.gr} is a list of
#'   multiple GRanges, or if \code{length(field) > 1}, a single dataframe is
#'   returned, which contains an additional column "sample", which contains
#'   individual sample names. If used with multi-width \code{regions.gr}, the
#'   resulting dataframe will only contain positions that are found within each
#'   respective region.
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
#' @importFrom parallel detectCores mcMap mclapply
#' @importFrom methods is
#' @importFrom GenomicRanges width
#' @importFrom IRanges subsetByOverlaps
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
                                 field = "score", NF = NULL, blacklist = NULL,
                                 NA_blacklisted = FALSE, melt = FALSE,
                                 ncores = detectCores()) {

    smw <- match.arg(simplify.multi.widths, c("error", "list",
                                              "pad 0", "pad NA"))

    # apply blacklisting; xy positions (matrix row, column) of blacklist hits;
    # delay application if multiwidth -> send down arg as xy.bl
    xy.bl <- if (smw == "error") NULL else NA_blacklisted

    if (!is.null(blacklist)) {
        dataset.gr <- .blacklist(dataset.gr, blacklist, ncores)
        if (NA_blacklisted & smw == "error") {
            blacklist <- GPos(reduce(blacklist))
            hits.bl <- findOverlaps(regions.gr, blacklist)
            xy.bl <- .get_positions_in_regions(hits.bl, blacklist, regions.gr)
        }
    }

    # function dispatch (for lists)
    if (is.list(dataset.gr) || is(dataset.gr, "GRangesList")) {
        NF <- .check_nfs(dataset.gr, NF, field)
        # subset to increase performance for large datasets
        dataset.gr <- mclapply(dataset.gr, subsetByOverlaps, regions.gr,
                               mc.cores = ncores)
        cl <- mcMap(.get_cbp, dataset.gr, list(regions.gr), binsize, list(FUN),
                    smw, field, NF, melt, list(blacklist), list(xy.bl),
                    ncores = 1, mc.cores = ncores)
        if (melt)  cl <- .dfList2df(cl, prepend = FALSE)
        return(cl)

    } else {
        # convenience: if field not in dataset.gr, use all fields that are
        field <- .check_fields(dataset.gr, field)
        NF <- .check_nfs(dataset.gr, NF, field)
        # increase performance for large datasets
        dataset.gr <- subsetByOverlaps(dataset.gr, regions.gr)
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
        if (smw != "error")
            # check is necessary because NA_blacklisting would fail otherwise,
            # but will give an error in every case for the sake of consistency
            stop(.nicemsg("simplify.multi.widths changed from default, but
                          regions.gr is not multiwidth"))
        .get_signal_mat(hits, dataset.gr, regions.gr, binsize, FUN, field, NF,
                        melt, blacklist, xy.bl)
    }
}



#' @importFrom GenomicRanges width resize
.get_cbp_mw <- function(hits, dataset.gr, regions.gr, binsize, FUN, smw, field,
                        NF, melt, blacklist, xy.bl) {
    if (smw == "error")
        stop(.nicemsg("regions.gr contains ranges with multiple widths, but
                      simplify.multi.widths is set to 'error'. Did you mean
                      to call getCountsByRegions instead?"))

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

    if (melt | smw == "list") {
        # list of vectors whose lengths determined by width of range i
        cl <- Map(function(row.i, nbins.i) mat[row.i, seq_len(nbins.i)],
                  seq_len(nrow(mat)), nbins)
        if (melt)  return(.meltmw(cl))
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
    df <- data.frame(region = rep(seq_along(cl), lens),
                     position = unlist(lapply(lens, seq_len)),
                     signal = unlist(cl))
    rownames(df) <- NULL
    df
}


#' @importFrom GenomicRanges start mcols
.get_signal_mat <- function(hits, dataset.gr, regions.gr, binsize, FUN, field,
                            NF, melt, blacklist, xy.bl) {

    # initialize signal matrix of dim = (region, position within region)
    rwidth <- width(regions.gr[1])
    if (length(rwidth) == 0)
        stop("Cannot make counts matrix because regions.gr is empty")

    mat <- matrix(0L, length(regions.gr), rwidth)

    # find (x, y) = (region, position)
    xy <- .get_positions_in_regions(hits, dataset.gr, regions.gr)
    mat[xy] <- mcols(dataset.gr)[[field]][hits@to]

    if (!is.null(xy.bl)) # apply NA_blacklisting
        mat[xy.bl] <- NA

    if (binsize > 1) {
        mat <- apply(mat, 1, function(x) .binVector(x, binsize = binsize,
                                                    FUN = match.fun(FUN)))
        mat <- t(mat) # apply will cbind rather than rbind
    }
    mat <- mat * NF # apply normalization
    if (melt) return(.meltmat(mat))
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
    rownames(df) <- NULL
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
#'   excluded from signal counting. If \code{length.normalize = TRUE},
#'   blacklisted positions will be excluded from length calculations. Users
#'   should take care to note if regions of interest substantially overlap
#'   blacklisted positions.
#' @param melt If \code{melt = TRUE}, a dataframe is returned containing a
#'   column for regions and another column for pausing indices. If multiple
#'   datasets are given (if \code{dataset.gr} is a list or if
#'   \code{length(field) > 1}), the output dataframe is melted to contain a
#'   third column indicating the sample names. (See section on return values
#'   below).
#' @param region_names If \code{melt = TRUE}, an optional vector of names for
#'   the regions in \code{regions.gr}. If left as \code{NULL}, indices of
#'   \code{regions.gr} are used instead.
#' @param ncores Multiple cores will only be used if \code{dataset.gr} is a list
#'   of multiple datasets, or if \code{length(field) > 1}.
#'
#' @return A vector parallel to the input genelist, unless \code{remove.empty =
#'   TRUE}, in which case the vector may be shorter. If \code{dataset.gr} is a
#'   list, or if \code{length(field) > 1}, a dataframe is returned, containing a
#'   column for each field. However, if \code{melt = TRUE}, dataframes contain
#'   one column to indicate regions (either by their indices, or by
#'   \code{region_names}, if given), another column to indicate signal, and a
#'   third column containing the sample name (unless \code{dataset.gr} is a
#'   single GRanges object).
#'
#' @author Mike DeBerardine
#' @seealso \code{\link[BRGenomics:getCountsByRegions]{getCountsByRegions}}
#' @export
#' @importFrom parallel detectCores mcMap
#' @importFrom methods is
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
                              melt = FALSE, region_names = NULL,
                              ncores = detectCores()) {

    if (length(promoters.gr) != length(genebodies.gr))
        stop(message = .nicemsg("Number of ranges in promoters.gr not equal to
                                number of ranges in genebodies.gr"))

    dwidth <- NULL # differences in width caused by blacklisting
    if (!is.null(blacklist)) {
        dataset.gr <- .blacklist(dataset.gr, blacklist, ncores)
        if (length.normalize) {
            dwidth <- mclapply(list(promoters.gr, genebodies.gr),
                               .get_dwidth, blacklist, mc.cores = ncores)
        }
    }

    # vectors for single dataset; dataframes for lists or multiple fields
    counts_pr <- getCountsByRegions(dataset.gr, promoters.gr, field = field,
                                    ncores = ncores)
    counts_gb <- getCountsByRegions(dataset.gr, genebodies.gr, field = field,
                                    ncores = ncores)

    if (is.list(dataset.gr) || is(dataset.gr, "GRangesList")) {
        .pidx_multi(counts_pr, counts_gb, promoters.gr, genebodies.gr,
                    names(dataset.gr), length.normalize, remove.empty, melt,
                    region_names, dwidth)

    } else if (length(field) > 1) {
        .pidx_multi(counts_pr, counts_gb, promoters.gr, genebodies.gr,
                    field, length.normalize, remove.empty, melt, region_names,
                    dwidth)

    } else {
        .pidx_single(counts_pr, counts_gb, promoters.gr, genebodies.gr,
                     length.normalize, remove.empty, melt, region_names, dwidth)
    }
}


.get_dwidth <- function(regions, blacklist) {
    dwidth <- rep(0L, length(regions))
    hits <- findOverlaps(regions, blacklist)
    bloverlap <- pintersect(regions[hits@from], blacklist[hits@to])
    dwidth[hits@from] <- width(bloverlap)
    dwidth
}

.pidx_single <- function(counts_pr, counts_gb, promoters.gr, genebodies.gr,
                         length.normalize, remove.empty, melt, region_names,
                         dwidth) {

    if (length.normalize) {
        pwidths <- GenomicRanges::width(promoters.gr)
        gwidths <- GenomicRanges::width(genebodies.gr)
        if (!is.null(dwidth)) {
            pwidths <- pwidths - dwidth[[1]]
            gwidths <- gwidths - dwidth[[2]]
        }
        counts_pr <- counts_pr / pwidths
        counts_gb <- counts_gb / gwidths
    }

    if (remove.empty) {
        idx <- which(counts_pr != 0)
        counts_pr <- counts_pr[idx]
        counts_gb <- counts_gb[idx]
    }

    pidx <- counts_pr / counts_gb

    if (melt) {
        if (remove.empty) {
            if (is.null(region_names)) {
                region_names <- idx
            } else {
                region_names <- region_names[idx]
            }
        }
        pidx <- .melt_counts(pidx, NULL, region_names)
        colnames(pidx) <- c("region", "pauseIndex")
    }
    pidx
}

.pidx_multi <- function(counts_pr, counts_gb, promoters.gr, genebodies.gr,
                        dnames, length.normalize, remove.empty, melt,
                        region_names, dwidth) {

    if (length.normalize) {
        if (is.null(dwidth))
            dwidth <- list(NULL, NULL)
        counts_pr <- .lnorm_multi(counts_pr, promoters.gr, dnames, dwidth[[1]])
        counts_gb <- .lnorm_multi(counts_gb, genebodies.gr, dnames, dwidth[[2]])
    }

    if (remove.empty) {
        # idx to drop
        idx <- lapply(counts_pr, function(x) which(x == 0))
        idx <- unique(unlist(idx))
        counts_pr <- counts_pr[-idx, ]
        counts_gb <- counts_gb[-idx, ]
    }

    pidx <- counts_pr / counts_gb

    if (melt) {
        if (remove.empty) {
            if (is.null(region_names)) {
                region_names <- seq_along(promoters.gr)[-idx]
            } else {
                region_names <- region_names[-idx]
            }
        }
        pidx <- .melt_counts(pidx, colnames(pidx), region_names)
        colnames(pidx) <- c("region", "pauseIndex", "sample")
    }
    pidx
}

#' @importFrom GenomicRanges width
.lnorm_multi <- function(counts, regions, dnames, dwidth.i) {
    widths.i <- width(regions)
    if (!is.null(dwidth.i))
        widths.i <- widths.i - dwidth.i
    counts <- as.data.frame(lapply(counts, "/", widths.i))
    names(counts) <- dnames
    counts
}
