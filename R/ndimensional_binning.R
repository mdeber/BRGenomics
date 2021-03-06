### ========================================================================= #
### N-dimensional binning by nbins
### ------------------------------------------------------------------------- #
###

#' Generating and Aggregating Data Within N-dimensional Bins
#'
#' Divide data along different dimensions into equally spaced bins, and
#' summarize the datapoints that fall into any of these n-dimensional bins.
#'
#' @param x The name of the dimension in \code{dims.df} to aggregate, or a
#'   separate numerical vector or dataframe of data to be aggregated. If
#'   \code{x} is a numerical vector, each value in \code{x} corresponds to a row
#'   of \code{dims.df}, and so \code{length(x)} must be equal to
#'   \code{nrow(dims.df)}. Likewise, if \code{x} is a dataframe, \code{nrow(x)}
#'   must equal \code{nrow(dims.df)}. Supplying a dataframe for \code{x} has the
#'   advantage of simultaneously aggregating different sets of data, and
#'   returning a single dataframe.
#' @param dims.df A dataframe containing one or more columns of numerical data
#'   for which bins will be generated.
#' @param nbins Either a number giving the number of bins to use for all
#'   dimensions (default = 10), or a vector containing the number of bins to use
#'   for each dimension of input data given.
#' @param use_bin_numbers A logical indicating if ordinal bin numbers should be
#'   returned (\code{TRUE}), or if in place of the bin number, the center value
#'   of that bin should be returned. For instance, if the first bin encompasses
#'   data from 1 to 3, with \code{use_bin_numbers = TRUE}, a 1 is returned, but
#'   when \code{FALSE}, 2 is returned.
#' @param FUN A function to use for aggregating data within each bin.
#' @param ... Additional arguments passed to \code{FUN}.
#' @param ignore.na Logical indicating if \code{NA} values of \code{x} should be
#'   ignored. Default is \code{TRUE}.
#' @param drop A logical indicating if empty bin combinations should be removed
#'   from the output. By default (\code{FALSE}), all possible combinations of
#'   bins are returned, and empty bins contain a value given by \code{empty}.
#' @param empty When \code{drop = FALSE}, the value returned for empty bins. By
#'   default, empty bins return \code{NA}. However, in many circumstances (e.g.
#'   if \code{FUN = sum}), the empty value should be \code{0}.
#' @param ncores Number of cores to use for computations.
#'
#' @return A dataframe.
#'
#' @details These functions take in data along 1 or more dimensions, and for
#'   each dimension the data is divided into evenly-sized bins from the minimum
#'   value to the maximum value. For instance, if each row of \code{dims.df}
#'   were a gene, the columns (the different dimensions) would be various
#'   quantitative measures of that gene, e.g. expression level, number of exons,
#'   length, etc. If plotted in cartesian coordinates, each gene would be a
#'   single datapoint, and each measurement would be a separate dimension.
#'
#'   \code{binNdimensions} returns the bin numbers themselves. The output
#'   dataframe has the same dimensions as the input \code{dims.df}, but each
#'   input data has been replaced by its bin number (an integer). If
#'   code{use_bin_numbers = FALSE}, the center points of the bins are returned
#'   instead of the bin numbers.
#'
#'   \code{aggregateByNdimBins} summarizes some input data \code{x} in each
#'   combination of bins, i.e. in each n-dimensional bin. Each row of the output
#'   dataframe is a unique combination of the input bins (i.e. each
#'   n-dimensional bin), and the output columns are identical to those in
#'   \code{dims.df}, with the addition of one or more columns containing the
#'   aggregated data in each n-dimensional bin. If the input \code{x} was a
#'   vector, the column is named "value"; if the input \code{x} was a dataframe,
#'   the column names from \code{x} are maintained.
#'
#'   \code{densityInNdimBins} returns a dataframe just like
#'   \code{aggregateByNdimBins}, except the "value" column contains the number
#'   of observations that fall into each n-dimensional bin.
#'
#' @author Mike DeBerardine
#' @export
#' @importFrom parallel mcMap mclapply
#' @examples
#' data("PROseq") # import included PROseq data
#' data("txs_dm6_chr4") # import included transcripts
#'
#' #--------------------------------------------------#
#' # find counts in promoter, early genebody, and near CPS
#' #--------------------------------------------------#
#'
#' pr <- promoters(txs_dm6_chr4, 0, 100)
#' early_gb <- genebodies(txs_dm6_chr4, 500, 1000, fix.end = "start")
#' cps <- genebodies(txs_dm6_chr4, -500, 500, fix.start = "end")
#'
#' df <- data.frame(counts_pr = getCountsByRegions(PROseq, pr),
#'                  counts_gb = getCountsByRegions(PROseq, early_gb),
#'                  counts_cps = getCountsByRegions(PROseq, cps))
#'
#' #--------------------------------------------------#
#' # divide genes into 20 bins for each measurement
#' #--------------------------------------------------#
#'
#' bin3d <- binNdimensions(df, nbins = 20, ncores = 1)
#'
#' length(txs_dm6_chr4)
#' nrow(bin3d)
#' bin3d[1:6, ]
#'
#' #--------------------------------------------------#
#' # get number of genes in each bin
#' #--------------------------------------------------#
#'
#' bin_counts <- densityInNdimBins(df, nbins = 20, ncores = 1)
#'
#' bin_counts[1:6, ]
#'
#' #--------------------------------------------------#
#' # get mean cps reads in bins of promoter and genebody reads
#' #--------------------------------------------------#
#'
#' bin2d_cps <- aggregateByNdimBins("counts_cps", df, nbins = 20,
#'                                  ncores = 1)
#'
#' bin2d_cps[1:6, ]
#'
#' subset(bin2d_cps, is.finite(counts_cps))[1:6, ]
#'
#' #--------------------------------------------------#
#' # get median cps reads for those bins
#' #--------------------------------------------------#
#'
#' bin2d_cps_med <- aggregateByNdimBins("counts_cps", df, nbins = 20,
#'                                      FUN = median, ncores = 1)
#'
#' bin2d_cps_med[1:6, ]
#'
#' subset(bin2d_cps_med, is.finite(counts_cps))[1:6, ]
binNdimensions <- function(dims.df, nbins = 10L, use_bin_numbers = TRUE,
                           ncores = getOption("mc.cores", 2L)) {

    # check input bins
    n_dim <- ncol(dims.df)
    if (length(nbins) > 1L & length(nbins) != n_dim)
        stop(.nicemsg("User input %d dimensions of data, but length(nbins) = %d.
                      length(nbins) must be equal to ncol(dims.df), or be a
                      single number.", n_dim, length(nbins)))

    # get bin divisions for each dimension, evenly spaced from min to max values
    getBreaks <- function(x, y) {
        xfin <- x[is.finite(x)]
        seq(min(xfin), max(xfin), length.out = y + 1L)
    }
    breaks <- mcMap(getBreaks, dims.df, nbins, mc.cores = ncores)

    # get bin indices for each datapoint along each dimension
    bin_idx <- mcMap(findInterval, dims.df, breaks, rightmost.closed = TRUE,
                     mc.cores = ncores)

    # address infinite values
    bin_idx <- lapply(bin_idx, function(x) {
        x[x == 0L] <- NA
        x
    })

    if (!use_bin_numbers) {
        # get center value for each bin in each dimensions
        get_centers <- function(x) vapply(seq_len(length(x) - 1L),
                                          function(i) mean(x[i:(i + 1L)]),
                                          numeric(1L))
        break_centers <- mclapply(breaks, get_centers, mc.cores = ncores)
        bin_idx <- mcMap("[", break_centers, bin_idx) # get values for data
    }

    names(bin_idx) <- paste0("bin.", names(dims.df))
    as.data.frame(bin_idx)
}



#' @rdname binNdimensions
#' @export
aggregateByNdimBins <- function(x, dims.df, nbins = 10L, FUN = mean, ...,
                                ignore.na = TRUE, drop = FALSE, empty = NA,
                                use_bin_numbers = TRUE,
                                ncores = getOption("mc.cores", 2L)) {
    if (is.character(x))  {
        colx <- which(names(dims.df) %in% x)
        x <- dims.df[colx]
        dims.df <- dims.df[-colx]
    }

    # if x is a dataframe, keep the original names
    if (is.vector(x)) {
        x.df <- data.frame(value = x)
    } else {
        x.df <- x
    }

    if (ignore.na) {
        idx <- apply(x.df, 1L, function(x) all(!is.na(x)))
        x.df <- x.df[idx, , drop = FALSE]
        dims.df <- dims.df[idx, , drop = FALSE]
    }

    # get bins for data
    bins.df <- binNdimensions(dims.df, nbins, use_bin_numbers, ncores)

    # only for this combination of arguments do we need to differentiate the
    # NAs that are returned by FUN from the NAs that result from empty bins
    if (!drop && !is.na(empty) && !ignore.na)
        return(.aggbins_sep_inout_na(x.df, bins.df, FUN, ..., empty = empty))

    ag.bins <- aggregate(x.df, by = bins.df, FUN = FUN, ..., drop = drop)
    if (!is.na(empty)) {
        dnames <- names(x.df)
        e.idx <- which(is.na(ag.bins[, dnames]), arr.ind = TRUE)
        ag.bins[, dnames][e.idx] <- empty
    }
    ag.bins
}


#' @importFrom S4Vectors pc
.aggbins_sep_inout_na <- function(x.df, bins.df, FUN, ..., empty) {
    # aggregate returns NA for empty bins; but FUN can also return NAs;
    # if we're not ignoring NA values returned by FUN, AND we're not dropping
    # empty bins, AND we're not setting empty bins to NA, we need to
    # independently identify which bins are empty, and set them to the value
    # given by 'empty' (because is.na() will identify bins which are empty,
    # as well as bins for which FUN returned NA)

    # aggregate in non-empty bins
    ag.bins <- aggregate(x.df, by = bins.df, FUN = FUN, ..., drop = TRUE)

    # get all possible combinations of bins
    bin_comb <- lapply(bins.df, function(x) sort(unique(x)))
    df <- expand.grid(bin_comb)

    ndim <- ncol(bins.df)
    # make list of vectors of dimension values for...
    usedbins <- do.call(pc, as.list(ag.bins[seq_len(ndim)])) # ...used bins
    allbins <- do.call(pc, as.list(df)) # ...all possible bins

    # get indices of non-empty bins (expand.grid and aggregate order the same)
    idx <- which(allbins %in% usedbins)

    # fill values
    datnames <- names(x.df)
    df[, datnames] <- empty
    df[idx, datnames] <- ag.bins[, datnames]
    df
}


#' @rdname binNdimensions
#' @export
densityInNdimBins <- function(dims.df, nbins = 10L, use_bin_numbers = TRUE,
                              ncores = getOption("mc.cores", 2L)) {
    # avoid unnecessary evaluations in the aggregate function
    x <- rep.int(0L, nrow(dims.df))
    bins.df <- binNdimensions(dims.df, nbins, use_bin_numbers, ncores)
    ag.bins <- aggregate(data.frame(value = x), by = bins.df, FUN = length,
                         drop = FALSE)
    ag.bins$value[is.na(ag.bins$value)] <- 0L
    ag.bins
}
