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
#'   input data has been replaced by its bin number (an integer).
#'
#'   \code{aggregateByNdimensionalBins} summarizes some input data \code{x} in
#'   each combination of bins, i.e. in each n-dimensional bin. Each row of the
#'   output dataframe is a unique combination of the input bins (i.e. each
#'   n-dimensional bin), and the output columns are identical to those in
#'   \code{dims.df}, with the addition of one or more columns containing the
#'   aggregated data in each n-dimensional bin. If the input \code{x} was a
#'   vector, the column is named "value"; if the input \code{x} was a dataframe,
#'   the column names from \code{x} are maintained.
#'
#'   \code{densityInNdimensionalBins} returns a dataframe just like
#'   \code{aggregateByNdimensionalBins}, except the "value" column contains the
#'   number of observations that fall into each n-dimensional bin.
#'
#' @author Mike DeBerardine
#' @export
#' @importFrom parallel detectCores mcMap mclapply
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
#' bin3d <- binNdimensions(df, nbins = 20, ncores = 2)
#'
#' length(txs_dm6_chr4)
#' nrow(bin3d)
#' bin3d[1:6, ]
#'
#' #--------------------------------------------------#
#' # get number of genes in each bin
#' #--------------------------------------------------#
#'
#' bin_counts <- densityInNdimensionalBins(df, nbins = 20,
#'                                         ncores = 2)
#'
#' bin_counts[1:6, ]
#'
#' #--------------------------------------------------#
#' # get mean cps reads in bins of promoter and genebody reads
#' #--------------------------------------------------#
#'
#' bin2d_cps <- aggregateByNdimensionalBins("counts_cps", df,
#'                                          nbins = 20, ncores = 2)
#'
#' bin2d_cps[1:6, ]
#'
#' subset(bin2d_cps, is.finite(value))[1:6, ]
#'
#' #--------------------------------------------------#
#' # get median cps reads for those bins
#' #--------------------------------------------------#
#'
#' bin2d_cps_med <- aggregateByNdimensionalBins("counts_cps", df, nbins = 20,
#'                                              FUN = median, ncores = 2)
#'
#' bin2d_cps_med[1:6, ]
#'
#' subset(bin2d_cps_med, is.finite(value))[1:6, ]
binNdimensions <- function(dims.df, nbins = 10, ncores = detectCores()) {

    # check input bins
    n_dim <- ncol(dims.df)
    if (length(nbins) > 1 & length(nbins) != n_dim)
        stop(.nicemsg("User input %d dimensions of data, but length(nbins) = %d.
                      length(nbins) must be equal to ncol(dims.df), or be a
                      single number.", n_dim, length(nbins)))

    # get bin divisions for each dimension, evenly spaced from min to max values
    getBreaks <- function(x, y) {
        xfin <- x[is.finite(x)]
        seq(min(xfin), max(xfin), length = y)
    }
    breaks <- mcMap(getBreaks, dims.df, nbins, mc.cores = ncores)

    # get bin indices for each datapoint along each dimension
    bin_idx <- mcMap(findInterval, dims.df, breaks, mc.cores = ncores)
    names(bin_idx) <- paste0("bin.", names(dims.df))
    as.data.frame(bin_idx)
}




#' @rdname binNdimensions
#' @importFrom parallel detectCores
#' @export
aggregateByNdimensionalBins <- function(x, dims.df, nbins = 10, FUN = mean, ...,
                                        ignore.na = TRUE, drop = FALSE,
                                        empty = NA, ncores = detectCores()) {
    if (is.character(x))  {
        colx <- which(names(dims.df) == x)
        x <- dims.df[, colx]
        dims.df <- dims.df[, -colx]
    }
    if (ignore.na) {
        idx <- !is.na(x)
        x <- x[idx]
        dims.df <- dims.df[idx, ]
    }

    # if x is a dataframe, keep the original names
    if (is.vector(x)) {
        x.df <- data.frame(value = x)
    } else {
        x.df <- x
    }

    # get bins for data
    bins.df <- binNdimensions(dims.df, nbins, ncores)

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

    # make list of vectors of dimension values for...
    ndim <- ncol(bins.df)
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
#' @importFrom parallel detectCores
#' @export
densityInNdimensionalBins <- function(dims.df, nbins = 10,
                                      ncores = detectCores()) {
    # avoid unnecessary evaluations in the aggregate function
    x <- rep(0L, nrow(dims.df))
    bins.df <- binNdimensions(dims.df, nbins, ncores)
    ag.bins <- aggregate(data.frame(value = x), by = bins.df, FUN = length,
                         drop = FALSE)
    ag.bins$value[is.na(ag.bins$value)] <- 0
    ag.bins
}
