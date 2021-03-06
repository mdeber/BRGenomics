### ========================================================================= #
### Bootstrapping mean signal counts
### ------------------------------------------------------------------------- #
###


#' Bootstrapping Mean Signal by Position for Metaplotting
#'
#' These functions perform bootstrap subsampling of mean readcounts at different
#' positions within regions of interest (\code{metaSubsample}), or, in the more
#' general case of \code{metaSubsampleMatrix}, column means of a matrix are
#' bootstrapped by sampling the rows. Mean signal counts can be calculated at
#' base-pair resolution, or over larger bins.
#'
#' @param dataset.gr A GRanges object in which signal is contained in metadata
#'   (typically in the \code{"score"} field), or a list of such GRanges objects.
#' @param regions.gr A GRanges object containing intervals over which to
#'   metaplot. All ranges must have the same width.
#' @param counts.mat A matrix over which to bootstrap column means by
#'   subsampling its rows. Typically, a matrix of readcounts with rows for genes
#'   and columns for positions within those genes.
#' @param binsize The size of bin (in basepairs, or number of columns for
#'   \code{metaSubsampleMatrix}) to use for counting signal. Especially
#'   important for counting signal over large or sparse regions.
#' @param first.output.xval The relative start position of the first bin, e.g.
#'   if \code{regions.gr} begins at 50 bases upstream of the TSS, set
#'   \code{first.output.xval = -50}. This number only affects the x-values that
#'   are returned, which are provided as a convenience.
#' @param sample.name Defaults to the name of the input dataset. This is
#'   included in the output as a convenience, as it allows row-binding outputs
#'   from different samples. If \code{length(field) > 1} and the default
#'   \code{sample.name} is left, the sample names will be inferred from the
#'   field names.
#' @param n.iter Number of random subsampling iterations to perform. Default is
#'   \code{1000}.
#' @param prop.sample The proportion of the ranges in \code{regions.gr} (e.g.
#'   the proportion of genes) or the proportion of rows in \code{counts.mat} to
#'   sample in each iteration. The default is \code{0.1} (10 percent).
#' @param lower,upper The lower and upper quantiles of subsampled signal means
#'   to return. The defaults, \code{0.125} and \code{0.875} (i.e. the 12.5th and
#'   85.5th percentiles) return a 75 percent confidence interval about the
#'   bootstrapped mean.
#' @param field One or more metadata fields of \code{dataset.gr} to be counted.
#' @param NF An optional normalization factor by which to multiply the counts.
#'   If given, \code{length(NF)} must be equal to \code{length(field)}.
#' @param remove.empty A logical indicating whether regions
#'   (\code{metaSubsample}) or rows (\code{metaSubsampleMatrix}) without signal
#'   should be removed from the analysis. Not recommended if using multiple
#'   fields, as the gene lists will no longer be equivalent.
#' @param blacklist An optional GRanges object containing regions that should be
#'   excluded from signal counting.
#' @param zero_blacklisted When set to \code{FALSE} (the default), signal at
#'   blacklisted sites is ignored from computations. If set to \code{TRUE},
#'   signal at blacklisted sites will be treated as equal to zero. For
#'   bootstrapping, the default behavior of ignoring signal at blacklisted sites
#'   should almost always be kept.
#' @param expand_ranges Logical indicating if ranges in \code{dataset.gr} should
#'   be treated as descriptions of single molecules (\code{FALSE}), or if ranges
#'   should be treated as representing multiple adjacent positions with the same
#'   signal (\code{TRUE}). See \code{\link[BRGenomics:getCountsByRegions]{
#'   getCountsByRegions}}.
#' @param ncores Number of cores to use for computations.
#'
#' @return A dataframe containing x-values, means, lower quantiles, upper
#'   quantiles, and the sample name (as a convenience for row-binding multiple
#'   of these dataframes). If a list of GRanges is given as input, or if
#'   multiple fields are given, a single, combined dataframe is returned
#'   containing data for all fields/datasets.
#'
#' @author Mike DeBerardine
#' @seealso \code{\link[BRGenomics:getCountsByPositions]{getCountsByPositions}}
#' @name bootstrap-signal-by-position
#'
#' @examples
#' data("PROseq") # import included PROseq data
#' data("txs_dm6_chr4") # import included transcripts
#'
#' # for each transcript, use promoter-proximal region from TSS to +100
#' pr <- promoters(txs_dm6_chr4, 0, 100)
#'
#' #--------------------------------------------------#
#' # Bootstrap average signal in each 5 bp bin across all transcripts,
#' # and get confidence bands for middle 30% of bootstrapped means
#' #--------------------------------------------------#
#'
#' set.seed(11)
#' df <- metaSubsample(PROseq, pr, binsize = 5,
#'                     lower = 0.35, upper = 0.65,
#'                     ncores = 1)
#' df[1:10, ]
#'
#' #--------------------------------------------------#
#' # Plot bootstrapped means with confidence intervals
#' #--------------------------------------------------#
#'
#' plot(mean ~ x, df, type = "l", main = "PROseq Signal",
#'      ylab = "Mean + 30% CI", xlab = "Distance from TSS")
#' polygon(c(df$x, rev(df$x)), c(df$lower, rev(df$upper)),
#'         col = adjustcolor("black", 0.1), border = FALSE)
#'
#'
#' #==================================================#
#' # Using a matrix as input
#' #==================================================#
#'
#' # generate a matrix of counts in each region
#' countsmat <- getCountsByPositions(PROseq, pr)
#' dim(countsmat)
#'
#' #--------------------------------------------------#
#' # bootstrap average signal in 10 bp bins across all transcripts
#' #--------------------------------------------------#
#'
#' set.seed(11)
#'  df <- metaSubsampleMatrix(countsmat, binsize = 10,
#'                            sample.name = "PROseq",
#'                            ncores = 1)
#' df[1:10, ]
#'
#' #--------------------------------------------------#
#' # the same, using a normalization factor, and changing the x-values
#' #--------------------------------------------------#
#'
#' set.seed(11)
#' df <- metaSubsampleMatrix(countsmat, binsize = 10,
#'                           first.output.xval = 0, NF = 0.75,
#'                           sample.name = "PROseq", ncores = 1)
#' df[1:10, ]
NULL


#' @rdname bootstrap-signal-by-position
#' @export
#' @importFrom parallel mcMap
#' @importFrom GenomicRanges width
metaSubsample <- function(
    dataset.gr, regions.gr, binsize = 1L, first.output.xval = 1L,
    sample.name = deparse(substitute(dataset.gr)), n.iter = 1000L,
    prop.sample = 0.1, lower = 0.125, upper = 0.875, field = "score",
    NF = NULL, remove.empty = FALSE, blacklist = NULL, zero_blacklisted = FALSE,
    expand_ranges = FALSE, ncores = getOption("mc.cores", 2L)
) {
    if (length(unique(width(regions.gr))) > 1L) # check ranges all same width
        stop(message = "Not all ranges in regions.gr are the same width")

    # Signal in each gene; matrix of dim = (ngenes, nbins), or list of matrices
    signal.bins <- getCountsByPositions(
        dataset.gr, regions.gr, binsize = binsize, field = field,
        blacklist = blacklist, NA_blacklisted = !zero_blacklisted,
        expand_ranges = expand_ranges, ncores = ncores
    )

    if (is.list(signal.bins)) { # (if length(dataset.gr) > 1)
        if (length(sample.name) != length(signal.bins))
            sample.name <- names(signal.bins)
        if (length(sample.name) == 0L)
            sample.name <- as.character(seq_along(signal.bins))
        if (remove.empty)
            warning("remove.empty set with multiple fields/datasets")
        ncores_in <- 1L
        ncores_out <- ncores
    } else {
        ncores_in <- ncores
        ncores_out <- 1L
        signal.bins <- list(signal.bins)
    }
    if (is.null(NF))  NF <- rep.int(1L, length(signal.bins))

    cl <- mcMap(metaSubsampleMatrix, signal.bins, binsize = 1L,
                sample.name = sample.name, n.iter = n.iter,
                prop.sample = prop.sample, lower = lower, upper = upper,
                NF = NF, remove.empty = remove.empty, ncores = ncores_in,
                mc.cores = ncores_out, mc.set.seed = FALSE)

    # fix x-values to match bins, and binsize-normalize the returned values
    if (length(cl) == 1L) {
        .fixbins(cl[[1L]], binsize, first.output.xval)
    } else {
        cl <- lapply(cl, .fixbins, binsize, first.output.xval)
        do.call(rbind, c(cl, make.row.names = FALSE))
    }
}


.fixbins <- function(df, binsize, first.output.xval) {
    nbins <- nrow(df)
    if (binsize == 1L) {
        df$x <- seq(0L, nbins - 1L) + first.output.xval
    } else {
        df$x <- .binxval(nbins, binsize, first.output.xval)
        y_vals <- c("mean", "lower", "upper")
        df[, y_vals] <- df[, y_vals] / binsize
    }
    rownames(df) <- NULL
    return(df)
}


#' @importFrom stats median
.binxval <- function(nbins, binsize, first.output.xval) {
    firstbin <- seq(first.output.xval, by = 1L, length.out = binsize)
    binstart <- median(firstbin) # center of first bin
    seq(binstart, by = binsize, length.out = nbins)
}



#' @rdname bootstrap-signal-by-position
#' @export
#' @importFrom parallel mclapply
#' @importFrom stats quantile
metaSubsampleMatrix <- function(counts.mat, binsize = 1L, first.output.xval = 1L,
                                sample.name = NULL, n.iter = 1000L,
                                prop.sample = 0.1, lower = 0.125, upper = 0.875,
                                NF = 1L, remove.empty = FALSE,
                                ncores = getOption("mc.cores", 2L)) {
    # Check that enough iterations are given for meaningful quantiles
    if (n.iter != 1L) .check_iter(n.iter, lower, upper)
    if (remove.empty)
        counts.mat <- counts.mat[rowSums(counts.mat) > 0L, , drop = FALSE]
    if (is.null(sample.name)) sample.name <- deparse(substitute(counts.mat))

    nbins <- floor(ncol(counts.mat) / binsize)
    ngenes <- nrow(counts.mat)
    nsample <- ceiling(prop.sample*ngenes)

    if (binsize > 1L) {
        counts.mat <- apply(counts.mat, 1L, .binVector, binsize = binsize)
        counts.mat <- if (is.array(counts.mat)) t(counts.mat) else
            matrix(counts.mat) # (for single column matrices)
    }
    # 1. Randomly subsample rows of the counts.mat
    #    -> List of length = n.iter, containing vectors of length = nsample
    # 2. For each iteration, get average signal in each bin
    #    -> List of length = n.iter containing vectors of length = nbins
    # 3. Collapse list into matrix
    #    -> Matrix of dim = (nbins, n.iter)
    idx.list <- replicate(n.iter, sample(ngenes, nsample), simplify = FALSE)
    binavg <- mclapply(idx.list,
                       function(idx) colMeans(counts.mat[idx, , drop = FALSE],
                                              na.rm = TRUE),
                       mc.cores = ncores)
    binavg.mat <- matrix(unlist(binavg), ncol = n.iter)

    # calculate final outputs
    if (n.iter == 1L) {
        message(.nicemsg("With n.iter = 1, output means and quantiles are not
                         bootstrapped"))
        idx <- unlist(idx.list)
        mean <- NF*binavg.mat
        lower <- NF*apply(counts.mat[idx, ], 2L, quantile, lower, na.rm = TRUE)
        upper <- NF*apply(counts.mat[idx, ], 2L, quantile, upper, na.rm = TRUE)
    } else {
        mean <- NF*apply(binavg.mat, 1L, median, na.rm = TRUE)
        lower <- NF*apply(binavg.mat, 1L, quantile, lower,
                          na.rm = TRUE, names = FALSE)
        upper <- NF*apply(binavg.mat, 1L, quantile, upper,
                          na.rm = TRUE, names = FALSE)
    }

    # Also return x-values, sample names for plotting; x-vals centered in bins
    if (binsize == 1L) x <- first.output.xval + seq(0L, nbins - 1L)
    if (binsize > 1L) x <- .binxval(nbins, binsize, first.output.xval)
    return(data.frame(x, mean, lower, upper, sample.name))
}


.check_iter <- function(n.iter, lower, upper) {
    checks <- c(round(n.iter*lower) == 0,
                round(n.iter*lower) == 0.5*n.iter,
                round(n.iter*upper) == 0.5*n.iter,
                round(n.iter*upper) == n.iter)
    if (any(checks))
        stop("Insufficient iterations to obtain distinct values of quantiles")
}
