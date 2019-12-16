### =========================================================================
### Bootstrapping mean signal counts
### -------------------------------------------------------------------------
###
### The essential boostrapping function performs subsampling on matrix rows

#' Iterative Subsampling for Metaplotting (On Count Matrices)
#'
#' In the most general sense, this function performs iterations of randomly
#' subsampling rows of a matrix, and
#' returns a summary of mean values calculated for each column. The typical
#' application is for generating
#' metaplots, with the typical input being a matrix in which each row is a
#' gene or other region of interest,
#' each column is a position within that gene (either a specific basepair or
#' a bin), and element values are
#' signal (e.g. read counts) within those positions.
#'
#' @param counts.mat A matrix of signal counts in which rows are regions of
#' interest and columns are
#' sites/bins in each region.
#' @param binsize The size of bin (number of columns, e.g. basepairs) to use
#' for metaplotting. Especially
#' important for metaplots over large/sparse regions.
#' @param first.bin.start The relative start position of the first bin, e.g.
#' if regions.gr begins at 50 bases
#' upstream of the TSS, set \code{first.bin.start = -50}. This number only
#' affects the x-values that are
#' returned, which are provided as a convenience.
#' @param sample_name Defaults to the name of \code{dataset.gr}.
#' @param n.iter Number of random subsampling iterations to perform.
#' Default is 1000.
#' @param prop.subsample The proportion of rows to subsample in each iteration.
#' The default is 0.1.
#' @param lower The lower quantile of subsampled signal means to return. The
#' default is
#' 0.125 (12.5th percentile).
#' @param upper The upper quantile of subsampled signal means to return. The
#' default is
#' 0.875 (85.5th percentile).
#' @param NF Optional normalization factor by which to multiply the counts.
#' @param ncores Number of cores to use for parallel computation. As of writing,
#' parallel processing doesn't
#' show any benefit for short computation times (e.g. <1 minute for our typical
#' experience on a laptop).
#'
#' @return Dataframe containing x-values, means, lower quantiles, upper
#' quantiles, and the sample name (as a
#' convenience for row-binding multiple of these dataframes).
#' @export
md.meta.subsample.matrix <- function(counts.mat,
                                     binsize = 1,
                                     first.bin.start = 1,
                                     sample_name = deparse(substitute(dataset.gr)),
                                     n.iter = 1000,
                                     prop.subsample = 0.1,
                                     lower = 0.125,
                                     upper = 0.875,
                                     NF = 1,
                                     ncores = 1) {

  # Check that enough iterations are given such that quantiles can be calculated
  if ((round(n.iter*lower) == 0) |
      (round(n.iter*lower) == 0.5*n.iter) |
      (round(n.iter*upper) == 0.5*n.iter) |
      (round(n.iter*upper) == n.iter)) {
    stop(message = "Insufficient iterations to obtain distinct values of quantiles")
    return(geterrmessage())
  }

  n.bins <- floor(ncol(counts.mat) / binsize)
  n.genes <- nrow(counts.mat)
  n.sample <- round(prop.subsample*n.genes)

  # (This section can be re-written using the findInterval function)
  if (binsize > 1) {
    # column indices for each bin's start and end column:
    bin.starts <- seq(1, binsize*n.bins, binsize)
    bin.ends <- seq(binsize, binsize*n.bins, binsize)
    # for each row (region) in counts.mat, reassign as a vector of length n.bins
    # each column the sum of signal within each bin
    counts.mat <- apply(counts.mat, 1, function(row.i) {
      vapply(seq_len(n.bins),
             function(bin.i) sum(row.i[ bin.starts[bin.i]:bin.ends[bin.i] ]),
             FUN.VALUE = numeric(1))
    })
    # (transpose, as sapply/vapply will cbind the input rows)
    counts.mat <- t(counts.mat) # final dim = (genes, bins)
  }

  if (ncores == 1) {
    # Randomly subsample rows of counts.mat
    # -> Matrix of dim = (n.sample, n.iter)
    idx.gene.iters <- replicate(n.iter,
                                sample(n.genes, size = n.sample),
                                simplify = "array")

    # Find average signal in each bin (column of counts.mat) across all subsampling iterations
    # -> Matrix of dim = (n.bins, n.iter)
    avgsignal.gene.iters <- apply(idx.gene.iters,
                                  2,
                                  function(idx.genes) colMeans(counts.mat[idx.genes, ]))
  } else {
    # Randomly subsample rows of the counts.mat
    # -> List of length = n.iter, containing vectors of length = n.sample
    idx.gene.iters <- replicate(n.iter,
                                sample(n.genes, size = n.sample),
                                simplify = FALSE)

    # Find average signal in each bin (column of counts.mat) across all subsampling iterations
    # -> List of length = n.iter containing vectors of length = n.bins
    avgsignal.gene.iters <- mclapply(idx.gene.iters,
                                     function(idx.genes) colMeans(counts.mat[idx.genes, ]),
                                     mc.cores = ncores)

    # Collapse list into matrix
    # -> Matrix of dim = (n.bins, n.iter)
    avgsignal.gene.iters <- matrix(unlist(avgsignal.gene.iters), ncol = n.iter)
  }

  # Calculate final outputs: quantiles of means across iterations
  mean <- NF * apply(avgsignal.gene.iters, 1,
                     function(bin.i) quantile(bin.i, 0.5, names = F))
  lower <- NF * apply(avgsignal.gene.iters, 1,
                      function(bin.i) quantile(bin.i, lower, names = F))
  upper <- NF * apply(avgsignal.gene.iters, 1,
                      function(bin.i) quantile(bin.i, upper, names = F))

  # Also return x-values and sample names for plotting; x-values centered in bins
  if (binsize == 1) {
    x <- seq(0, n.bins - 1) + first.bin.start
  } else {
    x <- seq(0.5*binsize, n.bins*binsize - 0.5*binsize, binsize) + first.bin.start
  }

  return(data.frame(x, mean, lower, upper, sample_name))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Bootstrapping signal means from GRanges
###


#' Iterative Subsampling for Metaplotting
#'
#' @param dataset.gr A GRanges object in which signal is contained in metadata
#' (typically in the "score" field).
#' @param regions.gr A GRanges object containing intervals over which to
#' metaplot. All ranges must have
#' the same width.
#' @param binsize The size of bin (number of columns, e.g. basepairs) to use for
#' metaplotting. Especially
#' important for metaplots over large/sparse regions.
#' @param first.bin.start The relative start position of the first bin, e.g. if
#' regions.gr begins at 50 bases
#' upstream of the TSS, set \code{first.bin.start = -50}. This number only
#' affects the x-values that are
#' returned, which are provided as a convenience.
#' @param sample_name Defaults to the name of \code{dataset.gr}.
#' @param n.iter Number of random subsampling iterations to perform. Default
#' is 1000.
#' @param prop.subsample The proportion of the ranges in regions.gr (e.g. the
#' proportion of genes) to
#' subsample in each iteration. The default is 0.1.
#' @param lower The lower quantile of subsampled signal means to return. The
#' default is
#' 0.125 (12.5th percentile).
#' @param upper The upper quantile of subsampled signal means to return. The
#' default is
#' 0.875 (85.5th percentile).
#' @param NF Optional normalization factor by which to multiply the counts.
#' @param field The metadata field of \code{dataset.gr} to be counted.
#' @param remove.empty A logical indicating whether regions without signal
#' should be removed from the analysis.
#' @param ncores Number of cores to use for parallel computation. As of writing,
#' parallel processing doesn't
#' show any benefit for short computation times (e.g. <1 minute for our typical
#' experience on a laptop).
#'
#' @return Dataframe containing x-values, means, lower quantiles, upper quantiles,
#'  and the sample name (as a
#' convenience for row-binding multiple of these dataframes).
#' @export
md.meta.subsample <- function(dataset.gr,
                              regions.gr,
                              binsize = 1,
                              first.bin.start = 1,
                              sample_name = deparse(substitute(dataset.gr)),
                              n.iter = 1000,
                              prop.subsample = 0.1,
                              lower = 0.125,
                              upper = 0.875,
                              NF = 1,
                              field = "score",
                              remove.empty = FALSE,
                              ncores = 1) {

  if (length(unique(width(regions.gr))) > 1) {
    stop(message = "Not all ranges in regions.gr are the same width")
    return(geterrmessage())
  }

  # Get signal in each bin of each gene
  # -> Matrix of dim = (n.genes, n.bins)
  signal.bins <- md.get.counts_matrix(dataset.gr = dataset.gr,
                                      regions.gr = regions.gr,
                                      binsize = binsize,
                                      field = field,
                                      remove.empty = remove.empty)

  metamat <- md.meta.subsample.matrix(counts.mat = signal.bins,
                                      binsize = 1,
                                      first.bin.start = first.bin.start,
                                      sample_name = sample_name,
                                      n.iter = n.iter,
                                      prop.subsample = prop.subsample,
                                      lower = lower,
                                      upper = upper,
                                      NF = NF,
                                      ncores = ncores)

  # fix x-values to match bins, and binsize-normalize the returned values
  if (binsize != 1) {
    n.bins <- nrow(metamat)
    metamat$x <- seq(0.5*binsize, n.bins*binsize - 0.5*binsize, binsize) + first.bin.start
    metamat[, c("mean", "lower", "upper")] <- metamat[, c("mean", "lower", "upper")] / binsize
  }

  return(metamat)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Bootstrapping signal means over variable-length ranges
###


#' Iterative Subsampling for Metaplotting (With Scaled Regions)
#'
#' This function can perform iterative metaplot subsampling in several configurations. All arguments for
#' regions-of-interest for metaplot subsampling are optional, as long as at least one region is supplied.
#' Unnamed arguments are regions-of-interest of variable "widths" (i.e. "lengths" in basepairs) over which to
#' perform length-scaled metaplot subsampling. Length-scaled metaplot subsampling involves dividing each
#' range (e.g. each region-of-interest) into \code{nbins_scaled} number of equally-sized bins, and obtaining
#' signal counts in each bin, divided by the size of the bin for that particular region-of-interest.
#' Non-length-scaled subsampling (as would be done using \code{md.meta.subsample}) is performed on the named
#' arguments \code{linear_regions_start.gr} and \code{linear_regions_end.gr}. The output is constructed in
#' the order \code{linear_regions_start.gr}, unnamed scaled regions (in order given), and then
#' \code{linear_regions_end.gr}, with x-values corresponding to the bin number.
#'
#' The user must be able to determine the correct meaning of the bin numbers in the final output, and for
#' that reason arguments for binning are always explicitly the number of bins, and not the size of the bins
#' (as would be possible for linear (un-scaled) regions). For example, if the user provides
#' \code{linear_regions_start.gr}, one unnamed GRanges for length-scaled subsampling, and
#' \code{linear_regions_end.gr}, the output x-values \code{1:nbins_linear_start} will correspond to
#' equally-sized bins in \code{linear_regions_start.gr}; the subsequent \code{nbins_scaled} x-values will
#' correspond to variably-sized bins in the unnamed GRanges object given, and the final
#' \code{nbins_linear_end} x-values will correspond to equally-sized bins in \code{linear_regions_end.gr}.
#'
#' @param dataset.gr A GRanges object in which signal is contained in metadata (typically in the "score" field).
#' @param ... 0 or more GRanges objects containing regions of interest over which to do length-scaled signal
#' counting and metaplot subsampling. The output x-positions will be determined by the order in which these
#' regions are supplied, the number of bins used for counting signal within variable length regions
#' (\code{nbins_scaled}), and whether or not a \code{linear_regions_start.gr} object is given.
#' @param linear_regions_start.gr Optional GRanges object containing regions of interest over which to do
#' linear (un-scaled) signal counting and metaplot subsampling. Because no length-scaling is performed,
#' all ranges must have the same width. These regions will be put before any supplied regions for
#' length-scaled metaplot subsampling, i.e. the first \code{nbins_linear_start} x-values will be from
#' subsampling \code{linear_regions_start.gr}.
#' @param linear_regions_end.gr Optional GRanges object containing regions of interest over which to do
#' linear (un-scaled) signal counting and metaplot subsampling. Because no length-scaling is performed,
#' all ranges must have the same width. These regions will be placed after any supplied regions for
#' length-scaled metaplot subsampling, i.e. the last \code{nbins_linear_end} x-values will be from
#' subsampling \code{linear_regions_end.gr}.
#' @param nbins_scaled The number of bins to use for length scaling signal counts.
#' @param nbins_linear_start The number of bins to use for counting signal within \code{linear_regions_start.gr}.
#' Defaults to the width of the regions, i.e. a binsize of 1 (no binning).
#' @param nbins_linear_end The number of bins to use for counting signal within \code{linear_regions_end.gr}.
#' Defaults to the width of the regions, i.e. a binsize of 1 (no binning).
#' @param sample_name Defaults to the name of \code{dataset.gr}.
#' @param n.iter Number of random subsampling iterations to perform. Default is 1000.
#' @param prop.subsample The proportion of the genelist (regions.gr) to subsample in each iteration.
#' The default is 0.1.
#' @param lower The lower quantile of subsampled signal means to return. The default is
#' 0.125 (12.5th percentile).
#' @param upper The upper quantile of subsampled signal means to return. The default is
#' 0.875 (85.5th percentile).
#' @param NF Optional normalization factor by which to multiply the counts.
#' @param field The metadata field of \code{dataset.gr} to be counted.
#' @param remove.empty A logical indicating whether regions without signal should be removed from the analysis.
#' @param ncores Number of cores to use for parallel computation. As of writing, parallel processing doesn't
#' show any benefit for short computation times (e.g. <1 minute for our typical experience on a laptop).
#'
#' @return Dataframe containing x-values, means, lower quantiles, upper quantiles, and the sample name (as a
#' convenience for row-binding multiple output dataframes). X-values correspond to bins based on the input
#' regions given and the specified binsizes to use.
#' @export
#'
#' @examples md.meta.scaled_subsample(my_proseq_data, genes.early_genebodies, genes.late_genebodies,
#' linear_regions_start.gr = genes.promoter_proximal, linear_regions_end.gr = genes.cps_proximal,
#' nbins_scaled = 500, nbins_linear_end = 500)
md.meta.scaled_subsample <- function(dataset.gr,
                                     ...,
                                     linear_regions_start.gr = NULL,
                                     linear_regions_end.gr = NULL,
                                     nbins_scaled = 100,
                                     nbins_linear_start = unique(width(linear_regions_start.gr)),
                                     nbins_linear_end = unique(width(linear_regions_end.gr)),
                                     sample_name = deparse(substitute(dataset.gr)),
                                     n.iter = 1000,
                                     prop.subsample = 0.1,
                                     lower = 0.125,
                                     upper = 0.875,
                                     NF = 1,
                                     field = "score",
                                     remove.empty = FALSE,
                                     ncores = 1) {

  scaled_regions.list <- list(...)

  # -------------------------------------------------- #
  # Check input
  # -------------------------------------------------- #

  # Check linear ranges (if given)

  if (!is.null(linear_regions_start.gr) &&
      length(unique(width(linear_regions_start.gr))) > 1) {
    stop(message = "Not all ranges in linear_regions_start.gr are the same width")
    return(geterrmessage())
  }

  if (!is.null(linear_regions_end.gr) &&
      length(unique(width(linear_regions_end.gr))) > 1) {
    stop(message = "Not all ranges in linear_regions_end.gr are the same width")
    return(geterrmessage())
  }

  # -------------------------------------------------- #

  # Check that all ranges given have equal lengths (number of genes)

  length_check <- c()

  if (!is.null(linear_regions_start.gr)) {
    length_check <- c(length_check, length(linear_regions_start.gr))
  }

  if (!is.null(linear_regions_end.gr)) {
    length_check <- c(length_check, length(linear_regions_end.gr))
  }

  if (length(scaled_regions.list) > 0) {
    for (k in seq_along(scaled_regions.list)) {
      length_check <- c(length_check, length(scaled_regions.list[[k]]))
    }
  }

  if (length(length_check) == 0) {
    stop(message = "No regions given as input")
    return(geterrmessage())
  } else if (length(unique(length_check)) > 1) {
    stop(message = "Not all regions have the same number of ranges")
    return(geterrmessage())
  }

  # -------------------------------------------------- #

  # Check quantiles

  if ((round(n.iter*lower) == 0) |
      (round(n.iter*lower) == 0.5*n.iter) |
      (round(n.iter*upper) == 0.5*n.iter) |
      (round(n.iter*upper) == n.iter)) {
    stop(message = "Insufficient iterations to obtain distinct values of quantiles")
    return(geterrmessage())
  }


  # -------------------------------------------------- #
  # Functions for generating count matrices
  # -------------------------------------------------- #

  get_linear_counts <- function(regions.gr, nbins) {
    binsize <- floor( unique(width(regions.gr)) / nbins )
    counts.mat <- md.get.counts_matrix(dataset.gr = dataset.gr,
                                       regions.gr = regions.gr,
                                       binsize = binsize,
                                       field = field,
                                       remove.empty = FALSE)
    counts.mat <- counts.mat / binsize
    return(counts.mat)
  }

  get_scaled_counts <- function(regions.gr, nbins) {

    # make tiles for each gene
    tiles <- as(tile(regions.gr, n = nbins), "GRangesList")

    # reverse minus strand tiles, and shift them toward their strand-specific starts (to match plus strand tiles)
    idx.minus.genes <- which( as.character(strand(regions.gr)) == "-" )
    tiles.m <- tiles[idx.minus.genes]
    tiles.m <- lapply(tiles.m, rev)
    widths.m <- width(regions.gr[idx.minus.genes])
    shifts.m <- widths.m %% nbins
    tiles.m <- lapply(seq_along(tiles.m), function(i) shift(tiles.m[[i]], shifts.m[i]))
    # check, but should replace with:
    # tiles.m <- mapply(shift, tiles.m, shifts.m, SIMPLIFY = FALSE)
    tiles[idx.minus.genes] <- as(tiles.m, "GRangesList")

    # make GRanges; every nbin-th range is the same bin in the next gene
    tiles <- unlist(tiles)

    hits.tiles <- findOverlaps(tiles, dataset.gr)
    signal.tiles <- aggregate(mcols(dataset.gr)[[field]][hits.tiles@to],
                              by = list(hits.tiles@from),
                              FUN = sum)
    names(signal.tiles) <- c("tile.idx", "signal")

    # Make signal (binsize-normalized) matrix of dim = (nbins, n.genes), then transpose
    signal.tiles.all <- rep(0, length(tiles)) # initialize
    signal.tiles.all[signal.tiles$tile.idx] <- signal.tiles$signal / width(tiles[signal.tiles$tile.idx])
    signal.tiles.mat <- matrix(signal.tiles.all,
                               nrow = nbins,
                               ncol = n.genes)

    # Final counts matrix: dim = (n.genes, nbins)
    return(t(signal.tiles.mat))
  }


  # -------------------------------------------------- #
  # Create and combine count matrices
  # -------------------------------------------------- #

  # Creates matrix of dim = (n.genes, total number of bins across all given regions)

  n.genes <- unique(length_check)
  n.sample <- round(prop.subsample*n.genes)
  counts.mat <- c()

  # For scaled counts, combine in order given
  if (length(scaled_regions.list) > 0) {
    for (k in seq_along(scaled_regions.list)) {
      counts.mat <- cbind(counts.mat, get_scaled_counts(scaled_regions.list[[k]],
                                                        nbins_scaled))
    }
  }

  # Prepend linear_regions_start counts
  if (!is.null(linear_regions_start.gr)) {
    counts.start <- get_linear_counts(regions.gr = linear_regions_start.gr,
                                      nbins = nbins_linear_start)
    counts.mat <- cbind(counts.start, counts.mat)
  }

  # Append linear_regions_end counts
  if (!is.null(linear_regions_end.gr)) {
    counts.end <- get_linear_counts(regions.gr = linear_regions_end.gr,
                                    nbins = nbins_linear_end)
    counts.mat <- cbind(counts.mat, counts.end)
  }

  # -------------------------------------------------- #

  if (remove.empty) {
    counts.mat <- counts.mat[rowSums(counts.mat) > 0, ]
  }

  # -------------------------------------------------- #
  # Iterative subsampling of final count matrix
  # -------------------------------------------------- #

  md.meta.subsample.matrix(counts.mat = counts.mat,
                           binsize = 1,
                           sample_name = sample_name,
                           n.iter = n.iter,
                           prop.subsample = prop.subsample,
                           lower = lower,
                           upper = upper,
                           NF = NF,
                           ncores = ncores)
}


