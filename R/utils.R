# ========================================================================= #
# Package documentation
# ------------------------------------------------------------------------- #

#' BRGenomics
#'
#' @description BRGenomics provides useful functions for analyzing genomics data
#'   at base-pair resolution, and for doing so in a way that maximizes
#'   compatibility with the wide array of packages available through
#'   Bioconductor.
#'
#'   For interactive documentation with code examples, see the online
#'   documentation: \url{https://mdeber.github.io/}
#'
#' @import GenomicRanges S4Vectors
#' @aliases BRGenomics
#' @name BRGenomics-package
NULL


# ========================================================================= #
# Interal helper functions
# ------------------------------------------------------------------------- #

## for output messages, removes newlines and leading spaces from strings

.nicemsg <- function(...) strwrap(sprintf(...), prefix = " ", initial = "")


## superimposes evenly spaced bins over a vector, and performs FUN on them
## (output length = nbins)

.binVector <- function(x, binsize = NULL, nbins = NULL, FUN = sum) {
    if (missing(nbins))  nbins <- floor(length(x) / binsize)
    if (missing(binsize))  binsize <- floor(length(x) / nbins)

    # bin_breaks <- seq(1, binsize * nbins, binsize)
    # bin_idx <- findInterval(seq_len(binsize * nbins), bin_breaks)
    # aggregate(x[seq_len(binsize*nbins)], by = list(bin_idx), FUN = FUN)[, 2]

    mat <- matrix(x[seq_len(nbins*binsize)], nrow = binsize)
    apply(mat, 2, FUN)
}

