### =========================================================================
### N-dimensional binning by quantiles
### -------------------------------------------------------------------------
###


#' N-dimensional binning of data by quantiles
#'
#' This function takes in data along 1 or more dimensions, and for each
#' dimension evenly divides the data in evenly-sized quantiles from the minimum
#' value to the maximum value. For each input data point, indices are returned
#' giving the bin in each dimension.
#'
#' @param ... A single dataframe, or any number of lists or
#'   vectors containing different measurements across the same samples. If lists
#'   or vectors are given, they must all have the same lengths.
#' @param quantiles Either a number giving the number of quantiles to use for
#'   all dimensions (default = 10), or a vector containing the number of
#'   quantiles to use for each dimension of input data given.
#'
#' @return A dataframe containing indices in \code{1:quantiles} for each
#'   datapoint in each dimension.
#'
#' @author Mike DeBerardine
#'
#' @export
#'
#' @examples
binNdimensions <- function(..., quantiles = 10) {

    data_in <- list(...)

    # check input data and convert to dataframe
    input_classes <- vapply(data_in, class, FUN.VALUE = character(1))
    if (all(input_classes) %in% c("list", "numeric")) {
        data <- as.data.frame( lapply(data_in, unlist) )

    } else if (any(input_classes == "data.frame")) {
        if (length(data_in) > 1) {
            stop(.nicemsg("If a dataframe is given as input, no additional
                          data objects can be given."))
            return(geterrmessage())
        }
        data <- data_in[[1]]

    } else {
        if (length(data_in) == 1) {
            warning(.nicemsg("Attempting to coerce input of class %s into a
                             dataframe. User should verify this operation is
                             valid.", input_classes[[1]]))
            data <- as.data.frame(data_in[[1]])
        } else {
            warning(.nicemsg("Attempting to coerce a list of input of class(es)
                             %s into a dataframe. User should verify this
                             operation is valid.",
                             Reduce(function(...) paste(..., sep = ","),
                                    input_classes)))
            data <- as.data.frame(lapply(data_in, as.vector))
        }
    }

    # check input quantiles
    n_dim <- ncol(data)
    if (length(quantiles) > 1 & length(quantiles) != n_dim) {
        stop(.nicemsg("User input %d dimensions of data, but length(quantiles)
                      = %d. Quantiles must match number of dimensions, or be a
                      single number.",
                      n_dim, length(quantiles)))
        return(geterrmessage())
    }

    # get bin divisions for each dimension, evenly spaced from min to max values
    seqRange <- function(x, y) seq(min(x), max(x), length = y)
    bin_seqs <- mapply(seqRange, data, quantiles, SIMPLIFY = FALSE)

    # get bin indices for each datapoint along each dimension
    bin_idx <- mapply(findInterval, data, bin_seqs, SIMPLIFY = FALSE)

    as.data.frame(bin_idx)
}


