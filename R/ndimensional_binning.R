### =========================================================================
### N-dimensional binning by quantiles
### -------------------------------------------------------------------------
###


#' N-dimensional binning of data by quantiles
#'
#' This function takes in data along 1 or more dimensions, and for each
#' dimension evenly divides the data in evenly-sized quantiles from the minimum
#' to the maximum value. For each input data point, indices are returned giving
#' its bin.
#'
#' @param ... Input data can be a single dataframe or any number of lists or
#' vectors containing different measurements across the same samples.
#' @param quantiles Either a number giving the number of quantiles to use for
#' all dimensions (default = 10), or a vector containing the number of quantiles
#' to use for each dimension of input data given.
#'
#' @return A dataframe containing indices in \code{1:quantiles} for each
#' datapoint in each dimension.
#' @export
#'
#' @examples
binNDimensions <- function(..., quantiles = 10) {

    data_in <- list(...)

    # check input data and convert to dataframe
    input_classes <- vapply(data_in, class, FUN.VALUE = character(1))
    if (all(input_classes) %in% c("list", "numeric")) {
        data_in <- lapply(data_in, unlist)
        data <- as.data.frame(data_in)

    } else if (any(input_classes == "data.frame")) {
        if (length(data_in) > 1) {
            stop(paste("If a dataframe is given as input,",
                       "no additional data objects can be supplied"))
            return(geterrmessage())
        }
        data <- data_in[[1]]

    } else {
        if (length(data_in) == 1) {
            warning(paste("Attempting to coerce input of class",
                          input_classes[[1]], "to dataframe.",
                          "User should verify this operation is valid."))
            data <- as.data.frame(data_in[[1]])
        } else {
            warning(paste("Attempting to coerce a list of input of class(es)",
                          input_classes, "into a dataframe.",
                          "User should verify this operation is valid."))
            data <- as.data.frame(lapply(data_in, as.vector))
        }
    }

    # check input quantiles
    n_dim <- ncol(data)
    if (length(quantiles) > 1 & length(quantiles) != n_dim) {
        stop(paste0("User input ", n_dim, " dimensions of data, ",
                    "but length(quantiles) = ", length(quantiles),
                    ". Quantiles must match number of dimensions, ",
                    "or be a single number."))
        return(geterrmessage())
    }

    # get bin divisions for each dimension, evenly spaced from min to max values
    bin_seqs <- mapply(function(x, y) seq(min(x), max(x), length = y),
                       data, quantiles, SIMPLIFY = FALSE)

    # get bin indices for each datapoint along each dimension
    bin_idx <- mapply(findInterval, data, bin_seqs, SIMPLIFY = FALSE)

    as.data.frame(bin_idx)
}


