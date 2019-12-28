### =========================================================================
### N-dimensional binning by quantiles
### -------------------------------------------------------------------------
###


#' N-dimensional binning of data by quantiles
#'
#' This function takes in data along 1 or more dimensions, and for each
#' dimension the data is divided into evenly-sized quantiles from the minimum
#' value to the maximum value, and bin numbers are returned. For instance, if
#' each index of the input data were a gene, the input dimensions would be
#' various quantitative measures of that gene, e.g. expression level, number of
#' exons, length, etc. If plotted in cartesian coordinates, each gene would be a
#' single datapoint, and each measurement would be a separate dimension. The bin
#' numbers for each datapoint in each dimension are returned in a dataframe,
#' with a column for each dimension and a row for each index.
#'
#' @param ... A single dataframe, or any number of lists or vectors containing
#'   different measurements across the same datapoints. If a dataframe is given,
#'   columns should correspond to measurements (dimensions). If lists or
#'   vectors are given, they must all have the same lengths. Other input classes
#'   will be coerced into a single dataframe.
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
    if (all(input_classes %in% c("list", "numeric"))) {
        data <- as.data.frame( lapply(data_in, unlist) )
        dim_names <- .get_unnamed_names()

    } else if (any(input_classes == "data.frame")) {
        if (length(data_in) > 1) {
            stop(.nicemsg("If a dataframe is given as input, no additional
                          data objects can be given."))
            return(geterrmessage())
        }
        data <- data_in[[1]]
        dim_names <- names(data)

    } else {
        if (length(data_in) == 1) {
            data <- as.data.frame(data_in[[1]])
            warning(.nicemsg("Coerced input of class %s into a
                             dataframe.", input_classes[[1]]))
            dim_names <- names(data)

        } else {
            data <- as.data.frame(lapply(data_in, as.vector))
            warning(.nicemsg("Coerced a list of input of class(es)
                             %s into a dataframe.",
                             Reduce(function(...) paste(..., sep = ","),
                                    input_classes)))
            dim_names <- .get_unnamed_names()
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
    seqRange <- function(x, y) seq(min(Filter(is.finite, x)),
                                   max(Filter(is.finite, x)),
                                   length = y)
    bin_seqs <- mapply(seqRange, data, quantiles, SIMPLIFY = FALSE)

    # get bin indices for each datapoint along each dimension
    bin_idx <- mapply(findInterval, data, bin_seqs, SIMPLIFY = FALSE)
    names(bin_idx) <- paste("bin", dim_names)

    as.data.frame(bin_idx)
}


.get_unnamed_names <- function(...) {
    # This function returns a character vector providing names for unnamed args
    #   passed to a parent function. If arguments were named by user (i.e.
    #   x = data_x, ...), it uses those names; but otherwise uses the names of
    #   the objects themselves.
    # This function excludes the named argument "quantiles"
    dim_names <- as.list( match.call(call = sys.call(sys.parent(1))) )[-1]
    if (!is.null(names(dim_names))) {
        dim_names <- dim_names[names(dim_names) != "quantiles"]
        # if quantiles the only named argument, or not all named, don't return
        if (all(nchar(names(dim_names)) > 0)) {
            return(names(dim_names))
        }
    }
    return(as.character(dim_names))
}


