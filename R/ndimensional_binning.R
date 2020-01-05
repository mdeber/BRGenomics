### ========================================================================= #
### N-dimensional binning by nbins
### ------------------------------------------------------------------------- #
###


#' N-dimensional binning
#'
#' This function takes in data along 1 or more dimensions, and for each
#' dimension the data is divided into evenly-sized bins from the minimum
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
#' @param nbins Either a number giving the number of bins to use for
#'   all dimensions (default = 10), or a vector containing the number of
#'   bins to use for each dimension of input data given.
#'
#' @return A dataframe containing indices in \code{1:nbins} for each
#'   datapoint in each dimension.
#'
#' @author Mike DeBerardine
#'
#' @export
#'
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
#' counts_pr <- getCountsByRegions(PROseq, pr)
#' counts_gb <- getCountsByRegions(PROseq, early_gb)
#' counts_cps <- getCountsByRegions(PROseq, cps)
#'
#' #--------------------------------------------------#
#' # divide genes into 20 bins for each measurement
#' #--------------------------------------------------#
#'
#' count_bins <- binNdimensions(counts_pr, counts_gb, counts_cps, nbins = 20)
#'
#' length(txs_dm6_chr4)
#' nrow(count_bins)
#' count_bins[1:10, ]
binNdimensions <- function(..., nbins = 10) {

    # check input data; convert to dataframe; give good dim. names
    in.names <- .get_unnamed_names() # names of args/objects in '...'
    data <- .get_data(list(...), in.names)

    # check input bins
    n_dim <- ncol(data)
    if (length(nbins) > 1 & length(nbins) != n_dim) {
        stop(.nicemsg("User input %d dimensions of data, but length(nbins) = %d.
                      nbins must match number of dimensions, or be a single
                      number.", n_dim, length(nbins)))
        return(geterrmessage())
    }

    # get bin divisions for each dimension, evenly spaced from min to max values
    seqRange <- function(x, y) seq(min(Filter(is.finite, x)),
                                   max(Filter(is.finite, x)),
                                   length = y)
    bin_seqs <- Map(seqRange, data, nbins)

    # get bin indices for each datapoint along each dimension
    bin_idx <- Map(findInterval, data, bin_seqs)
    names(bin_idx) <- paste0("bin.", names(data))

    as.data.frame(bin_idx)
}


.get_data <- function(in.data, in.names) {
    classes.in <- vapply(in.data, class, FUN.VALUE = character(1))
    if (all(classes.in %in% c("list", "numeric", "integer"))) {
        data <- as.data.frame( lapply(in.data, unlist) )
        names(data) <- in.names
        return(data)

    } else if (any(classes.in == "data.frame")) {
        if (length(in.data) > 1) {
            stop(.nicemsg("If a dataframe is given as input, no additional data
                          objects can be given."))
            return(geterrmessage())
        }
        return(in.data[[1]])

    } else if (length(in.data) == 1) {
        warning(.nicemsg("Coercing input of class %s into adataframe...",
                         classes.in[[1]]), immediate. = TRUE)
        return(as.data.frame(in.data[[1]]))

    } else {
        data <- as.data.frame(lapply(in.data, as.vector))
        warning(.nicemsg("Coerced a list of input of class(es) %s into a
                         dataframe.",
                         Reduce(function(...) paste(..., sep = ","),
                                classes.in)))
        names(data) <- in.names
        return(data)
    }
}


.get_unnamed_names <- function(...) {
    # This function returns a character vector providing names for unnamed args
    #   passed to a parent function. If arguments were named by user (i.e.
    #   x = data_x, ...), it uses those names; but otherwise uses the names of
    #   the objects themselves.
    # This function excludes the named argument "nbins"
    dim_names <- as.list( match.call(call = sys.call(sys.parent(1))) )[-1]
    if (!is.null(names(dim_names))) {
        dim_names <- dim_names[names(dim_names) != "nbins"]
        # if nbins the only named argument, or not all named, don't return
        if (all(nchar(names(dim_names)) > 0)) {
            return(names(dim_names))
        }
    }
    return(as.character(dim_names))
}


