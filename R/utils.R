# ========================================================================= #
# Interal helper functions
# ------------------------------------------------------------------------- #


.nicemsg <- function(...) {
    # format output msgs; removes newlines and leading spaces from strings
    strwrap(sprintf(...), prefix = " ", initial = "")
}

.binVector <- function(x, binsize = NULL, nbins = NULL, FUN = sum) {
    # superimpose evenly spaced bins over a vector, and perform FUN on them
    #   (output length = nbins)
    if (missing(nbins))
        nbins <- floor(length(x) / binsize)
    if (missing(binsize))
        binsize <- floor(length(x) / nbins)
    mat <- matrix(x[seq_len(nbins*binsize)], nrow = binsize)
    apply(mat, 2, FUN)
}


.check_xor_args <- function(arg1, arg2) {
    # check two mutually exclusive arguments, which default to NULL
    name1 <- deparse(substitute(arg1))
    name2 <- deparse(substitute(arg2))
    if (!xor(is.null(arg1), is.null(arg2)))
        stop(sprintf("provide either %s or %s, but not both", name1, name2))
}


.dfList2df <- function(df.list, col_name = "sample", prepend = TRUE) {
    # convert a list of dataframes into a single dataframe, with a column
    # added for sample names
    df.list <- lapply(seq_along(df.list), function(i) {
        df <- df.list[[i]]
        df[[col_name]] <- names(df.list)[i]
        if (prepend)
            df <- df[, c(ncol(df), seq(1, ncol(df) - 1))]
        df
    })
    do.call(rbind, df.list)
}
