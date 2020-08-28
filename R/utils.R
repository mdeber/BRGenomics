# ========================================================================= #
# Interal helper functions
# ------------------------------------------------------------------------- #

mcMap <- function(f, ...) {
    # Current windows fxn passes directly to Map, putting mc.cores in the dots;
    # -> Contacted R-for-Windows devs & this will be fixed in next release
    f <- match.fun(f)
    mcmapply(f, ..., SIMPLIFY = FALSE, mc.silent = TRUE)
}

.nicemsg <- function(...) {
    # format output msgs; removes newlines and leading spaces from strings
    strwrap(sprintf(...), prefix = " ", initial = "")
}

.close_int <- function(x) {
    # for a numeric vector x, can it safely be coerced to integer?
    all( abs(round(x) - x) < 10e-12 )
}

.binVector <- function(x, binsize = NULL, nbins = NULL, FUN = sum) {
    # superimpose evenly spaced bins over a vector, and perform FUN on them
    #   (output length = nbins)
    # (all conditions for safety)
    if (is.null(nbins)) {
        binsize <- as.integer(binsize)
        nbins <- as.integer(length(x) / binsize)
    } else if (is.null(binsize)) {
        nbins <- as.integer(nbins)
        binsize <- as.integer(length(x) / nbins)
    } else {
        nbins <- as.integer(nbins)
        binsize <- as.integer(binsize)
    }
    mat <- matrix(x[seq_len(nbins*binsize)], nrow = binsize)
    apply(mat, 2L, FUN)
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
            df <- df[, c(ncol(df), seq(1L, ncol(df) - 1L))]
        df
    })
    do.call(rbind, df.list)
}
