### ========================================================================= #
### DESeq2 helper functions
### ------------------------------------------------------------------------- #
###


#' Get DESeqDataSet objects for downstream analysis
#'
#' This is a convenience function for generating \code{DESeqDataSet} objects,
#' but this function also adds support for counting reads across non-contiguous
#' regions.
#'
#' @param dataset.list An object containing GRanges datasets that can be passed
#'   to \code{\link[BRGenomics:getCountsByRegions]{getCountsByRegions}},
#'   typically a list of GRanges objects, or a
#'   \code{\link[BRGenomics:mergeGRangesData]{ multiplexed GRanges}} object (see
#'   details below).
#' @param regions.gr A GRanges object containing regions of interest.
#' @param sample_names Names for each dataset in \code{dataset.list} are
#'   required. By default (\code{sample_names = NULL}), if \code{dataset.list}
#'   is a list, the names of the list elements are used; for a multiplexed
#'   GRanges object, the field names are used. The names must each contain the
#'   string "_rep#", where "#" is a single character (usually a number)
#'   indicating the replicate. Sample names across different replicates must be
#'   otherwise identical.
#' @param gene_names An optional character vector giving gene names, or any
#'   other identifier over which reads should be counted. Gene names are
#'   required if counting is to be performed over non-contiguous ranges, i.e. if
#'   any genes have multiple ranges. If supplied, gene names are added to the
#'   resulting \code{DESeqDataSet} object.
#' @param sizeFactors DESeq2 \code{sizeFactors} can be optionally applied in to
#'   the \code{DESeqDataSet} object in this function, or they can be applied
#'   later on, either by the user or in a call to \code{getDESeqResults}.
#'   Applying the \code{sizeFactors} later is useful if multiple sets of factors
#'   will be explored, although \code{sizeFactors} can be overwritten at any
#'   time. Note that DESeq2 \code{sizeFactors} are \emph{not} the same as
#'   normalization factors defined elsewhere in this package. See details below.
#' @param field Argument passed to \code{getCountsByRegions}. Can be used to
#'   specify fields in a single multiplexed GRanges object, or individual fields
#'   for each GRanges object in \code{dataset.list}.
#' @param blacklist An optional GRanges object containing regions that should be
#'   excluded from signal counting. Use of this argument is distinct from the
#'   use of non-contiguous gene regions (see details below), and the two can be
#'   used simultaneously. Blacklisting doesn't affect the ranges returned as
#'   rowRanges in the output DESeqDataSet object (unlike the use of
#'   non-contiguous regions).
#' @param expand_ranges Logical indicating if ranges in \code{dataset.gr} should
#'   be treated as descriptions of single molecules (\code{FALSE}), or if ranges
#'   should be treated as representing multiple adjacent positions with the same
#'   signal (\code{TRUE}). See \code{\link[BRGenomics:getCountsByRegions]{
#'   getCountsByRegions}}.
#' @param ncores Number of cores to use for read counting across all samples. By
#'   default, all available cores are used.
#' @param quiet If \code{TRUE}, all output messages from call to
#'   \code{\link[DESeq2:DESeqDataSet]{DESeqDataSet}} will be suppressed.
#'
#' @return A \code{DESeqData} object in which \code{rowData} are given as
#'   \code{rowRanges}, which are equivalent to \code{regions.gr}, unless there
#'   are non-contiguous gene regions (see note below). Samples (as seen in
#'   \code{colData}) are factored so that samples are grouped by
#'   \code{replicate} and \code{condition}, i.e. all non-replicate samples are
#'   treated as distinct, and the DESeq2 design = \code{~condition}.
#'
#' @section Use of non-contiguous gene regions: In DESeq2, genes must be defined
#'   by single, contiguous chromosomal locations. In contrast, this function
#'   allows individual genes to be encompassed by multiple distinct ranges in
#'   \code{regions.gr}. To use non-contiguous gene regions, provide
#'   \code{gene_names} in which some names are duplicated. For each unique gene
#'   in \code{gene_names}, this function will generate counts across all ranges
#'   for that gene, but be aware that it will only keep the largest range for
#'   each gene in the resulting \code{DESeqDataSet} object's \code{rowRanges}.
#'   If the desired use is to blacklist certain sites in a genelist, note that
#'   the \code{blacklist} argument can be used.
#'
#' @section A note on DESeq2 sizeFactors: DESeq2 \code{sizeFactors} are
#'   sample-specific normalization factors that are applied by division, i.e.
#'   \eqn{counts_{norm,i}=counts_i / sizeFactor_i}{normcounts_i = counts_i /
#'   sizeFactor_i}. This is in contrast to normalization factors as defined in
#'   this package (and commonly elsewhere), which are applied by multiplication.
#'   Also note that DESeq2's "\code{normalizationFactors}" are not sample
#'   specific, but rather gene specific factors used to correct for
#'   ascertainment bias across different genes (e.g. as might be relevant for
#'   GSEA or Go analysis).
#'
#' @section On gene names and unexpected errors: Certain gene names can cause
#' this function to return an error. We've never encountered errors using
#' conventional, systematic naming schemes (e.g. ensembl IDs), but we have
#' seen errors when using Drosophila (Flybase) "symbols". We expect this is due
#' to the unconventional use of non-alphanumeric characters in some Drosophila
#' gene names.
#'
#' @export
#' @importFrom DESeq2 DESeqDataSet sizeFactors<-
#'
#' @author Mike DeBerardine
#' @seealso \code{\link[DESeq2:DESeqDataSet]{DESeq2::DESeqDataSet}},
#'   \code{\link[BRGenomics:getDESeqResults]{getDESeqResults}}
#'
#' @examples
#' suppressPackageStartupMessages(require(DESeq2))
#' data("PROseq") # import included PROseq data
#' data("txs_dm6_chr4") # import included transcripts
#'
#' # divide PROseq data into 6 toy datasets
#' ps_a_rep1 <- PROseq[seq(1, length(PROseq), 6)]
#' ps_b_rep1 <- PROseq[seq(2, length(PROseq), 6)]
#' ps_c_rep1 <- PROseq[seq(3, length(PROseq), 6)]
#'
#' ps_a_rep2 <- PROseq[seq(4, length(PROseq), 6)]
#' ps_b_rep2 <- PROseq[seq(5, length(PROseq), 6)]
#' ps_c_rep2 <- PROseq[seq(6, length(PROseq), 6)]
#'
#' ps_list <- list(A_rep1 = ps_a_rep1, A_rep2 = ps_a_rep2,
#'                 B_rep1 = ps_b_rep1, B_rep2 = ps_b_rep2,
#'                 C_rep1 = ps_c_rep1, C_rep2 = ps_c_rep2)
#'
#' # make flawed dataset (ranges in txs_dm6_chr4 not disjoint)
#' #    this means there is double-counting
#' # also using discontinuous gene regions, as gene_ids are repeated
#' dds <- getDESeqDataSet(ps_list,
#'                        txs_dm6_chr4,
#'                        gene_names = txs_dm6_chr4$gene_id,
#'                        quiet = TRUE,
#'                        ncores = 1)
#' dds
getDESeqDataSet <- function(dataset.list, regions.gr, sample_names = NULL,
                            gene_names = NULL, sizeFactors = NULL,
                            field = "score", blacklist = NULL,
                            expand_ranges = FALSE,
                            ncores = getOption("mc.cores", 2L), quiet = FALSE) {

    # Get counts dataframe for all samples in each range of regions.gr
    counts.df <- getCountsByRegions(dataset.list, regions.gr, field = field,
                                    blacklist = blacklist,
                                    expand_ranges = expand_ranges,
                                    ncores = ncores)

    # check for valid sample_names
    if (is.null(sample_names))
        sample_names <- names(counts.df)
    .checkdsnames(ncol(counts.df), sample_names, length(regions.gr), gene_names)

    # if given, check gene_names match regions.gr, and if they match multiple
    discont.genes <- FALSE
    if (!is.null(gene_names))
        discont.genes <- length(unique(gene_names)) != length(gene_names)

    # Make colData for SummarizedExperiment
    coldat <- data.frame(condition = factor(sub("_rep.*", "", sample_names)),
                         replicate = factor(sub(".*rep", "rep", sample_names)),
                         row.names = sample_names)

    # Make SummarizedExperiment object
    se <- .get_se(counts.df, regions.gr, gene_names, discont.genes, coldat)

    if (quiet) {
        suppressMessages(dds <- DESeqDataSet(se, design = ~condition))
    } else {
        dds <- DESeqDataSet(se, design = ~condition)
    }

    if (!is.null(sizeFactors))
        sizeFactors(dds) <- sizeFactors
    return(dds)
}


.checkdsnames <- function(ns, sample_names, nr, gene_names) {
    # ns = number of samples; nr = number of regions
    if (length(sample_names) != ns)
        stop(.nicemsg("sample_names are required, and a name is required for
                      each element of dataset.list"))

    if (any(!grepl("_rep.", sample_names)))
        stop(.nicemsg("all sample_names must contain strings naming replicates
                      as such: 'rep1', 'rep2', etc."))

    if (!is.null(gene_names) && length(gene_names) != nr)
        stop(.nicemsg("gene_names given are not the same length as regions.gr;
                      gene_names must correspond 1:1 with the ranges in
                      regions.gr"))
}


#' @importFrom parallel mclapply
#' @importFrom GenomicRanges width
#' @importFrom SummarizedExperiment SummarizedExperiment
.get_se <- function(counts.df, regions.gr, gene_names, discont.genes, coldat) {
    if (!discont.genes) {
        rownames(counts.df) <- gene_names

    } else {
        # aggregate counts by gene
        counts.df <- aggregate(counts.df, by = list(gene_names), FUN = sum)
        rownames(counts.df) <- counts.df[, 1]
        counts.df <- counts.df[, -1] # alphabatized by aggregate

        # for rowRanges, use longest range for each gene
        idx.by.width <- order(width(regions.gr), decreasing = TRUE)
        gnames.by.width <- gene_names[idx.by.width] # sort by range width
        idx <- which(!duplicated(gnames.by.width)) # keep only non-duplicates

        gnames.sort <- gnames.by.width[idx]
        regions.gr <- regions.gr[idx.by.width][idx]

        # sort rowRanges & counts.df to match input order
        # -> counts.df already alphabetized, so use input ordering directly
        input_ordering <- rank(unique(gene_names))
        counts.df <- counts.df[input_ordering, ]

        # -> for rowRanges, map to alphabet then to input order
        map_current_to_input <- order(gnames.sort)[input_ordering]
        regions.gr <- regions.gr[map_current_to_input]
    }

    SummarizedExperiment(assays = list(counts = as.matrix(counts.df)),
                         rowRanges = regions.gr, colData = coldat)
}


### ========================================================================= #
### Get DESeq2 Results from DESeqDataSet
### ------------------------------------------------------------------------- #
###


#' Get DESeq2 results using reduced dispersion matrices
#'
#' This function calls \code{\link[DESeq2:DESeq]{DESeq2::DESeq}} and
#' \code{\link[DESeq2:results]{DESeq2::results}} on a pre-existing
#' \code{DESeqDataSet} object and returns a \code{DESeqResults} table for one or
#' more pairwise comparisons. However, unlike a standard call to
#' \code{DESeq2::results} using the \code{contrast} argument, this function
#' subsets the dataset so that DESeq2 only estimates dispersion for the samples
#' being compared, and not for all samples present.
#'
#' @param dds A DESeqDataSet object, produced using either
#'   \code{\link[BRGenomics:getDESeqDataSet]{getDESeqDataSet}} from this package
#'   or \code{\link[DESeq2:DESeqDataSet]{DESeqDataSet}} from \code{DESeq2}. If
#'   \code{dds} was not created using \code{getDESeqDataSet}, \code{dds} must be
#'   made with \code{design = ~condition} such that a unique \code{condition}
#'   level exists for each sample/treatment condition.
#' @param contrast.numer A string naming the \code{condition} to use as the
#'   numerator in the DESeq2 comparison, typically the perturbative condition.
#' @param contrast.denom A string naming the \code{condition} to use as the
#'   denominator in the DESeq2 comparison, typically the control condition.
#' @param comparisons As an optional alternative to supplying a single
#'   \code{contrast.numer} and \code{contrast.denom}, users can supply a list of
#'   character vectors containing numerator-denominator pairs, e.g.
#'   \code{list(c("B", "A"), c("C", "A"), c("C", "B"))}. \code{comparisons} can
#'   also be a dataframe in which each row is a comparison, the first column
#'   contains the numerators, and  the second column contains the denominators.
#' @param sizeFactors A vector containing DESeq2 \code{sizeFactors} to apply to
#'   each sample. Each sample's readcounts are \emph{divided} by its respective
#'   DESeq2 \code{sizeFactor}. A warning will be generated if the
#'   \code{DESeqDataSet} already contains \code{sizeFactors}, and the previous
#'   \code{sizeFactors} will be over-written.
#' @param alpha The significance threshold passed to \code{DESeqResults}, which
#'   is used for independent filtering of results (see DESeq2 documentation).
#' @param lfcShrink Logical indicating if log2FoldChanges and their standard
#'   errors should be shrunk using \code{\link[DESeq2:lfcShrink]{lfcShrink}}.
#'   LFC shrinkage is very useful for making fold-change values meaningful, as
#'   low-expression/high variance genes are given low fold-changes.
#'   Set to \code{FALSE} by default.
#' @param args.DESeq Additional arguments passed to
#'   \code{\link[DESeq2:DESeq]{DESeq}}, given as a list of argument-value pairs,
#'   e.g. \code{list(fitType = "local", useT = TRUE)}. All arguments given here
#'   will be passed to \code{DESeq} except for \code{object} and
#'   \code{parallel}. If no arguments are given, all defaults will be used.
#' @param args.results Additional arguments passed to
#'   \link[DESeq2:results]{DESeq2::results}, given as a list of argument-value
#'   pairs, e.g. \code{list(altHypothesis = "greater", lfcThreshold = 1.5)}. All
#'   arguments given here will be passed to \code{results} except for
#'   \code{object}, \code{contrast}, \code{alpha}, and \code{parallel}. If no
#'   arguments are given, all defaults will be used.
#' @param args.lfcShrink Additional arguments passed to
#'   \code{\link[DESeq2:lfcShrink]{lfcShrink}}, given as a list of
#'   argument-value pairs. All arguments given here will be passed to
#'   \code{lfcShrink} except for \code{dds}, \code{coef}, \code{contrast}, and
#'   \code{parallel}. If no arguments are given, all defaults will be used.
#' @param ncores The number of cores to use for parallel processing. Multicore
#'   processing is only used if more than one comparison is being made (i.e.
#'   argument \code{comparisons} is used), and the number of cores utilized will
#'   not be greater than the number of comparisons being performed.
#' @param quiet If \code{TRUE}, all output messages from calls to \code{DESeq}
#'   and \code{results} will be suppressed, although passing option \code{quiet}
#'   in \code{args.DESeq} will supersede this option for the call to
#'   \code{DESeq}.
#'
#' @return For a single comparison, the output is the \code{DESeqResults} result
#'   table. If a \code{comparisons} is used to make multiple comparisons, the
#'   output is a named list of \code{DESeqResults} objects, with elements named
#'   following the pattern \code{"X_vs_Y"}, where \code{X} is the name of the
#'   numerator condition, and \code{Y} is the name of the denominator condition.
#'
#' @section Errors when \code{ncores > 1}: If this function returns an error,
#'   set \code{ncores = 1}. Whether or not this occurs can depend on whether
#'   users are using alternative BLAS libraries (e.g. OpenBLAS or Apple's
#'   Accelerate framework) and/or how DESeq2 was installed. This is because some
#'   DESeq2 functions (e.g. \code{\link[DESeq2:nbinomWaldTest]{
#'   nbinomWaldTest}}) use C code that can be compiled to use parallelization,
#'   and this conflicts with our use of process forking (via the
#'   \code{\link[parallel:parallel-package]{parallel package}}) when
#'   \code{ncores > 1}.
#'
#' @author Mike DeBerardine
#' @seealso \code{\link[BRGenomics:getDESeqDataSet]{getDESeqDataSet}},
#'   \code{\link[DESeq2:results]{DESeq2::results}}
#' @export
#' @importFrom DESeq2 sizeFactors sizeFactors<-
#' @importFrom parallel mclapply
#'
#' @examples
#' #--------------------------------------------------#
#' # getDESeqDataSet
#' #--------------------------------------------------#
#' suppressPackageStartupMessages(require(DESeq2))
#' data("PROseq") # import included PROseq data
#' data("txs_dm6_chr4") # import included transcripts
#'
#' # divide PROseq data into 6 toy datasets
#' ps_a_rep1 <- PROseq[seq(1, length(PROseq), 6)]
#' ps_b_rep1 <- PROseq[seq(2, length(PROseq), 6)]
#' ps_c_rep1 <- PROseq[seq(3, length(PROseq), 6)]
#'
#' ps_a_rep2 <- PROseq[seq(4, length(PROseq), 6)]
#' ps_b_rep2 <- PROseq[seq(5, length(PROseq), 6)]
#' ps_c_rep2 <- PROseq[seq(6, length(PROseq), 6)]
#'
#' ps_list <- list(A_rep1 = ps_a_rep1, A_rep2 = ps_a_rep2,
#'                 B_rep1 = ps_b_rep1, B_rep2 = ps_b_rep2,
#'                 C_rep1 = ps_c_rep1, C_rep2 = ps_c_rep2)
#'
#' # make flawed dataset (ranges in txs_dm6_chr4 not disjoint)
#' #    this means there is double-counting
#' # also using discontinuous gene regions, as gene_ids are repeated
#' dds <- getDESeqDataSet(ps_list,
#'                        txs_dm6_chr4,
#'                        gene_names = txs_dm6_chr4$gene_id,
#'                        ncores = 1)
#'
#' dds
#'
#' #--------------------------------------------------#
#' # getDESeqResults
#' #--------------------------------------------------#
#'
#' res <- getDESeqResults(dds, "B", "A")
#'
#' res
#'
#' reslist <- getDESeqResults(dds, comparisons = list(c("B", "A"), c("C", "A")),
#'                            ncores = 1)
#' names(reslist)
#'
#' reslist$B_vs_A
#'
#' # or using a dataframe
#' reslist <- getDESeqResults(dds, comparisons = data.frame(num = c("B", "C"),
#'                                                          den = c("A", "A")),
#'                            ncores = 1)
#' reslist$B_vs_A
getDESeqResults <- function(dds, contrast.numer, contrast.denom,
                            comparisons = NULL, sizeFactors = NULL,
                            alpha = 0.1, lfcShrink = FALSE,
                            args.DESeq = NULL, args.results = NULL,
                            args.lfcShrink = NULL,
                            ncores = getOption("mc.cores", 2L), quiet = FALSE) {

    # check inputs
    comparisons <- .check_args(match.call(), comparisons, quiet)

    ## if length(sizeFactors) matches dds, apply them ("apply early");
    ## else, hold on to them and try to apply after subsetting dds
    when_sf <- .when_sf(dds, sizeFactors) # early, late, or never
    .msgs_early_sf(dds, comparisons, when_sf, quiet)
    if (when_sf == "early") {
        sizeFactors(dds) <- sizeFactors
        sizeFactors <- NULL # prevent re-application
    }

    if (is.null(comparisons)) {
        res <- .get_deseq_results(
            dds, contrast.numer, contrast.denom, sizeFactors = sizeFactors,
            alpha = alpha, lfcShrink = lfcShrink, args.DESeq = args.DESeq,
            args.results = args.results, args.lfcShrink = args.lfcShrink,
            quiet = quiet
        )
        return(res)

    } else {
        args.DESeq <- args.DESeq[names(args.DESeq) != "quiet"]
        results.out <- mclapply(comparisons, function(x) {
            .get_deseq_results(
                dds, x[1], x[2], sizeFactors = sizeFactors, alpha = alpha,
                lfcShrink = lfcShrink, args.DESeq = args.DESeq,
                args.results = args.results, args.lfcShrink = args.lfcShrink,
                quiet = TRUE
            )}, mc.cores = ncores)

        names(results.out) <- vapply(comparisons,
                                     function(x) paste0(x[1], "_vs_", x[2]),
                                     FUN.VALUE = character(1))
        return(results.out)
    }
}


.check_args <- function(args, comparisons, quiet) {
    args <- as.list(args)[-1]
    num <- "contrast.numer" %in% names(args)
    denom <- "contrast.denom" %in% names(args)
    clist <- !is.null(comparisons)

    if (clist) comparisons <- .check_clist(comparisons)

    if (!xor(clist, num & denom))
        stop(.nicemsg("Either provide both contrast.numer and contrast.denom,
                      or provide comparisons, but not both"))
    return(comparisons)
}

.check_clist <- function(comparisons) {
    if (is.data.frame(comparisons)) {
        comparisons <- as.data.frame(t(comparisons), stringsAsFactors = FALSE)
        comparisons <- as.list(comparisons)
    }
    class_ok <- if (!is.list(comparisons)) FALSE else {
        all(vapply(comparisons, is.character, logical(1)))
    }
    lengths_ok <- all(lengths(comparisons) == 2)

    if (!(class_ok & lengths_ok))
        stop(.nicemsg("comparisons provided as input, but it's not a list of
                      length = 2 character vectors, or a dataframe of
                      characters with 2 columns"))
    comparisons
}


.when_sf <- function(dds, sizeFactors) {
    if (is.null(sizeFactors)) return("never")
    if (length(sizeFactors) == nrow(dds@colData)) return("early")
    return("late")
}

.msgs_early_sf <- function(dds, comparisons, when_sf, quiet) {

    already_sf <- !is.null(sizeFactors(dds))

    if (when_sf == "early" && already_sf && !quiet)
        warning("Overwriting previous sizeFactors", immediate. = TRUE)

    if (when_sf == "late" && length(comparisons) > 1)
        stop(message = .nicemsg("Length of sizeFactors not equal to number of
                                samples in dds"))
}


#' @importFrom DESeq2 DESeq results
.get_deseq_results <- function(dds, contrast.numer, contrast.denom, sizeFactors,
                               alpha, lfcShrink, args.DESeq, args.results,
                               args.lfcShrink, quiet) {

    # Subset for pairwise comparison
    dds <- dds[, dds$condition %in% c(contrast.numer, contrast.denom)]
    # drop & sort levels (needed for using apeglm shrinkage)
    dds$condition <- factor(dds$condition, levels = c(contrast.denom,
                                                      contrast.numer))

    # try to apply sizeFactors that weren't the same size as original dds
    dds <- .apply_sf_late(dds, sizeFactors, quiet)

    # =============== Call DESeq2::DESeq() =============== #
    # Get args; only use parent function 'quiet' arg if not in args.DESeq
    args.DESeq <- .merge_args(rqd = expression(object = dds, parallel = FALSE),
                              usr = args.DESeq,
                              exclude = c("object", "parallel"))
    if (!"quiet" %in% names(args.DESeq))
        args.DESeq$quiet <- quiet

    dds <- do.call(DESeq2::DESeq, args.DESeq)

    # ============== Call DESeq2::results() ============== #
    # Get args
    rqd = expression(object = dds, alpha = alpha,
                     contrast = c("condition", contrast.numer, contrast.denom))
    args.results <- .merge_args(rqd = rqd, usr = args.results,
                                exclude = c("object", "contrast",
                                            "alpha", "parallel"))
    if (!quiet) {
        res <- do.call(DESeq2::results, args.results)
    } else {
        res <- suppressWarnings(suppressMessages(
            do.call(DESeq2::results, args.results)
        ))
    }

    # ============= Call DESeq2::lfcShrink() ============= #
    if (lfcShrink) {
        # Get args
        rqd = expression(dds = dds,
                         coef = paste0("condition_", contrast.numer,
                                       "_vs_", contrast.denom),
                         res = res)
        args.lfcShrink <- .merge_args(rqd = rqd, usr = args.lfcShrink,
                                      exclude = c("dds", "coef", "contrast",
                                                  "res", "parallel"))
        if (!"quiet" %in% names(args.lfcShrink))
            args.lfcShrink$quiet <- quiet

        res <- do.call(DESeq2::lfcShrink, args.lfcShrink)
    }
    return(res)
}


#' @importFrom DESeq2 sizeFactors sizeFactors<-
.apply_sf_late <- function(dds, sizeFactors, quiet) {
    when_sf <- .when_sf(dds, sizeFactors)
    already_sf <- !is.null(sizeFactors(dds))

    if (when_sf == "late")
        stop(.nicemsg("Length of sizeFactors not equal to number of samples in
                      dds nor the number of samples in comparison group"))

    if (when_sf == "early") {
        if (already_sf & !quiet)
            warning("Overwriting previous sizeFactors", immediate. = TRUE)
        sizeFactors(dds) <- sizeFactors
    }
    dds
}


.merge_args <- function(rqd, usr, exclude = NULL) {
    # function to combine required args with optional user args
    # exclude is an optional character vector of user args to remove
    if (is.null(usr))
        return(as.list(rqd))

    if (!class(usr) %in% c("list", "expression") || is.null(names(usr)))
        stop(.nicemsg("If given, args.DESeq and args.results must be named
                      lists or R expressions containing argument names and
                      values. See documentation"))
    usr <- as.expression(usr)
    usr <- usr[!names(usr) %in% exclude]
    as.list(c(rqd, usr))
}
