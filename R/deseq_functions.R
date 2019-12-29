### ========================================================================= #
### DESeq2 helper functions
### ------------------------------------------------------------------------- #
###

# [Note that supplying gene_names can cause an error in call to `DESeqDataSet`.
# I notice it can happen if I try to use symbols, but not for flybase unique IDs].
#   may have fixed by stopping setting names(regions.gr) <- gene_names;
#   replaced with rownames(counts.mat)

# can speed up this up using the duplicated function that will pull
# out the repeated gene names; so can avoid calling function on whole genelist


#' Get DESeqDataSet objects for downstream analysis
#'
#' This is a convenience function for generating \code{DESeqDataSet} objects,
#' but this function also adds support for counting reads across non-contiguous
#' regions.
#'
#' @param dataset.list A list of GRanges datasets that can be individually
#'   passed to \code{\link[BRGenomics:getCountsByRegions]{getCountsByRegions}}.
#' @param regions.gr A GRanges object containing regions of interest.
#' @param sample_names Names for each dataset in \code{dataset.list} are
#'   required, and by default the names of the list elements are used. The names
#'   must each contain the string "_rep#", where "#" is a single character
#'   (usually a number) indicating the replicate. Sample names across different
#'   replicates must be otherwise identical.
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
#'   time.
#' @param field Argument passed to \code{getCountsByRegions}.
#' @param ncores Number of cores to use for read counting across all samples.
#'   Default is the total number of cores available.
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
#'   by single, contiguous chromosomal locations. This function allows
#'   individual genes to be encompassed by multiple distinct ranges in
#'   \code{regions.gr}. To use non-contiguous gene regions, provide
#'   \code{gene_names} in which some names are duplicated. For each unique gene
#'   in \code{gene_names}, this function will generate counts across all ranges
#'   for that gene, but be aware that it will only keep the largest range for
#'   each gene in the resulting \code{DESeqDataSet} object's \code{rowRanges}.
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
#' @export
#' @importFrom stats aggregate
#' @importFrom parallel detectCores mclapply
#' @importFrom SummarizedExperiment SummarizedExperiment
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
#' ps_list <- list(A_rep1 = ps_a_rep1,
#'                 A_rep2 = ps_a_rep2,
#'                 B_rep1 = ps_b_rep1,
#'                 B_rep2 = ps_b_rep2,
#'                 C_rep1 = ps_c_rep1,
#'                 C_rep2 = ps_c_rep2)
#'
#' # make flawed dataset (ranges in txs_dm6_chr4 not disjoint)
#' #    this means there is double-counting
#' # also using discontinuous gene regions, as gene_ids are repeated
#' dds <- getDESeqDataSet(ps_list,
#'                        txs_dm6_chr4,
#'                        gene_names = txs_dm6_chr4$gene_id,
#'                        quiet = TRUE,
#'                        ncores = 2)
#' dds
getDESeqDataSet <- function(dataset.list, # assumes names end in "_rep#"
                            regions.gr,
                            sample_names = names(dataset.list),
                            gene_names = NULL,
                            sizeFactors = NULL,
                            field = "score",
                            ncores = detectCores(),
                            quiet = FALSE) {

    if (is.null(sample_names)) {
        stop(message = .nicemsg("sample_names are required, but none were
                                found"))
        return(geterrmessage())
    }

    if (any(!grepl("_rep.", sample_names))) {
        stop(message = .nicemsg("all sample_names must contain strings naming
                                replicates as such: 'rep1', 'rep2', etc."))
        return(geterrmessage())
    }

    # check gene_names, and check for non-contiguous genes (>1 range per gene)
    discont.genes <- FALSE
    if (!is.null(gene_names)) {
        if (length(gene_names) != length(regions.gr)) {
            stop(message = .nicemsg("gene_names given are not the same length
                                    as regions.gr; gene_names must correspond
                                    1:1 with the ranges in regions.gr"))
            return(geterrmessage())
        }

        if (length(unique(gene_names)) != length(gene_names))
            discont.genes <- TRUE
    }

    # Get counts matrix for all samples in regions.gr
    counts.regions <- mclapply(dataset.list,
                               getCountsByRegions,
                               regions.gr = regions.gr,
                               field = field,
                               mc.cores = ncores)

    if (discont.genes) {
        counts.regions <- lapply(counts.regions, function(x) {
            aggregate(x, by = list(gene_names), FUN = sum)[,2] })
        counts.mat <- as.data.frame(counts.regions)
        # (aggregate outputs in sorted order)
        rownames(counts.mat) <- sort(unique(gene_names))
        colnames(counts.mat) <- sample_names
    } else {
        counts.mat <- as.data.frame(counts.regions)
        rownames(counts.mat) <- gene_names
        colnames(counts.mat) <- sample_names
    }

    # Make column data (colData) for SummarizedExperiment
    #   Factors are for each condition, only grouping by replicates
    # replicate will be "rep1", "rep2", etc.
    coldat <- data.frame(condition = sub("_rep.", "", sample_names),
                         replicate = sub(".*rep", "rep", sample_names),
                         row.names = sample_names)

    if (discont.genes) {
        # for setting rowRanges, will use the largest range for each gene
        idx <- mclapply(rownames(counts.mat), function(i) {
            idx_i <- which(gene_names == i)
            idx_i[which.max(GenomicRanges::width(regions.gr[idx_i]))]
        }, mc.cores = ncores)
        regions.gr <- regions.gr[unlist(idx)]

        # sort regions.gr to be in input order, and sort counts.mat to match
        reorder.idx <- rank(unique(gene_names))
        regions.gr <- regions.gr[reorder.idx]
        counts.mat <- counts.mat[reorder.idx, ]
    }

    counts.mat <- as.matrix(counts.mat)
    counts.se <- SummarizedExperiment(assays = counts.mat,
                                      rowRanges = regions.gr,
                                      colData = coldat)

    if (quiet) {
        suppressMessages(
            counts.dds <- DESeq2::DESeqDataSet(counts.se, design = ~condition)
        )
    } else {
        counts.dds <- DESeq2::DESeqDataSet(counts.se, design = ~condition)
    }


    if (!is.null(sizeFactors)) DESeq2::sizeFactors(counts.dds) <- sizeFactors

    return(counts.dds)
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
#' @param comparisons.list As an optional alternative to supplying a single
#'   \code{contrast.numer} and \code{contrast.denom}, users can supply a list of
#'   character vectors containing numerator-denominator pairs, e.g.
#'   \code{list(c("B", "A"), c("C", "A"), c("C", "B"))}.
#' @param sizeFactors A vector containing DESeq2 \code{sizeFactors} to apply to
#'   each sample. Each sample's readcounts are \emph{divided} by its respective
#'   DESeq2 \code{sizeFactor}. A warning will be generated if the
#'   \code{DESeqDataSet} already contains \code{sizeFactors}, and the previous
#'   \code{sizeFactors} will be over-written.
#' @param alpha The significance threshold passed to \code{DESeqResults}. This
#'   won't affect the output results, but is used as a performance optimization
#'   by DESeq2.
#' @param args.DESeq Additional arguments passed to
#'   \code{\link[DESeq2:DESeq]{DESeq}}, given as a list of argument-value pairs,
#'   e.g. \code{list(test = "LRT", fitType = "local")}. All arguments given here
#'   will be passed to \code{DESeq} except for \code{object} and
#'   \code{parallel}. If no arguments are given, all defaults will be used.
#' @param args.results Additional arguments passed to
#'   \link[DESeq2:results]{DESeq2::results}, given as a list of argument-value
#'   pairs, e.g. \code{list(altHypothesis = "greater", lfcThreshold = 1.5)}. All
#'   arguments given here will be passed to \code{results} except for
#'   \code{object}, \code{contrast}, \code{alpha}, and \code{parallel}. If no
#'   arguments are given, all defaults will be used.
#' @param ncores The number of cores to use for parallel processing. Multicore
#'   processing is only used if more than one comparison is being made (i.e.
#'   argument \code{comparisons.list} is used), and the number of cores utilized
#'   will not be greater than the number of comparisons being performed.
#' @param quiet If \code{TRUE}, all output messages from calls to \code{DESeq}
#'   and \code{results} will be suppressed, although passing option \code{quiet}
#'   in \code{args.DESeq} will supersede this option for the call to
#'   \code{DESeq}.
#'
#' @return For a single comparison, the output is the \code{DESeqResults} result
#'   table. If a \code{comparisons.list} is used to make multiple comparisons,
#'   the output is a named list of \code{DESeqResults} objects, with elements
#'   named following the pattern \code{"X_vs_Y"}, where \code{X} is the name of
#'   the numerator condition, and \code{Y} is the name of the denominator
#'   condition.
#'
#' @author Mike DeBerardine
#' @seealso \code{\link[BRGenomics:getDESeqDataSet]{getDESeqDataSet}},
#'   \code{\link[DESeq2:results]{DESeq2::results}}
#' @export
#' @importFrom parallel detectCores mclapply
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
#' ps_list <- list(A_rep1 = ps_a_rep1,
#'                 A_rep2 = ps_a_rep2,
#'                 B_rep1 = ps_b_rep1,
#'                 B_rep2 = ps_b_rep2,
#'                 C_rep1 = ps_c_rep1,
#'                 C_rep2 = ps_c_rep2)
#'
#' # make flawed dataset (ranges in txs_dm6_chr4 not disjoint)
#' #    this means there is double-counting
#' # also using discontinuous gene regions, as gene_ids are repeated
#' dds <- getDESeqDataSet(ps_list,
#'                        txs_dm6_chr4,
#'                        gene_names = txs_dm6_chr4$gene_id,
#'                        ncores = 2)
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
#' reslist <- getDESeqResults(dds,
#'                            comparisons.list = list(c("B", "A"), c("C", "A")),
#'                            ncores = 1)
#' names(reslist)
#'
#' reslist[[1]]
getDESeqResults <- function(dds,
                            contrast.numer,
                            contrast.denom,
                            comparisons.list = NULL,
                            sizeFactors = NULL,
                            alpha = 0.1,
                            args.DESeq = NULL,
                            args.results = NULL,
                            ncores = detectCores(),
                            quiet = FALSE) {

    # check only one set of contrast args used
    # and, if given, check format of contrasts.list
    .check_contrast_args(as.list(match.call()[-1]))
    if (!is.null(comparisons.list)) .check_list(comparisons.list)

    # determine if/when to apply sizeFactors, and whether dds already has them
    when_sf <- .when_sf(dds, sizeFactors)
    exist_sf <- .exist_sizeFactors(dds)

    # apply early NFs
    if (when_sf == "early") {
        if (exist_sf & !quiet)
            warning("Overwriting previous sizeFactors", immediate. = TRUE)
        DESeq2::sizeFactors(dds) <- sizeFactors
        sizeFactors <- NULL
    }

    if (is.null(comparisons.list)) {

        res <- .get_deseq_results(
            dds, contrast.numer, contrast.denom, sizeFactors = sizeFactors,
            alpha = alpha, args.DESeq = args.DESeq, args.results = args.results,
            ncores = ncores, quiet = quiet
        )
        return(res)

    } else {

        # if multiple comparisons, return error if trying to apply NFs late
        if (length(comparisons.list) > 1)  {
            if (when_sf == "late") {
                stop(message = .nicemsg("Length of sizeFactors not equal to
                                        number of samples in dds"))
                return(geterrmessage())
            }
            sizeFactors <- NULL
        }

        args.DESeq <- args.DESeq[names(args.DESeq) != "quiet"]
        comparisons <- mclapply(comparisons.list, function(x) {
            .get_deseq_results(
                dds, x[1], x[2], sizeFactors = sizeFactors, alpha = alpha,
                args.DESeq = args.DESeq, args.results = args.results,
                ncores = 1, quiet = TRUE
            )}, mc.cores = ncores)

        names(comparisons) <- vapply(comparisons.list,
                                     function(x) paste0(x[1], "_vs_", x[2]),
                                     FUN.VALUE = character(1))
        return(comparisons)
    }
}

## Helper functions

.check_contrast_args <- function(args) {
    num <- "contrast.numer" %in% names(args)
    denom <- "contrast.denom" %in% names(args)
    clist <- !is.null(args$comparisons.list)

    if (!xor(clist, num & denom)) {
        stop(message = .nicemsg("Either provide both contrast.numer and
                                contrast.denom, or provide comparisons.list,
                                but not both"))
        return(geterrmessage())
    }
}

.check_list <- function(comparisons.list) {
    class_ok <- .class_check(comparisons.list)
    lengths_ok <- all(lengths(comparisons.list) == 2)
    if (!(class_ok & lengths_ok)) {
        stop(message = .nicemsg("comparisons.list provided as input, but
                                it's not a list of length = 2 character
                                vectors"))
        return(geterrmessage())
    }
}

.class_check <- function(comparisons) {
    if (!is.list(comparisons)) return(FALSE)
    classes <- vapply(comparisons, class, FUN.VALUE = character(1))
    if (all(classes == "character")) return(TRUE)
    return(FALSE)
}

.when_sf <- function(dds, sizeFactors) {
    if (is.null(sizeFactors)) return("never")
    if (length(sizeFactors) == nrow(dds@colData)) return("early")
    return("late")
}

#' @importFrom DESeq2 sizeFactors
.exist_sizeFactors <- function(dds) {
    if (!is.null(sizeFactors(dds))) return(TRUE)
    return(FALSE)
}

.merge_args <- function(rqd_args, user_args, exclude = NULL) {
    # rqd_args/user_args are expressions or lists of expressions
    # exclude is an optional character vector of user_args to remove
    if (is.null(user_args))  return(as.list(rqd_args))

    if (!class(user_args) %in% c("list", "expression") |
        is.null(names(user_args))) {
        stop(message = .nicemsg("If given, args.DESeq and args.results must be
                                named lists or R expressions containing argument
                                names and values. See documentation"))
        return(geterrmessage())
    }
    user_args <- as.expression(user_args)
    user_args <- user_args[!names(user_args) %in% exclude]
    return(as.list(c(rqd_args, user_args)))
}

# Function dispatched by getDESeqResults
.get_deseq_results <- function(dds,
                               contrast.numer,
                               contrast.denom,
                               sizeFactors,
                               alpha,
                               args.DESeq,
                               args.results,
                               ncores,
                               quiet) {

    # Subset for pairwise comparison
    dds <- dds[, dds$condition %in% c(contrast.numer, contrast.denom)]
    dds$condition <- factor(dds$condition) # remove unused levels

    # try to apply sizeFactors that weren't the same size as original dds
    when_sf <- .when_sf(dds, sizeFactors)
    exist_sf <- .exist_sizeFactors(dds)

    if (when_sf == "late") {
        stop(message = .nicemsg("Length of sizeFactors not equal to number of
                                samples in dds nor the number of samples in
                                comparison group"))
        return(geterrmessage())
    }

    if (when_sf == "early") {
        if (exist_sf & !quiet)
            warning("Overwriting previous sizeFactors", immediate. = TRUE)
        DESeq2::sizeFactors(dds) <- sizeFactors
    }

    #---------------#
    # Call DESeq()
    #---------------#

    # Get args; only use parent function 'quiet' arg if not in args.DESeq
    args.DESeq <- .merge_args(expression(object = dds, parallel = FALSE),
                              user_args = args.DESeq,
                              exclude = c("object", "parallel"))
    if (!"quiet" %in% names(args.DESeq))
        args.DESeq <- .merge_args(args.DESeq, list(quiet = quiet))

    dds <- do.call(DESeq2::DESeq, args.DESeq)

    #---------------#
    # Call results()
    #---------------#

    # Get args
    args.results <- .merge_args(expression(object = dds,
                                           contrast = c("condition",
                                                        contrast.numer,
                                                        contrast.denom),
                                           alpha = alpha),
                                args.results,
                                exclude = c("object", "contrast",
                                            "alpha", "parallel"))
    if (quiet) {
        suppressWarnings(suppressMessages(
            res <- do.call(DESeq2::results, args.results)
        ))
    } else {
        res <- do.call(DESeq2::results, args.results)
    }
    return(res)
}


### ========================================================================= #
### Get Batches of DESeq2 Results from DESeqDataSet
### ------------------------------------------------------------------------- #
###


# #' Automate batch calls to getDESeqResults
# #'
# #' This function can automate the generation numerous pairwise DESeq2
# #' comparisons using several logical schemes.
# #'
# #'
# #' @param dds
# #' @param sizeFactors
# #' @param alpha
# #' @param anchor
# #' @param permutations
# #' @param additional_comparisons
# #' @param ncores
# #'
# #' @return
# #' @export
# #'
# #' @examples
# getDESeqResultsInBatch <- function(dds,
#                                    sizeFactors = NULL,
#                                    alpha = 0.05,
#                                    anchor = NULL,
#                                    permutations = FALSE,
#                                    additional_comparisons = NULL,
#                                    ncores = detectCores()) {
#
#     if (!is.null(sizeFactors))  sizeFactors(dds) <- sizeFactors
#
#     condition_names <- as.character(levels(dds$condition))
#
#     # If anchor is given, compare everything to every sample name that pattern
#     # matches anchor
#     if (!is.null(anchor)) {
#         idx.anchor <- grep(anchor, condition_names)
#         all_pairs <- data.frame(Var1 = rep(idx.anchor,
#                                            each = length(condition_names)),
#                                 Var2 = seq_along(condition_names))
#         all_pairs <- subset(all_pairs, Var1 != Var2)
#     }
#
#     # If no anchor given, do combinations of comparisons
#     if (is.null(anchor)) {
#         if (permutations) {
#             # all pairwise comparisons, order matters
#             #   (i.e. get 1-vs-2 and 2-vs-1)
#             all_pairs <- expand.grid( rep(list(seq_along(condition_names)), 2) )
#             all_pairs <- subset(all_pairs, Var1 != Var2)
#             all_pairs <- all_pairs[, 2:1]
#         } else {
#             # combinations, order doesn't matter
#             #   (i.e. if get 1-vs-2, don't get 2-vs-1)
#             all_pairs <- combn(seq_along(condition_names), 2)
#             all_pairs <- t(all_pairs)
#         }
#     }
#
#     # if additional_comparisons given, add them to the end
#     if (!is.null(additional_comparisons)) {
#         add.numer <- vapply(additional_comparisons, "[[", 1, character(1))
#         add.denom <- vapply(additional_comparisons, "[[", 2, character(1))
#         idx.add <- data.frame(
#             Var1 = vapply(add.denom,
#                           function(i) which(condition_names == i),
#                           FUN.VALUE = integer(1)),
#             Var2 = vapply(add.numer,
#                           function(i) which(condition_names == i),
#                           FUN.VALUE = integer(1)),
#             row.names = NULL
#         )
#         all_pairs <- rbind(all_pairs, idx.add)
#     }
#
#     # call DESeq2::results for all comparisons
#     comparisons <- mclapply(seq_len(nrow(all_pairs)), function(i) {
#         numer_i <- all_pairs[i, 2]
#         denom_i <- all_pairs[i, 1]
#         md.getDESeqResults(dds = dds,
#                            contrast.numer = condition_names[numer_i],
#                            contrast.denom = condition_names[denom_i],
#                            alpha = alpha)
#     }, mc.cores = ncores)
#
#     names(comparisons) <- vapply(seq_along(comparisons), function(i) {
#         numer_i <- all_pairs[i, 2]
#         denom_i <- all_pairs[i, 1]
#         paste0(condition_names[numer_i], "_vs_", condition_names[denom_i])
#     }, FUN.VALUE = character(1))
#
#     return(comparisons)
# }


