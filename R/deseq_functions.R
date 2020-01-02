### ========================================================================= #
### DESeq2 helper functions
### ------------------------------------------------------------------------- #
###

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
#' @importFrom parallel detectCores
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
getDESeqDataSet <- function(dataset.list, regions.gr,
                            sample_names = names(dataset.list),
                            gene_names = NULL, sizeFactors = NULL,
                            field = "score", ncores = detectCores(),
                            quiet = FALSE) {
    # check for valid sample_names; check gene_names match regions.gr
    .check_snames(dataset.list, sample_names)

    # if given, check gene_names match regions.gr, and if they match multiple
    discont.genes <- FALSE
    if (!is.null(gene_names)) {
        .check_gnames(regions.gr, gene_names)
        discont.genes <- length(unique(gene_names)) != length(gene_names)
    }

    # Make column data (colData) for SummarizedExperiment
    coldat <- .get_coldat(sample_names)

    # Get counts matrix for all samples in each range of regions.gr
    counts.df <- as.data.frame(mclapply(dataset.list, getCountsByRegions,
                                        regions.gr = regions.gr, field = field,
                                        mc.cores = ncores))

    # Make SummarizedExperiment object
    se <- .get_se(counts.df, regions.gr, gene_names, discont.genes, coldat)

    if (quiet) {
        suppressMessages(dds <- DESeqDataSet(se, design = ~condition))
    } else {
        dds <- DESeqDataSet(se, design = ~condition)
    }

    if (!is.null(sizeFactors)) sizeFactors(dds) <- sizeFactors
    return(dds)
}


.check_snames <- function(dataset.list, sample_names) {
    if (length(sample_names) != length(dataset.list)) {
        stop(message = .nicemsg("sample_names are required, and a name is
                                required for each element of dataset.list"))
        return(geterrmessage())
    }

    if (any(!grepl("_rep.", sample_names))) {
        stop(message = .nicemsg("all sample_names must contain strings naming
                                replicates as such: 'rep1', 'rep2', etc."))
        return(geterrmessage())
    }
}

.check_gnames <- function(regions.gr, gene_names) {
    if (length(gene_names) != length(regions.gr)) {
        stop(message = .nicemsg("gene_names given are not the same length
                                as regions.gr; gene_names must correspond
                                1:1 with the ranges in regions.gr"))
        return(geterrmessage())
    }
}

.get_coldat <- function(sample_names) {
    data.frame(condition = sub("_rep.", "", sample_names),
               replicate = sub(".*rep", "rep", sample_names),
               row.names = sample_names)
}

#' @importFrom parallel mclapply
#' @importFrom stats aggregate
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

    SummarizedExperiment(assays = as.matrix(counts.df), rowRanges = regions.gr,
                         colData = coldat)
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
#'   argument \code{comparisons} is used), and the number of cores utilized will
#'   not be greater than the number of comparisons being performed.
#' @param quiet If \code{TRUE}, all output messages from calls to \code{DESeq}
#'   and \code{results} will be suppressed, although passing option \code{quiet}
#'   in \code{args.DESeq} will supersede this option for the call to
#'   \code{DESeq}.
#'
#' @return For a single comparison, the output is the \code{DESeqResults} result
#'   table. If a \code{comparisons} is used to make multiple comparisons,
#'   the output is a named list of \code{DESeqResults} objects, with elements
#'   named following the pattern \code{"X_vs_Y"}, where \code{X} is the name of
#'   the numerator condition, and \code{Y} is the name of the denominator
#'   condition.
#'
#' @author Mike DeBerardine
#' @seealso \code{\link[BRGenomics:getDESeqDataSet]{getDESeqDataSet}},
#'   \code{\link[DESeq2:results]{DESeq2::results}}
#' @export
#' @importFrom DESeq2 sizeFactors sizeFactors<-
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
#'                            comparisons = list(c("B", "A"), c("C", "A")),
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
                            alpha = 0.1, args.DESeq = NULL, args.results = NULL,
                            ncores = detectCores(), quiet = FALSE) {

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
            alpha = alpha, args.DESeq = args.DESeq, args.results = args.results,
            quiet = quiet
        )
        return(res)

    } else {
        args.DESeq <- args.DESeq[names(args.DESeq) != "quiet"]
        results.out <- mclapply(comparisons, function(x) {
            .get_deseq_results(
                dds, x[1], x[2], sizeFactors = sizeFactors, alpha = alpha,
                args.DESeq = args.DESeq, args.results = args.results,
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

    if (!xor(clist, num & denom)) {
        stop(message = .nicemsg("Either provide both contrast.numer and
                                contrast.denom, or provide comparisons,
                                but not both"))
        return(geterrmessage())
    }
    return(comparisons)
}

.check_clist <- function(comparisons) {
    if (is.data.frame(comparisons)) {
        comparisons <- as.data.frame(t(comparisons),
                                     stringsAsFactors = FALSE)
        comparisons <- as.list(comparisons)
    }
    class_ok <- .class_check(comparisons)
    lengths_ok <- all(lengths(comparisons) == 2)

    if (!(class_ok & lengths_ok)) {
        stop(message = .nicemsg("comparisons provided as input, but it's not a
                                list of length = 2 character vectors, or a
                                dataframe of characters with 2 columns"))
        return(geterrmessage())
    }
    return(comparisons)
}

.class_check <- function(comparisons) {
    if (!is.list(comparisons)) return(FALSE)
    classes <- vapply(comparisons, is.character, FUN.VALUE = logical(1))
    if (all(classes)) return(TRUE)
    return(FALSE)
}


.when_sf <- function(dds, sizeFactors) {
    if (is.null(sizeFactors)) return("never")
    if (length(sizeFactors) == nrow(dds@colData)) return("early")
    return("late")
}

.msgs_early_sf <- function(dds, comparisons, when_sf, quiet) {

    already_sf <- !is.null(sizeFactors(dds))

    if (when_sf == "early" & already_sf & !quiet) {
        warning("Overwriting previous sizeFactors", immediate. = TRUE)
    }

    if (when_sf == "late" & length(comparisons) > 1)  {
        stop(message = .nicemsg("Length of sizeFactors not equal to
                                number of samples in dds"))
        return(geterrmessage())
    }
}


#' @importFrom DESeq2 DESeq results
.get_deseq_results <- function(dds, contrast.numer, contrast.denom,
                               sizeFactors, alpha, args.DESeq, args.results,
                               quiet) {

    # Subset for pairwise comparison
    dds <- dds[, dds$condition %in% c(contrast.numer, contrast.denom)]
    dds$condition <- factor(dds$condition) # remove unused levels

    # try to apply sizeFactors that weren't the same size as original dds
    dds <- .apply_sf_late(dds, sizeFactors, quiet)

    # ==== Call DESeq2::DESeq()
    # Get args; only use parent function 'quiet' arg if not in args.DESeq
    args.DESeq <- .merge_args(rqd = expression(object = dds, parallel = FALSE),
                              usr = args.DESeq,
                              exclude = c("object", "parallel"))

    if (!"quiet" %in% names(args.DESeq)) args.DESeq$quiet <- quiet
    dds <- do.call(DESeq2::DESeq, args.DESeq)

    # ==== Call DESeq2::results()
    # Get args
    rqd = expression(object = dds, alpha = alpha,
                     contrast = c("condition", contrast.numer, contrast.denom))
    args.results <- .merge_args(rqd = rqd, usr = args.results,
                                exclude = c("object", "contrast",
                                            "alpha", "parallel"))

    if (!quiet) return( do.call(DESeq2::results, args.results) )
    suppressWarnings(suppressMessages( do.call(DESeq2::results, args.results) ))
}


#' @importFrom DESeq2 sizeFactors sizeFactors<-
.apply_sf_late <- function(dds, sizeFactors, quiet) {
    when_sf <- .when_sf(dds, sizeFactors)
    already_sf <- !is.null(sizeFactors(dds))

    if (when_sf == "late") {
        stop(message = .nicemsg("Length of sizeFactors not equal to number of
                                samples in dds nor the number of samples in
                                comparison group"))
        return(geterrmessage())
    }

    if (when_sf == "early") {
        if (already_sf & !quiet)
            warning("Overwriting previous sizeFactors", immediate. = TRUE)
        sizeFactors(dds) <- sizeFactors
    }

    return(dds)
}


.merge_args <- function(rqd, usr, exclude = NULL) {
    # function to combine required args with optional user args
    # exclude is an optional character vector of user args to remove
    if (is.null(usr))  return(as.list(rqd))

    if (!class(usr) %in% c("list", "expression") |
        is.null(names(usr))) {
        stop(message = .nicemsg("If given, args.DESeq and args.results must be
                                named lists or R expressions containing argument
                                names and values. See documentation"))
        return(geterrmessage())
    }
    usr <- as.expression(usr)
    usr <- usr[!names(usr) %in% exclude]
    return(as.list(c(rqd, usr)))
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
#             all_pairs <- expand.grid(rep(list(seq_along(condition_names)), 2))
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


