#-----------------------------------------------------# #-----------------------------------------------------#
# Functions for working with DESeq2
#-----------------------------------------------------# #-----------------------------------------------------#

# [Note that supplying gene_names can cause an error in call to `DESeqDataSet`.
# I notice it can happen if I try to use symbols, but not for flybase unique IDs].

# in documentation, state that sizeFactors can also be added later

# removed setting names(regions.gr) <- gene_names; replaced with rownames(counts.mat)

# it's possible to have discontinuous genes, if gene_names includes the same gene multiple times
# note that the DESeq2 counts matrix will have the same number of rows as there are unique gene names,
# and there won't be any rowRanges... (?).
#   the DESeq object's rowRanges will contain only the longest range for that gene
#   runs slower if discontinuous genes present

# can speed up this function by using the repeated function (or w/e it's called) that will pull out the
# repeated gene names

md.getDESeqDataSet <- function(dataset.list,
                               regions.gr,
                               sample_names = names(dataset.list),
                               gene_names = NULL,
                               sizeFactors = NULL,
                               ncores = detectCores()) {

    # check for discontinuous genes (single genes with multiple ranges)
    if (!is.null(gene_names) & length(unique(gene_names)) != length(gene_names)) {
        discont.genes <- TRUE
    } else {
        discont.genes <- FALSE
    }

    # Get counts matrix for all samples in regions.gr
    counts.regions <- mclapply(dataset.list,
                               md.get.counts,
                               regions.gr = regions.gr,
                               mc.cores = ncores)

    if (discont.genes) {
        counts.regions <- lapply(counts.regions,
                                 function(x) aggregate(x, by = list(gene_names), FUN = sum)[,2])
        counts.mat <- as.data.frame(counts.regions)
        rownames(counts.mat) <- sort(unique(gene_names)) # aggregate outputs in sorted order
        colnames(counts.mat) <- sample_names
    } else {
        counts.mat <- as.data.frame(counts.regions)
        rownames(counts.mat) <- gene_names
        colnames(counts.mat) <- sample_names
    }

    # Make column data (colData) for SummarizedExperiment
    #   Factors are for each condition, only grouping by replicates

    coldat <- data.frame(treatment = sub("_rep.", "", sample_names),
                         replicate = sub(".*rep", "rep", sample_names), # will be "rep1", "rep2", etc.
                         row.names = sample_names)

    if (discont.genes) {
        # for setting rowRanges, will use the largest range for each gene
        idx <- mclapply(rownames(counts.mat), function(i) {
            idx_i <- which(gene_names == i)
            idx_i[which.max(width(regions.gr[idx_i]))]
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

    suppressMessages(
        counts.dds <- DESeqDataSet(counts.se, design = ~treatment)
    )

    if (!is.null(sizeFactors))  sizeFactors(counts.dds) <- sizeFactors

    return(counts.dds)
}


#-----------------------------------------------------# #-----------------------------------------------------#


# only call DESeq on the subsetted matrix for each pairwise comparison
#  gene-wise dispersion calculated only for the pair!
# results table returned in full, but as per DESeq, the eventual p-value cutoff should
#   be supplied here
# sizeFactors are added before subsetting (so must be for all samples in dds)
#   I should do a check for length(sizeFactors), and apply late if length = number of samples in treatments
#     (number can vary depending on the contrast/number of replicates)
# detectCores only used if comparisons.list given

# comparisons list list(c(contrast.numer, contrast.denom), ...)

# make clear that sizeFactors shouldn't be included twice!
#  issue warning if sizeFactors are defined (or aren't all 1 or whatever it does)

md.getDESeqResults <- function(dds,
                               contrast.numer,
                               contrast.denom,
                               comparisons.list = NULL,
                               sizeFactors = NULL,
                               alpha = 0.05,
                               ncores = detectCores(),
                               quiet = TRUE) {

    if (!is.null(sizeFactors)) {
        if ( !is.null(sizeFactors(dds)) & !(all(sizeFactors(dds) == 1)) ) {
            warning("Overwriting previous sizeFactors", call. = FALSE, immediate. = TRUE)
        }
        sizeFactors(dds) <- sizeFactors
    }

    if (is.null(comparisons.list)) {

        # Run DESeq and return results
        dds <- dds[, dds$treatment %in% c(contrast.numer, contrast.denom)]
        dds$treatment <- factor(dds$treatment) # remove unused levels

        if (quiet) {
            suppressMessages(suppressWarnings(
                dds <- DESeq(dds)
            ))
        } else {
            dds <- DESeq(dds)
        }

        res <- results(dds,
                       contrast = c("treatment", contrast.numer, contrast.denom),
                       alpha = alpha)

        return(res)

    } else {

        # recursive call
        comparisons <- mclapply(comparisons.list,
                                function(x) md.getDESeqResults(dds = dds,
                                                               contrast.numer = x[1],
                                                               contrast.denom = x[2],
                                                               comparisons.list = NULL,
                                                               sizeFactors = NULL,
                                                               alpha = alpha,
                                                               ncores = 1,
                                                               quiet = TRUE),
                                mc.cores = ncores,
                                mc.allow.recursive = FALSE)

        names(comparisons) <- vapply(comparisons.list,
                                     function(x) paste0(x[1], "_vs_", x[2]),
                                     FUN.VALUE = character(1))

        return(comparisons)
    }
}

#-----------------------------------------------------# #-----------------------------------------------------#


md.getDESeqResults.All <- function(dds,
                                   sizeFactors = NULL,
                                   alpha = 0.05,
                                   anchor = NULL,
                                   permutations = FALSE,
                                   additional_comparisons = NULL,
                                   ncores = detectCores()) {

    if (!is.null(sizeFactors))  sizeFactors(dds) <- sizeFactors

    treatment_names <- as.character(levels(dds$treatment))

    # If anchor is given, compare everything to every sample name that pattern matches anchor
    if (!is.null(anchor)) {
        idx.anchor <- grep(anchor, treatment_names)
        all_pairs <- data.frame(Var1 = rep(idx.anchor, each = length(treatment_names)),
                                Var2 = seq_along(treatment_names))
        all_pairs <- subset(all_pairs, Var1 != Var2)
    }

    # If no anchor given, do combinations of comparisons
    if (is.null(anchor)) {
        if (permutations) {
            # all pairwise comparisons, order matters (i.e. get 1-vs-2 and 2-vs-1)
            all_pairs <- expand.grid( rep(list(seq_along(treatment_names)), 2) )
            all_pairs <- subset(all_pairs, Var1 != Var2)
            all_pairs <- all_pairs[, 2:1]
        } else {
            # combinations, order doesn't matter (i.e. if get 1-vs-2, don't get 2-vs-1)
            all_pairs <- combn(seq_along(treatment_names), 2)
            all_pairs <- t(all_pairs)
        }
    }

    # if additional_comparisons given, add them to the end
    if (!is.null(additional_comparisons)) {
        add.numer <- vapply(additional_comparisons, "[[", 1, FUN.VALUE = character(1))
        add.denom <- vapply(additional_comparisons, "[[", 2, FUN.VALUE = character(1))
        idx.add <- data.frame(Var1 = vapply(add.denom,
                                            function(i) which(treatment_names == i),
                                            FUN.VALUE = integer(1)),
                              Var2 = vapply(add.numer,
                                            function(i) which(treatment_names == i),
                                            FUN.VALUE = integer(1)),
                              row.names = NULL)
        all_pairs <- rbind(all_pairs, idx.add)
    }

    # call DESeq2::results for all comparisons
    comparisons <- mclapply(seq_len(nrow(all_pairs)), function(i) {
        numer_i <- all_pairs[i, 2]
        denom_i <- all_pairs[i, 1]
        md.getDESeqResults(dds = dds,
                           contrast.numer = treatment_names[numer_i],
                           contrast.denom = treatment_names[denom_i],
                           alpha = alpha)
    }, mc.cores = ncores)

    names(comparisons) <- vapply(seq_along(comparisons), function(i) {
        numer_i <- all_pairs[i, 2]
        denom_i <- all_pairs[i, 1]
        paste0(treatment_names[numer_i], "_vs_", treatment_names[denom_i])
    }, FUN.VALUE = character(1))

    return(comparisons)
}


