### ========================================================================= #
### Data Import
### ------------------------------------------------------------------------- #
###

# if scores are whole numbers, coerce them to integers
.try_int_score <- function(gr) {
    if (all( round(score(gr) %% 1, 3) == 0 ))
        score(gr) <- as.integer(score(gr))
    gr
}



#' Remove odd chromosomes from GRanges objects
#'
#' This convenience function removes non-standard, mitochondrial, and/or sex
#' chromosomes from any GRanges object. For the chromosomes being removed, any
#' ranges found on those chromosomes are removed, and the chromosomes are also
#' removed from \code{seqinfo}. Standard chromosomes are defined using the
#' \code{\link[GenomeInfoDb:standardChromosomes]{standardChromosomes}} function
#' from the \code{GenomeInfoDb} package.
#'
#' @param gr Any GRanges object, however the object should have a standard
#'   genome set, e.g. \code{genome(gr) <- "hg38"}
#' @param keep.X,keep.Y,keep.M,keep.nonstandard Logicals indicating which
#'   non-autosomes should be kept. By default, sex chromosomes are kept, but
#'   mitochondrial and non-standard chromosomes are removed.
#'
#' @author Mike DeBerardine
#' @seealso
#'    \code{\link[GenomeInfoDb:standardChromosomes]{
#'    GenomeInfoDb::standardChromosomes}},
#'
#' @export
tidyChromosomes <- function(gr,
                            keep.X = TRUE,
                            keep.Y = TRUE,
                            keep.M = FALSE,
                            keep.nonstandard = FALSE) {

    chrom <- standardChromosomes(gr)

    if (keep.nonstandard) chrom <- seqlevels(gr)
    if (!keep.X)  chrom <- chrom[ chrom != "chrX" ]
    if (!keep.Y)  chrom <- chrom[ chrom != "chrY" ]
    if (!keep.M)  chrom <- chrom[ (chrom != "chrM") & (chrom != "chrMT") ]
    gr <- keepSeqlevels(gr, chrom, pruning.mode = "tidy")
    sortSeqlevels(gr)
}



#' Import PRO-seq (or similar) bigWig files
#'
#' This function imports plus/minus pairs of bigWig files containing
#' basepair-resolution data, e.g. PRO-seq or PRO-cap data.
#'
#' @param plus_bw Path of plus strand bigWig file.
#' @param minus_bw Path of minus strand bigWig file.
#' @param genome Optional string for UCSC reference genome, e.g. "hg38". If
#'   given, non-standard chromosomes are trimmed, and options for sex and
#'   mitochondrial chromosomes are applied.
#' @param keep.X Logical indicating whether the X chromosome should be kept.
#' @param keep.Y Logical indicating whether the Y chromosome should be kept.
#' @param keep.M Logical indicating whether mitochondrial chromosomes should be
#'   kept.
#'
#' @return Imports a GRanges object containing base-pair resolution data, with
#'   the \code{score} metadata column indicating readcounts at each base. All
#'   ranges are of width = 1.
#' @author Mike DeBerardine
#' @export
import.PROseq <- function(plus_file,
                          minus_file,
                          genome = NULL,
                          keep.X = TRUE,
                          keep.Y = TRUE,
                          keep.M = FALSE,
                          keep.nonstandard = FALSE) {

    # make possible to import only plus or minus

    # import bw as GRanges objects
    p_gr <- rtracklayer::import.bw(plus_file)
    m_gr <- rtracklayer::import.bw(minus_file)
    score(m_gr) <- abs(score(m_gr)) # scores = reads; make all positive

    strand(p_gr) <- "+"
    strand(m_gr) <- "-"
    suppressWarnings( gr <- c(p_gr, m_gr) )

    # scores are imported as doubles by default; if whole numbers, make integers
    gr <- .try_int_score(gr)

    # Make the width of each range equal to 1
    gr <- makeGRangesBPres(gr)

    if (!is.null(genome)) {
        genome(gr) <- genome
        gr <- tidyChromosomes(gr,
                              keep.X = keep.X,
                              keep.Y = keep.Y,
                              keep.M = keep.M,
                              keep.nonstandard = keep.nonstandard)
    }

    return(sort(gr))
}



#' Import CoPRO (or similar) bedGraph files
#'
#' This function imports plus/minus pairs of bedGraph files. This function is
#' useful for when both 5'- and 3'-end information is to be maintained for each
#' sequenced molecule.
#'
#' @param plus.file Path of plus strand bedGraph file.
#' @param minus.file Path of minus strand bedGraph file.
#' @param genome Optional string for UCSC reference genome, e.g. "hg38". If
#'   given, non-standard chromosomes are trimmed.
#' @param keep.X Logical indicating whether the X chromosome should be kept.
#' @param keep.Y Logical indicating whether the Y chromosome should be kept.
#' @param keep.M Logical indicating whether mitochondrial chromosomes should be
#'   kept.
#'
#' @return Imports a GRanges object containing entire strand-specific reads.
#'   Each range is unique, and the \code{score} metadata column indicates the
#'   number of identical reads (which share the same 5' and 3' ends).
#' @author Mike DeBerardine
#' @export
import.CoPRO <- function(plus_file,
                         minus_file,
                         genome = NULL,
                         keep.X = TRUE,
                         keep.Y = TRUE,
                         keep.M = FALSE,
                         keep.nonstandard = FALSE) {
    # import bw as GRanges objects
    p_bg <- rtracklayer::import.bedGraph(plus_file)
    m_bg <- rtracklayer::import.bedGraph(minus_file)
    score(m_bg) <- abs(score(m_bg)) # scores = reads; make all positive

    strand(p_bg) <- "+"
    strand(m_bg) <- "-"
    suppressWarnings(gr <- c(p_bg, m_bg)) # combine into 1 GRanges object

    # scores are imported as doubles by default; if whole numbers, make integers
    gr <- .try_int_score(gr)

    if (!is.null(genome)) {
        genome(gr) <- genome
        gr <- tidyChromosomes(gr,
                              keep.X = keep.X,
                              keep.Y = keep.Y,
                              keep.M = keep.M,
                              keep.nonstandard = keep.nonstandard)
    }
    return(sort(gr))
}


#' Import bigWig files (general)
#'
#' General function for importing a single bigWig file as a GRanges object. The
#' added functionality over \code{rtracklayer::import.bw} is in trimming odd
#' chromosomes.
#'
#' @param file Path of a bigWig file (non-stranded).
#' @param genome Optional string for UCSC reference genome, e.g. "hg38". If
#'   given, non-standard chromosomes are trimmed.
#' @param keep.X Logical indicating whether the X chromosome should be kept.
#' @param keep.Y Logical indicating whether the Y chromosome should be kept.
#' @param keep.M Logical indicating whether mitochondrial chromosomes should be
#'   kept.
#'
#' @return Imports a GRanges object
#' @author Mike DeBerardine
#' @export
import.bw_trim <- function(file,
                           genome = NULL,
                           keep.X = TRUE,
                           keep.Y = TRUE,
                           keep.M = FALSE,
                           keep.nonstandard = FALSE) {

    gr <- rtracklayer::import.bw(file)

    # scores are imported as doubles by default; if whole numbers, make integers
    gr <- .try_int_score(gr)

    if (!is.null(genome)) {
        genome(gr) <- genome
        gr <- tidyChromosomes(gr,
                              keep.X = keep.X,
                              keep.Y = keep.Y,
                              keep.M = keep.M,
                              keep.nonstandard = keep.nonstandard)
    }

    return(sort(gr))
}


# #' Import transcripts from UCSC
# #'
# #' Imports all annotated transcripts by calling
# #' \code{\link[GenomicFeatures:transcripts]{GenomicFeatures::transcripts}, which
# #' includes alternative isoforms. Currently supports hg38, hg19, mm10, mm9, dm6,
# #' and dm3. This script filters non-standard and mitochondrial chromosomes.
# #'
# #' @param genome A string indicating the UCSC reference genome, e.g. "hg38".
# #' @param keep.X A logical indicating if X chromosome transcripts should be
# #'   kept.
# #' @param keep.Y A logical indicating if Y chromosome transcripts should be
# #'   kept.
# #'
# #' @return A GRanges object, including metadata for unique identifiers.
# #' @author Mike DeBerardine
# #' @export
# importTxsUCSC <- function(genome,
#                           keep.X = TRUE,
#                           keep.Y = TRUE,
#                           keep.M = FALSE,
#                           keep.nonstandard = FALSE) {
#
#     db.txs <- .getTxDb(genome)
#     txs <- GenomicFeatures::transcripts(db.txs)
#
#     # keep tx_name (unique ids); remove tx_id (numbered starting at 1);
#     mcols(txs) <- data.frame(tx_name = txs$tx_name, stringsAsFactors = F)
#
#     # add gene_id
#     suppressMessages(
#         gene_names <- AnnotationDbi::select(db.txs,
#                                            keys = keys(db.txs),
#                                            columns = c("TXNAME"),
#                                            keytype = "GENEID")
#     )
#     txs <- txs[order(txs$tx_name)] # pre-sort for speed
#     gene_names <- gene_names[order(gene_names$TXNAME), ]
#
#     if ( all(txs$tx_name == gene_names$TXNAME) )
#         txs$gene_id <- gene_names$GENEID
#
#     # remove non-standard & mitochondrial chromosomes
#     txs <- tidyChromosomes(txs,
#                        keep.X = keep.X,
#                        keep.Y = keep.Y,
#                        keep.M = keep.M,
#                        keep.nonstandard = keep.nonstandard)
#
#     return(sort(txs))
# }
#
#
#
# #' Import genes from UCSC
# #'
# #' Imports all annotated genes by calling
# #' \code{\link[GenomicFeatures:genes]{GenomicFeatures::genes}}, which provides a
# #' single range for each annotated gene. Currently supports hg38, hg19, mm10,
# #' mm9, dm6, and dm3. This script filters non-standard and mitochondrial
# #' chromosomes.
# #'
# #' @param genome A string indicating the UCSC reference genome, e.g. "hg38".
# #' @param keep.X A logical indicating if X chromosome transcripts should be
# #'   kept.
# #' @param keep.Y A logical indicating if Y chromosome transcripts should be
# #'   kept.
# #'
# #' @return A GRanges object, including metadata for unique identifiers.
# #' @author Mike DeBerardine
# #' @export
# importGenesUCSC <- function(genome,
#                             keep.X = TRUE,
#                             keep.Y = TRUE,
#                             keep.M = FALSE,
#                             keep.nonstandard = FALSE) {
#
#     db.txs <- .getTxDb(genome)
#     genes <- GenomicFeatures::genes(db.txs)
#     genes <- tidyChromosomes(txs,
#                          keep.X = keep.X,
#                          keep.Y = keep.Y,
#                          keep.M = keep.M,
#                          keep.nonstandard = keep.nonstandard)
#     return(sort(genes))
# }
#
