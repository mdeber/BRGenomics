### ========================================================================= #
### Data Import
### ------------------------------------------------------------------------- #
###

#' @importFrom GenomicRanges score score<-
.try_int_score <- function(gr) {
    # if scores are whole numbers, coerce them to integers
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
#'    GenomeInfoDb::standardChromosomes}}
#'
#' @export
#' @importFrom GenomeInfoDb standardChromosomes seqlevels keepSeqlevels
#'   sortSeqlevels
#' @examples
#' # make a GRanges
#' chrom <- c("chr2", "chr3", "chrX", "chrY", "chrM", "junk")
#' gr <- GRanges(seqnames = chrom,
#'               ranges = IRanges(start = 2*(1:6), end = 3*(1:6)),
#'               strand = "+",
#'               seqinfo = Seqinfo(chrom))
#' genome(gr) <- "hg38"
#'
#' gr
#'
#' tidyChromosomes(gr)
#'
#' tidyChromosomes(gr, keep.M = TRUE)
#'
#' tidyChromosomes(gr, keep.M = TRUE, keep.Y = FALSE)
#'
#' tidyChromosomes(gr, keep.nonstandard = TRUE)
tidyChromosomes <- function(gr, keep.X = TRUE, keep.Y = TRUE, keep.M = FALSE,
                            keep.nonstandard = FALSE) {

    chrom <- standardChromosomes(gr)

    if (keep.nonstandard) chrom <- seqlevels(gr)
    if (!keep.X)  chrom <- chrom[ chrom != "chrX" ]
    if (!keep.Y)  chrom <- chrom[ chrom != "chrY" ]
    if (!keep.M)  chrom <- chrom[ (chrom != "chrM") & (chrom != "chrMT") ]
    gr <- keepSeqlevels(gr, chrom, pruning.mode = "tidy")
    sortSeqlevels(gr)
}


#' Import basepair-resolution files
#'
#' @param plus_file,minus_file Paths for strand-specific input files.
#' @param genome Optional string for UCSC reference genome, e.g. "hg38". If
#'   given, non-standard chromosomes are trimmed, and options for sex and
#'   mitochondrial chromosomes are applied.
#' @param keep.X,keep.Y,keep.M,keep.nonstandard Logicals indicating which
#'   non-autosomes should be kept. By default, sex chromosomes are kept, but
#'   mitochondrial and non-standard chromosomes are removed.
#'
#' @details Imports a GRanges object containing base-pair resolution data, with
#'   the \code{score} metadata column indicating the number of reads represented
#'   by each range.
#'
#'   \code{import_bedGraph} is useful for when both 5'- and 3'-end information
#'   is to be maintained for each sequenced molecule. It effectively imports the
#'   entire read.
#'
#'   For \code{import_bigWig}, all ranges are of width = 1.
#'
#' @author Mike DeBerardine
#' @seealso \code{\link[BRGenomics:tidyChromosomes]{tidyChromosomes}},
#'   \code{\link[rtracklayer:import]{rtracklayer::import}}
#' @name import-functions
#' @examples
#' #--------------------------------------------------#
#' # Import PRO-seq bigWigs -> coverage of 3' bases
#' #--------------------------------------------------#
#'
#' # get local address for included bigWig files
#' p.bw <- system.file("extdata", "PROseq_dm6_chr4_plus.bw",
#'                     package = "BRGenomics")
#' m.bw <- system.file("extdata", "PROseq_dm6_chr4_minus.bw",
#'                     package = "BRGenomics")
#'
#' # import bigWigs
#' PROseq <- import_bigWig(p.bw, m.bw, genome = "dm6")
#' PROseq
#'
#' #--------------------------------------------------#
#' # Import PRO-seq bedGraphs -> whole reads (matched 5' and 3' ends)
#' #--------------------------------------------------#
#'
#' # get local address for included bedGraph files
#' p.bg <- system.file("extdata", "PROseq_dm6_chr4_plus.bedGraph",
#'                     package = "BRGenomics")
#' m.bg <- system.file("extdata", "PROseq_dm6_chr4_minus.bedGraph",
#'                     package = "BRGenomics")
#'
#' # import bedGraphs
#' PROseq_paired <- import_bedGraph(p.bg, m.bg, genome = "dm6")
#' PROseq_paired
NULL


#' @rdname import-functions
#' @export
#' @importFrom GenomicRanges score score<- strand<-
#' @importFrom GenomeInfoDb genome<-
import_bigWig <- function(plus_file, minus_file, genome = NULL,
                          keep.X = TRUE, keep.Y = TRUE, keep.M = FALSE,
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
    gr <- makeGRangesBRG(gr)

    if (!is.null(genome)) {
        genome(gr) <- genome
        gr <- tidyChromosomes(gr, keep.X = keep.X, keep.Y = keep.Y,
                              keep.M = keep.M,
                              keep.nonstandard = keep.nonstandard)
    }

    return(sort(gr))
}


#' @rdname import-functions
#' @export
#' @importFrom GenomicRanges score score<- strand<-
#' @importFrom GenomeInfoDb genome<-
import_bedGraph <- function(plus_file, minus_file, genome = NULL,
                            keep.X = TRUE, keep.Y = TRUE, keep.M = FALSE,
                            keep.nonstandard = FALSE) {
    # import bedgraph as GRanges objects
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
        gr <- tidyChromosomes(gr, keep.X = keep.X, keep.Y = keep.Y,
                              keep.M = keep.M,
                              keep.nonstandard = keep.nonstandard)
    }
    return(sort(gr))
}


### convenience function for getting TxDb objects
#
# .getTxDb <- function(genome) {
#     if (genome == "hg38") {
#         db.txs <-
#             TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
#     } else if (genome == "hg19") {
#         db.txs <-
#             TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
#     } else if (genome == "mm10") {
#         db.txs <-
#             TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
#     } else if (genome == "mm9") {
#         db.txs <-
#             TxDb.Mmusculus.UCSC.mm9.knownGene::TxDb.Mmusculus.UCSC.mm9.knownGene
#     } else if (genome == "dm6") {
#         db.txs <-
#             TxDb.Dmelanogaster.UCSC.dm6.ensGene::TxDb.Dmelanogaster.UCSC.dm6.ensGene
#     } else if (genome == "dm3") {
#         db.txs <-
#             TxDb.Dmelanogaster.UCSC.dm3.ensGene::TxDb.Dmelanogaster.UCSC.dm3.ensGene
#     }
# }


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
