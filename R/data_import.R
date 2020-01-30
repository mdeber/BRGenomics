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
#' chromosomes from any GRanges object.
#'
#' @param gr Any GRanges object, however the object should have a standard
#'   genome set, e.g. \code{genome(gr) <- "hg38"}
#' @param keep.X,keep.Y,keep.M,keep.nonstandard Logicals indicating which
#'   non-autosomes should be kept. By default, sex chromosomes are kept, but
#'   mitochondrial and non-standard chromosomes are removed.
#'
#' @return A GRanges object in which both ranges and \code{seqinfo} associated
#'   with trimmed chromosomes have been removed.
#'
#' @details Standard chromosomes are defined using the
#' \code{\link[GenomeInfoDb:standardChromosomes]{standardChromosomes}} function
#' from the \code{GenomeInfoDb} package.
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
#' Import functions for plus/minus pairs of \code{bigWig} or \code{bedGraph}
#' files.
#'
#' @param plus_file,minus_file Paths for strand-specific input files.
#' @param genome Optional string for UCSC reference genome, e.g. "hg38". If
#'   given, non-standard chromosomes are trimmed, and options for sex and
#'   mitochondrial chromosomes are applied.
#' @param keep.X,keep.Y,keep.M,keep.nonstandard Logicals indicating which
#'   non-autosomes should be kept. By default, sex chromosomes are kept, but
#'   mitochondrial and non-standard chromosomes are removed.
#'
#' @return Imports a GRanges object containing base-pair resolution data, with
#'   the \code{score} metadata column indicating the number of reads represented
#'   by each range.
#'
#' @details For \code{import_bigWig}, the output GRanges is formatted by
#'   \code{\link[BRGenomics:makeGRangesBRG]{makeGRangesBRG}}, such that all
#'   ranges are disjoint and have width = 1, and the \code{score} is single-base
#'   coverage, i.e. the number of reads for each position.
#'
#'   \code{import_bedGraph} is useful for when both 5'- and 3'-end information
#'   is to be maintained for each sequenced molecule. It effectively imports the
#'   entire read, and the \code{score} represents the number of reads sharing
#'   identical 5' and 3' ends.
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


#' Import bam files
#'
#' Import single-end or paired-end bam files as GRanges objects, with various
#' processing options.
#'
#' @param file Path of a bam file.
#' @param mapq Filter reads by a minimum MAPQ score. This is the correct way to
#'   filter multi-aligners.
#' @param revcomp Logical indicating if aligned reads should be
#'   reverse-complemented.
#' @param shift Either an integer giving the number of bases by which to shift
#'   the entire read upstream or downstream, or a pair of integers indicating
#'   shifts to be applied to the 5' and 3' ends of the reads, respectively.
#'   Shifting is strand-specific, with negative numbers shifting the reads
#'   upstream, and positive numbers shiftem them downstream. This option is
#'   applied \emph{after} the \code{revcomp}, but before \code{trim.to} and
#'   \code{ignore.strand} options are applied.
#' @param trim.to Option for selecting specific bases from the reads, applied
#'   after the \code{revcomp} and \code{shift} options. By default, the entire
#'   read is maintained. Other options are to take only the 5' base, only the 3'
#'   base, or the only the center base of the read.
#' @param ignore.strand Logical indicating if the strand information should be
#'   discarded. If \code{TRUE}, strand information is discarded \emph{after}
#'   \code{revcomp}, \code{trim.to}, and \code{shift} options are applied.
#' @param field Metadata field name to use for readcounts, usually "score". If
#'   set to \code{NULL}, identical reads (or identical positions if
#'   \code{trim.to} options applied) are not combined, and the length of the
#'   output GRanges will be equal to the number of input reads.
#' @param paired_end Logical indicating if reads should be treated as paired end
#'   reads. When set to \code{NULL} (the default), a test is performed to
#'   determine if the bam file contains any paired-end reads.
#' @param yieldSize The number of bam file records to process simultaneously,
#'   e.g. the "chunk size". Setting a higher chunk size will use more memory,
#'   which can increase speed if there is enough memory available. If chunking
#'   is not necessary, set to \code{NA}.
#'
#' @details If function produces an error, make the \code{paired_end} parameter
#'   explicit, i.e. \code{TRUE} or \code{FALSE}.
#'
#' @return A GRanges object.
#' @author Mike DeBerardine & Nate Tippens
#' @export
#' @importFrom Rsamtools BamFile testPairedEndBam ScanBamParam
#' @importFrom GenomicAlignments readGAlignments readGAlignmentPairs
#' @importFrom GenomicFiles reduceByYield
#' @importFrom GenomicRanges GRanges strand resize shift
#'
#' @examples
#' # get local address for included bam file
#' ps.bam <- system.file("extdata", "PROseq_dm6_chr4.bam",
#'                       package = "BRGenomics")
#'
#' #--------------------------------------------------#
#' # Import entire reads
#' #--------------------------------------------------#
#'
#' # Note that PRO-seq reads are sequenced as reverse complement
#' import_bam(ps.bam, revcomp = TRUE, paired_end = FALSE)
#'
#' #--------------------------------------------------#
#' # Import entire reads, 1 range per read
#' #--------------------------------------------------#
#'
#' import_bam(ps.bam, revcomp = TRUE, field = NULL,
#'            paired_end = FALSE)
#'
#' #--------------------------------------------------#
#' # Import PRO-seq reads at basepair-resolution
#' #--------------------------------------------------#
#'
#' # the typical manner to import PRO-seq data:
#' import_bam(ps.bam, revcomp = TRUE, trim.to = "3p",
#'            paired_end = FALSE)
#'
#' #--------------------------------------------------#
#' # Import PRO-seq reads, removing the run-on base
#' #--------------------------------------------------#
#'
#' # the best way to import PRO-seq data; removes the
#' # most 3' base, which was added in the run-on
#' import_bam(ps.bam, revcomp = TRUE, trim.to = "3p",
#'            shift = -1, paired_end = FALSE)
#'
#' #--------------------------------------------------#
#' # Import 5' ends of PRO-seq reads
#' #--------------------------------------------------#
#'
#' # will include bona fide TSSes as well as hydrolysis products
#' import_bam(ps.bam, revcomp = TRUE, trim.to = "5p",
#'            paired_end = FALSE)
import_bam <- function(file, mapq = 20, revcomp = FALSE, shift = 0L,
                       trim.to = c("whole", "5p", "3p", "center"),
                       ignore.strand = FALSE, field = "score",
                       paired_end = NULL, yieldSize = 2.5e5) {
    trim.to <- match.arg(trim.to, c("whole", "5p", "3p", "center"))

    # Load bam file
    bf <- BamFile(file, yieldSize = yieldSize)
    if (is.null(paired_end)) paired_end <- testPairedEndBam(bf)
    param <- ScanBamParam(mapqFilter = mapq)

    yield_fun <- .get_yield_fun(paired_end, param)

    gr <- GenomicFiles::reduceByYield(X = bf, YIELD = yield_fun, REDUCE = c)

    # Apply Options
    if (revcomp | shift != 0)  is_plus <- as.character(strand(gr)) == "+"
    if (revcomp) {
        strand(gr) <- "+"
        strand(gr)[is_plus] <- "-"
        is_plus <- !is_plus
    }
    if (shift != 0)  gr <- .shift_gr(gr, is_plus, shift)
    if (trim.to != "whole") {
        opt <- paste0("opt.", trim.to)
        opt.arg <- list(opt.5p = "start", opt.3p = "end", opt.center = "center")
        gr <- resize(gr, width = 1, fix = opt.arg[[opt]])
    }
    if (ignore.strand)  strand(gr) <- "*"

    gr <- sort(gr)
    if (!is.null(field))  gr <- .collapse_reads(gr, field)
    return(gr)
}

.get_yield_fun <- function(paired, param) {
    if (paired) {
        function(x) GRanges(readGAlignmentPairs(x, use.names = FALSE,
                                                param = param))
    } else {
        function(x) GRanges(readGAlignments(x, use.names = FALSE,
                                            param = param))
    }
}

.shift_gr <- function(gr, is_plus, shift) {
    if (length(shift) == 1) {
        shifts <- rep(-shift, length(is_plus)) # minus strand
        shifts[is_plus] <- shift
        return(shift(gr, shifts))
    } else if (length(shift) == 2) {
        return(genebodies(gr, shift[1], shift[2], min.window = 0))
    } else {
        stop(message = .nicemsg("shift argument must be a numeric, or a numeric
                                vector of length 2"))
        return(geterrmessage())
    }
}

#' @importFrom GenomicRanges countOverlaps mcols
.collapse_reads <- function(gr, field) {
    gr.out <- unique(gr)
    mcols(gr.out)[field] <- countOverlaps(gr.out, gr, type = "equal")
    return(gr.out)
}



# ---------------------- #


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
#         gene_names <- AnnotationDbi::select(
#             db.txs, keys = keys(db.txs),
#             columns = c("TXNAME"), keytype = "GENEID"
#         )
#     )
#     txs <- txs[order(txs$tx_name)] # pre-sort for speed
#     gene_names <- gene_names[order(gene_names$TXNAME), ]
#
#     if ( all(txs$tx_name == gene_names$TXNAME) )
#         txs$gene_id <- gene_names$GENEID
#
#     # remove non-standard & mitochondrial chromosomes
#     txs <- tidyChromosomes(txs, keep.X = keep.X, keep.Y = keep.Y,
#                            keep.M = keep.M,
#                            keep.nonstandard = keep.nonstandard)
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
