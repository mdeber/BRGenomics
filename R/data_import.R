### ========================================================================= #
### Data Import
### ------------------------------------------------------------------------- #
###

#' @importFrom GenomicRanges score score<-
.try_int_score <- function(gr) {
    # if scores are whole numbers, coerce them to integers
    if (.close_int(score(gr)))
        score(gr) <- as.integer(score(gr))
    gr
}


#' Remove odd chromosomes from GRanges objects
#'
#' This convenience function removes non-standard, mitochondrial, and/or sex
#' chromosomes from any GRanges object.
#'
#' @param gr Any GRanges object, or any another object with associated
#'   \code{seqinfo} (or a \code{Seqinfo} object itself). The object should
#'   typically have a standard genome associated with it, e.g. \code{genome(gr)
#'   <- "hg38"}. \code{gr} can also be a list of such GRanges objects.
#' @param keep.X,keep.Y,keep.M,keep.nonstandard Logicals indicating which
#'   non-autosomes should be kept. By default, sex chromosomes are kept, but
#'   mitochondrial and non-standard chromosomes are removed.
#' @param genome An optional string that, if supplied, will be used to set the
#'   genome of \code{gr}.
#'
#' @return A GRanges object in which both ranges and \code{seqinfo} associated
#'   with trimmed chromosomes have been removed.
#'
#' @details Standard chromosomes are defined using the
#'   \code{\link[GenomeInfoDb:seqlevels-wrappers]{standardChromosomes}} function
#'   from the \code{GenomeInfoDb} package.
#'
#' @author Mike DeBerardine
#' @seealso \code{\link[GenomeInfoDb:seqlevels-wrappers]{
#'   GenomeInfoDb::standardChromosomes}}
#'
#' @export
#' @importFrom GenomeInfoDb standardChromosomes seqlevels keepSeqlevels
#'   sortSeqlevels
#' @importFrom methods is
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
                            keep.nonstandard = FALSE, genome = NULL) {

    if (is.list(gr) || is(gr, "GRangesList"))
        return(lapply(gr, tidyChromosomes, keep.X, keep.Y, keep.M,
                      keep.nonstandard, genome))

    if (!is.null(genome))
        genome(gr) <- genome

    chrom <- standardChromosomes(gr)

    if (keep.nonstandard) chrom <- seqlevels(gr)
    if (!keep.X)  chrom <- chrom[ chrom != "chrX" ]
    if (!keep.Y)  chrom <- chrom[ chrom != "chrY" ]
    if (!keep.M)  chrom <- chrom[ (chrom != "chrM") & (chrom != "chrMT") ]

    if (is(gr, "Seqinfo")) {
        gr <- keepSeqlevels(gr, chrom)
    } else {
        gr <- keepSeqlevels(gr, chrom, pruning.mode = "tidy")
    }
    sortSeqlevels(gr)
}


#' Import basepair-resolution files
#'
#' Import functions for plus/minus pairs of \code{bigWig} or \code{bedGraph}
#' files.
#'
#' @param plus_file,minus_file Paths for strand-specific input files, or a
#'   vector of such paths. If vectors are given, the user should take care that
#'   the orders match!
#' @param genome Optional string for UCSC reference genome, e.g. "hg38". If
#'   given, non-standard chromosomes are trimmed, and options for sex and
#'   mitochondrial chromosomes are applied.
#' @param keep.X,keep.Y,keep.M,keep.nonstandard Logicals indicating which
#'   non-autosomes should be kept. By default, sex chromosomes are kept, but
#'   mitochondrial and non-standard chromosomes are removed.
#' @param ncores Number of cores to use, if importing multiple objects
#'   simultaneously.
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
#'   is to be maintained for each sequenced molecule. For this purpose, one use
#'   bedGraphs to store entire reads, with the \code{score} representing the
#'   number of reads sharing identical 5' and 3' ends. However,
#'   \code{import_bedGraph} doesn't modify the information in the bedGraph
#'   files. If the bedGraph file represents basepair-resolution coverage data,
#'   then users can coerce it to a basepair-resolution GRanges object by using
#'   \code{\link[BRGenomics:getStrandedCoverage]{getStrandedCoverage}} followed
#'   by \code{\link[BRGenomics:makeGRangesBRG]{makeGRangesBRG}}.
#'
#' @author Mike DeBerardine
#' @seealso \code{\link[BRGenomics:tidyChromosomes]{tidyChromosomes}},
#'   \code{\link[rtracklayer:export]{rtracklayer::import}}
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
#' # import bigWigs (not supported on windows)
#' if (.Platform$OS.type == "unix") {
#'     PROseq <- import_bigWig(p.bw, m.bw, genome = "dm6")
#'     PROseq
#' }
#'
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
#' @param makeBRG If \code{TRUE} (the default), the output ranges are made
#'   single-width using \code{\link[BRGenomics:makeGRangesBRG]{makeGRangesBRG}}
#' @export
#' @import rtracklayer
#' @importFrom GenomicRanges score score<- strand<-
#' @importFrom GenomeInfoDb genome<-
#' @importFrom parallel mcMap
import_bigWig <- function(plus_file = NULL, minus_file = NULL, genome = NULL,
                          keep.X = TRUE, keep.Y = TRUE, keep.M = FALSE,
                          keep.nonstandard = FALSE, makeBRG = TRUE,
                          ncores = getOption("mc.cores", 2L)) {

    if (length(plus_file) > 1 || length(minus_file) > 1) {
        if (!is.null(plus_file) && !is.null(minus_file))
            if (length(plus_file) != length(minus_file))
                stop("plus_file and minus_file are not the same length")
        if (is.null(plus_file))  plus_file <- list(NULL)
        if (is.null(minus_file))  minus_file <- list(NULL)
        if (is.null(genome))  genome <- list(NULL)
        return(mcMap(import_bigWig, plus_file, minus_file, genome, keep.X,
                     keep.Y, keep.M, keep.nonstandard, ncores = 1,
                     mc.cores = ncores))
    }

    p_gr <- m_gr <- GRanges() # initialize
    if (!is.null(plus_file)) {
        p_gr <- import.bw(plus_file)
        strand(p_gr) <- "+"
    }
    if (!is.null(minus_file)) {
        m_gr <- import.bw(minus_file)
        score(m_gr) <- abs(score(m_gr)) # make scores positive
        strand(m_gr) <- "-"
    }
    suppressWarnings( gr <- c(p_gr, m_gr) )

    # scores are imported as doubles by default; if whole numbers, make integers
    gr <- .try_int_score(gr)

    # Make the width of each range equal to 1
    if (makeBRG) gr <- makeGRangesBRG(gr)

    if (!is.null(genome)) {
        genome(gr) <- genome
        gr <- tidyChromosomes(gr, keep.X = keep.X, keep.Y = keep.Y,
                              keep.M = keep.M,
                              keep.nonstandard = keep.nonstandard)
    }
    sort(gr)
}


#' @rdname import-functions
#' @export
#' @import rtracklayer
#' @importFrom GenomicRanges score score<- strand<-
#' @importFrom GenomeInfoDb genome<-
#' @importFrom parallel mcMap
import_bedGraph <- function(plus_file = NULL, minus_file = NULL, genome = NULL,
                            keep.X = TRUE, keep.Y = TRUE, keep.M = FALSE,
                            keep.nonstandard = FALSE,
                            ncores = getOption("mc.cores", 2L)) {

    if (length(plus_file) > 1 || length(minus_file) > 1) {
        if (!is.null(plus_file) && !is.null(minus_file))
            if (length(plus_file) != length(minus_file))
                stop("plus_file and minus_file are not the same length")
        if (is.null(plus_file))  plus_file <- list(NULL)
        if (is.null(minus_file))  minus_file <- list(NULL)
        if (is.null(genome))  genome <- list(NULL)
        return(mcMap(import_bedGraph, plus_file, minus_file, genome, keep.X,
                     keep.Y, keep.M, keep.nonstandard, ncores = 1,
                     mc.cores = ncores))
    }

    p_gr <- m_gr <- GRanges() # initialize
    if (!is.null(plus_file)) {
        p_gr <- import.bedGraph(plus_file)
        strand(p_gr) <- "+"
    }
    if (!is.null(minus_file)) {
        m_gr <- import.bedGraph(minus_file)
        score(m_gr) <- abs(score(m_gr)) # make scores positive
        strand(m_gr) <- "-"
    }
    suppressWarnings( gr <- c(p_gr, m_gr) )

    # scores are imported as doubles by default; if whole numbers, make integers
    gr <- .try_int_score(gr)

    if (!is.null(genome)) {
        genome(gr) <- genome
        gr <- tidyChromosomes(gr, keep.X = keep.X, keep.Y = keep.Y,
                              keep.M = keep.M,
                              keep.nonstandard = keep.nonstandard)
    }
    sort(gr)
}


#' Import bam files
#'
#' Import single-end or paired-end bam files as GRanges objects, with various
#' processing options. It is highly recommend to index the BAM file first.
#'
#' @param file Path of a bam file, or a vector of paths.
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
#'   reads. When set to \code{NULL} (the default), the first 100k reads are
#'   checked.
#' @param yieldSize The number of bam file records to process simultaneously,
#'   e.g. the "chunk size". Setting a higher chunk size will use more memory,
#'   which can increase speed if there is enough memory available. If chunking
#'   is not necessary, set to \code{NA}.
#' @param ncores Number of cores to use for importing bam files. Currently,
#'   multicore is only implemented for simultaneously importing multiple bam
#'   files. For smaller datasets or machines with higher memory, this can
#'   increase performance, but can otherwise lead to substantial performance
#'   penalties.
#'
#' @references Hojoong Kwak, Nicholas J. Fuda, Leighton J. Core, John T. Lis
#'   (2013). Precise Maps of RNA Polymerase Reveal How Promoters Direct
#'   Initiation and Pausing. \emph{Science} 339(6122): 950–953.
#'   \url{https://doi.org/10.1126/science.1229386}
#'
#'   Jason D. Buenrostro, Paul G. Giresi, Lisa C. Zaba, Howard Y. Chang, William
#'   J. Greenleaf (2013). Transposition of native chromatin for fast and
#'   sensitive epigenomic profiling of open chromatin, dna-binding proteins and
#'   nucleosome position. \emph{Nature Methods} 10: 1213–1218.
#'   \url{https://doi.org/10.1038/nmeth.2688}
#'
#' @return A GRanges object.
#' @author Mike DeBerardine
#' @export
#' @importFrom GenomicRanges GRanges strand strand<- resize
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
                       paired_end = NULL, yieldSize = NA,
                       ncores = 1) {

    trim.to <- match.arg(trim.to, c("whole", "5p", "3p", "center"))

    if (length(file) > 1) {
        if (is.null(paired_end))  paired_end <- list(NULL)
        return(mclapply(file, import_bam, mapq, revcomp, shift, trim.to,
                        ignore.strand, field, paired_end, yieldSize,
                        ncores = 1, mc.cores = ncores))
    }

    gr <- .import_bam(file, paired_end, yieldSize, mapq)

    # Apply Options
    if (revcomp || !all(shift == 0))
        is_plus <- as.character(strand(gr)) == "+"
    if (revcomp) {
        strand(gr) <- "+"
        strand(gr)[is_plus] <- "-"
        is_plus <- !is_plus
    }
    if (!all(shift == 0))
        suppressWarnings(gr <- .shift_gr(gr, is_plus, shift))
    if (trim.to != "whole") {
        opt <- paste0("opt.", trim.to)
        opt.arg <- list(opt.5p = "start", opt.3p = "end", opt.center = "center")
        gr <- resize(gr, width = 1, fix = opt.arg[[opt]])
    }
    if (ignore.strand)
        strand(gr) <- "*"

    gr <- sort(gr)
    if (!is.null(field))
        gr <- .collapse_reads(gr, field)
    return(gr)
}

# # @importFrom GenomicAlignments readGAlignmentPairs readGAlignments

#' @import Rsamtools GenomicAlignments
#' @importFrom GenomicRanges GRanges
.import_bam <- function(file, paired_end, yield_size, mapq) {
    ## This function avoids any use of bpiterate(), as we've had problems;
    ## this also affects other functions like GenomicFiles::reduceByYield

    # Get parameters
    bf <- BamFile(file, yieldSize = yield_size)
    param <- ScanBamParam(mapqFilter = mapq)
    if (is.null(paired_end))
        paired_end <- .quick_check_paired(file)

    fxn <- if (paired_end) readGAlignmentPairs else readGAlignments
    yfun <- function(x) fxn(x, param = param)

    if (is.na(yield_size)) {
        # No chunking
        return(sort(GRanges(yfun(bf))))

    } else {
        # With chunking
        open(bf)
        on.exit(close(bf))
        finished <- function(x) length(x) == 0L || is.null(x)

        reads <- yfun(bf)
        if (finished(reads)) # early exit if no reads found
            return(GRanges())
        reads <- list(reads)

        repeat {
            reads.i <- yfun(bf)
            if (finished(reads.i))
                break
            reads <- append(reads, reads.i)
        }
        return(sort(GRanges(do.call(c, reads))))
    }
}

#' @import Rsamtools
.quick_check_paired <- function(file) {
    # check the first 100k reads
    bfq <- BamFile(file, yieldSize = 1e5)
    open(bfq)
    on.exit(close(bfq))
    param <- ScanBamParam(what = "flag")
    flag <- scanBam(bfq, param = param)[[1]]$flag
    any(bamFlagTest(flag, "isPaired"))
}

#' @importFrom GenomicRanges shift
.shift_gr <- function(gr, is_plus, shift) {
    if (length(shift) == 1) {
        shift(gr, ifelse(is_plus, shift, -shift))
    } else if (length(shift) == 2) {
        genebodies(gr, shift[1], shift[2], min.window = 0)
    } else {
        stop(.nicemsg("shift argument must be a single number, or a numeric
                      vector of length 2"))
    }
}

#' @importFrom GenomicRanges countOverlaps mcols<-
.collapse_reads <- function(gr, field) {
    gr.out <- unique(gr)
    mcols(gr.out)[field] <- countOverlaps(gr.out, gr, type = "equal")
    return(gr.out)
}



#' @rdname import_bam
#' @export
import_bam_PROseq <- function(file, mapq = 20, revcomp = TRUE, shift = -1L,
                              trim.to = "3p", ignore.strand = FALSE,
                              field = "score", paired_end = NULL,
                              yieldSize = NA, ncores = 1) {

    import_bam(file = file, mapq = mapq, revcomp = revcomp, shift = shift,
               trim.to = trim.to, ignore.strand = ignore.strand, field = field,
               paired_end = paired_end, yieldSize = yieldSize, ncores = ncores)
}

#' @rdname import_bam
#' @export
import_bam_PROcap <- function(file, mapq = 20, revcomp = FALSE, shift = 0L,
                              trim.to = "5p", ignore.strand = FALSE,
                              field = "score", paired_end = NULL,
                              yieldSize = NA, ncores = 1) {

    import_bam(file = file, mapq = mapq, revcomp = revcomp, shift = shift,
               trim.to = trim.to, ignore.strand = ignore.strand, field = field,
               paired_end = paired_end, yieldSize = yieldSize, ncores = ncores)
}

# must also mention that people should not use these options if their ATAC-seq
# pipeline has already trimmed the read overhangs down

#' @rdname import_bam
#' @export
import_bam_ATACseq <- function(file, mapq = 20, revcomp = FALSE,
                               shift = c(4, -5), trim.to = "whole",
                               ignore.strand = TRUE, field = "score",
                               paired_end = TRUE, yieldSize = NA,
                               ncores = 1) {

    import_bam(file = file, mapq = mapq, revcomp = revcomp, shift = shift,
               trim.to = trim.to, ignore.strand = ignore.strand, field = field,
               paired_end = paired_end, yieldSize = yieldSize, ncores = ncores)
}
