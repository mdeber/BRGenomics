#' Drosophila S2 cell PRO-seq, chr4
#'
#' GEO Accession: GSM1032758 Run: SRR611828
#'
#' Hojoong Kwak, Nicholas J. Fuda, Leighton J. Core, John T. Lis (2013). Precise
#' Maps of RNA Polymerase Reveal How Promoters Direct Initiation and Pausing.
#' Science, 339(6122), 950–953. https://doi.org/10.1126/science.1229386
#'
#' Base-pair resolution PRO-seq was originally sequenced on a GAIIx by the
#' authors. Unprocessed reads were downloaded from GEO, adapters were trimmed
#' using Cutadapt, and trimmed reads were aligned using bowtie2 in
#' high-sensitivity end-to-end alignment mode.
#'
#' To generate bigWig and bedGraph files, Bam files were processed in R using
#' Rsamtools and GenomeAlignments, and reads were processed using GenomicRanges.
#' For all processing, reads were first reverse complemented by inverting the
#' strand of the GRanges object.
#'
#' bedGraph files with paired 5'- and 3'-end information were created in which
#' each range represents a unique pairing of 5' and 3' ends in the sequencing,
#' and the score metadata column represents the number of times that pairing was
#' observed. bigWig files containing only coverage of the 3' ends were created
#' by truncating reads to their 3' ends, and IRanges was used to calculate
#' coverage. For both the paired bedGraph files and 3' end bigWig files,
#' rtracklayer was used to export only the reads mapping to chromosome 4.

# infile <- BamFile(paste0(opt$input_dir, bam_file))
# param <- ScanBamParam(mapqFilter = opt$quality)
# alignment = readGAlignments(infile,
#                             use.names = F,
#                             param = param)
# reads <- GRanges(alignment)
#
# # take chromosome 4 only
# reads <- subset(reads, seqnames == "chr4")
#
# # PRO-seq reads must be reverse complemented
# strand(reads) = ifelse(strand(reads) == '+', '-', '+')
#
# # Get 3' ends for PRO-seq
# reads.3p = GenomicRanges::resize(reads, width = 1, fix = "end")
#
# # Collapse identical molecules (for paired and 3'-ends-only)
# # score for number of identical bits of information
# collapse_reads <- function(reads) {
#     reads <- sort(reads)
#     reads.unique <- unique(reads)
#     score(reads.unique) <- countOverlaps(reads.unique, reads, type = "equal")
#     return(reads.unique)
# }
# reads <- collapse_reads(reads)
# reads.3p <- collapse_reads(reads.3p)
#
# # separate by strand
# reads.p <- reads[which(strand(reads)=="+")]
# reads.m <- reads[which(strand(reads)=="-")]
#
# reads.3p.p <- reads[which(strand(reads.3p)=="+")]
# reads.3p.m <- reads[which(strand(reads.3p)=="-")]
#
# # make scores for minus strand reads negative
# score(reads.m) <- -(score(reads.m))
# score(reads.3p.m) <- -(score(reads.3p.m))
#
# rtracklayer::export.bedGraph(reads.p, "PROseq_dm6_chr4_plus.bedGraph")
# rtracklayer::export.bedGraph(reads.m, "PROseq_dm6_chr4_plus.bedGraph")
#
# rtracklayer::export.bw(reads.3p.p, "PROseq_dm6_chr4_plus.bw")
# rtracklayer::export.bw(reads.3p.m, "PROseq_dm6_chr4_plus.bw")






