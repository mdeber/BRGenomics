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
#' using fastx clipper, and trimmed reads were aligned using bowtie2 in
#' end-to-end alignment mode, with options set to "sensitive".
#'
#' After sorting, the bam file was modified to contain only the first 100 plus
#' strand reads (biological plus strand, i.e. the reverse of the bam flag, as
#' PRO-seq is sequenced as the reverse complement), and the last 100 minus
#' strand reads.
#'
#' To generate bigWig and bedGraph files, the bam file was processed in R using
#' BRGenomics (Rsamtools and GenomeAlignments). bedGraph files with paired 5'-
#' and 3'-end information were created in which each range represents a unique
#' pairing of 5' and 3' ends in the sequencing, and the score metadata column
#' represents the number of times that pairing was observed. bigWig files
#' containing only coverage of the 3' ends were created by truncating reads to
#' their 3' ends. For both the paired bedGraph files and 3' end bigWig files,
#' rtracklayer was used to export only the reads mapping to chromosome 4.

# using BRGenomics import_bam and rtracklayer export functions:
# library(BRGenomics)
# reads <- import_bam(path_to_bam_file, mapq = 20,
#                     revcomp = TRUE, shift = c(0, -1))
# export.bedGraph(subset(reads, strand == "+"), "PROseq_dm6_chr4_plus.bedGraph")
# export.bedGraph(subset(reads, strand == "-"), "PROseq_dm6_chr4_minus.bedGraph")
#
# ps <- import_bam_PROseq(path_to_bam_file, mapq = 20)
# export.bedGraph(subset(ps, strand == "+"), "PROseq_dm6_chr4_plus.bw")
# export.bedGraph(subset(ps, strand == "-"), "PROseq_dm6_chr4_minus.bw")
