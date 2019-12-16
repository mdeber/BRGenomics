#' Drosophila S2 cell PRO-seq, chr2L
#'
#' Original data obtained from:
#'
#' Hojoong Kwak, Nicholas J. Fuda, Leighton J. Core, John T. Lis (2013).
#' Precise Maps of RNA Polymerase Reveal How Promoters Direct Initiation and
#' Pausing. Science, 339(6122), 950–953. https://doi.org/10.1126/science.1229386
#'
#' PRO-seq was originally sequenced on a GAIIx by the authors. Reads were
#' trimmed using Cutadapt; re-aligned using bowtie2; Rsamtools and
#' GenomeAlignments were used to truncate alignments to only the 3' ends of the
#' nascent RNA; base-pair resolution coverage was calculated using IRanges;
#' bigWig files were exported using rtracklayer; and this package was used
#' to import the bigWigs and re-export bigWigs containing only signal for
#' chromosome 2L.
