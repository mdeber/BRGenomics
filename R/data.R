#' PRO-seq data from Drosophila S2 cells
#'
#' PRO-seq data of Drosophila S2 cells, chromosome 4.
#'
#' Hojoong Kwak, Nicholas J. Fuda, Leighton J. Core, John T. Lis (2013).
#' Precise Maps of RNA Polymerase Reveal How Promoters Direct Initiation and
#' Pausing. Science, 339(6122), 950–953.
#' \url{https://doi.org/10.1126/science.1229386}
#'
#' @format A disjoint GRanges object with 47533 ranges with 1 metadata column:
#' \describe{
#'   \item{score}{coverage of PRO-seq read 3'-ends}
#'   ...
#' }
#' @source GEO Accession GSM1032758, run SRR611828.
"PROseq"


#' Paired PRO-seq data from Drosophila S2 cells
#'
#' PRO-seq data of Drosophila S2 cells, chromosome 4. Entire mapped reads kept.
#'
#' Hojoong Kwak, Nicholas J. Fuda, Leighton J. Core, John T. Lis (2013).
#' Precise Maps of RNA Polymerase Reveal How Promoters Direct Initiation and
#' Pausing. Science, 339(6122), 950–953.
#' \url{https://doi.org/10.1126/science.1229386}
#'
#' @format A GRanges object with 52464 ranges with 1 metadata column:
#' \describe{
#'   \item{score}{number of reads sharing the same mapped 5' and 3' ends}
#'   ...
#' }
#' @source GEO Accession GSM1032758, run SRR611828.
"PROseq_paired"
