#' PRO-seq data from Drosophila S2 cells
#'
#' PRO-seq data from chromosome 4 of Drosophila S2 cells.
#'
#' @references Hojoong Kwak, Nicholas J. Fuda, Leighton J. Core, John T. Lis
#'   (2013). Precise Maps of RNA Polymerase Reveal How Promoters Direct
#'   Initiation and Pausing. Science, 339(6122), 950–953.
#'   \url{https://doi.org/10.1126/science.1229386}
#'
#' @format A disjoint GRanges object with 47533 ranges with 1 metadata column:
#' \describe{
#'   \item{score}{coverage of PRO-seq read 3'-ends}
#'   ...
#' }
#' @source GEO Accession GSM1032758, run SRR611828.
#' @importFrom utils data
#'
#' @name PROseq-data

#' @rdname PROseq-data
#' @usage data(PROseq)
"PROseq"

#' @rdname PROseq-data
#' @usage data(PROseq_paired)
"PROseq_paired"


#' Ensembl transcripts for Drosophila melanogaster, dm6, chromosome 4.
#'
#' Transcripts obtained from annotation package
#' TxDb.Dmelanogaster.UCSC.dm6.ensGene, which was in turn made by the
#' Bioconductor Core Team from UCSC resources on 2019-04-25. Metadata columns
#' were obtained from "TXNAME" and "GENEID" columns. Data exported from
#' the TxDb package using GenomicFeatures version 1.35.11 on 2019-12-19.
#'
#' @format A GRanges object with 339 ranges and 2 metadata columns:
#' \describe{
#'   \item{tx_name}{Flybase unique identifiers for transcripts}
#'   \item{gene_id}{FLybase unique identifiers for the associated genes}
#' }
#' @source TxDb.Dmelanogaster.UCSC.dm6.ensGene version 3.4.6
#' @usage data(txs_dm6_chr4)
#' @importFrom utils data
"txs_dm6_chr4"
