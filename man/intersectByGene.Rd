% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/roi_functions.R
\name{intersectByGene}
\alias{intersectByGene}
\alias{reduceByGene}
\title{Intersect or reduce ranges according to gene names}
\usage{
intersectByGene(regions.gr, gene_names)

reduceByGene(regions.gr, gene_names, disjoin = FALSE)
}
\arguments{
\item{regions.gr}{A GRanges object containing regions of interest. If
\code{regions.gr} has the class \code{list}, \code{GRangesList}, or
\code{CompressedGRangesList}, it will be treated as if each list element is
a gene, and the GRanges within are the ranges associated with that gene.}

\item{gene_names}{A character vector with the same length as
\code{regions.gr}.}

\item{disjoin}{Logical. If \code{disjoin = TRUE}, the output GRanges is
disjoint, and each output range will match a single gene name. If
\code{FALSE}, segments from different genes can overlap.}
}
\value{
A GRanges object whose individual ranges are named for the associated
  gene.
}
\description{
These functions divide up regions of interest according to associated names,
and perform an inter-range operation on them. \code{intersectByGene} returns
the "consensus" segment that is common to all input ranges, and returns no
more than one range per gene. \code{reduceByGene} collapses the input ranges
into one or more non-overlapping ranges that encompass all segments from the
input ranges.
}
\details{
These functions modify regions of interest that have associated
  names, such that several ranges share the same name, e.g. transcripts with
  associated gene names. Both functions "combine" the ranges on a
  gene-by-gene basis.

  \strong{intersectByGene}

  \emph{For each unique gene}, the segment that overlaps \emph{all} input
  ranges is returned. If no single range can be constructed that overlaps all
  input ranges, no range is returned for that gene (i.e. the gene is
  effectively filtered).

  In other words, for all the ranges associated with a gene, the
  most-downstream start site is selected, and the most upstream end site is
  selected.

  \strong{reduceByGene}

  \emph{For each unique gene}, the associated ranges are
  \code{\link[IRanges:inter-range-methods]{reduced}} to produce one or
  more non-overlapping ranges. The output range(s) are effectively a
  \code{union} of the input ranges, and cover every input base.

  With \code{disjoin = FALSE}, no genomic segment is overlapped by more than
  one range \emph{of the same gene}, but ranges from different genes can
  overlap. With \code{disjoin = TRUE}, the output ranges are disjoint, and no
  genomic position is overlapped more than once. Any segment that overlaps
  more than one gene is removed, but any segment (i.e. any section of an
  input range) that overlaps only one gene is still maintained.
}
\section{Typical Uses}{


  A typical use for \code{intersectByGene} is to avoid transcript isoform
  selection, as the returned range is found in every isoform.

  \code{reduceByGene} can be used to count any and all reads that overlap any
  part of a gene's annotation, but without double-counting any of them. With
  \code{disjoin = FALSE}, no reads will be double-counted for the same gene,
  but the same read can be counted for multiple genes. With \code{disjoin =
  TRUE}, no read can be double-counted.
}

\examples{
# Make example data:
#  Ranges 1 and 2 overlap,
#  Ranges 3 and 4 are adjacent
gr <- GRanges(seqnames = "chr1",
              ranges = IRanges(start = c(1, 3, 7, 10),
                               end = c(4, 5, 9, 11)))
gr

#--------------------------------------------------#
# intersectByGene
#--------------------------------------------------#

intersectByGene(gr, c("A", "A", "B", "B"))

intersectByGene(gr, c("A", "A", "B", "C"))

gr2 <- gr
end(gr2)[1] <- 10
gr2

intersectByGene(gr2, c("A", "A", "B", "C"))

intersectByGene(gr2, c("A", "A", "A", "C"))

#--------------------------------------------------#
# reduceByGene
#--------------------------------------------------#

# For a given gene, overlapping/adjacent ranges are combined;
#  gaps result in multiple ranges for that gene
gr

reduceByGene(gr, c("A", "A", "A", "A"))

# With disjoin = FALSE, ranges from different genes can overlap
gnames <- c("A", "B", "B", "B")
reduceByGene(gr, gnames)

# With disjoin = TRUE, segments overlapping >1 gene are removed as well
reduceByGene(gr, gnames, disjoin = TRUE)

# Will use one more example to demonstrate how all
#  unambiguous segments are identified and returned
gr2

gnames
reduceByGene(gr2, gnames, disjoin = TRUE)

#--------------------------------------------------#
# reduceByGene, then aggregate counts by gene
#--------------------------------------------------#

# Consider if you did getCountsByRegions on the last output,
#  you can aggregate those counts according to the genes
gr2_redux <- reduceByGene(gr2, gnames, disjoin = TRUE)
counts <- c(5, 2, 3) # if these were the counts-by-regions
aggregate(counts ~ names(gr2_redux), FUN = sum)

# even more convenient if using a melted dataframe
df <- data.frame(gene = names(gr2_redux),
                 reads = counts)
aggregate(reads ~ gene, df, FUN = sum)

# can be extended to multiple samples
df <- rbind(df, df)
df$sample <- rep(c("s1", "s2"), each = 3)
df$reads[4:6] <- c(3, 1, 2)
df

aggregate(reads ~ sample*gene, df, FUN = sum)
}
\author{
Mike DeBerardine
}
