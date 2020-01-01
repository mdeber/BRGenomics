
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BRGenomics

<!-- badges: start -->

<!-- badges: end -->

Efficient tools for the analysis of high-resolution genomics data in R.

Explore the documentation: <https://mdeber.github.io/>

## Installation

Install development version from
[GitHub](https://github.com/mdeber/BRGenomics):

``` r
# install.packages("devtools")
devtools::install_github("mdeber/BRGenomics")
```

If you’re using Windows, [Rtools for
Windows](https://cran.rstudio.com/bin/windows/Rtools/) is required.

## Features

A full introductory vignette is in the works, but some quick highlights
include:

  - Simple data import for bigWig, bedGraph, and bam files
  - Read counting by genes and positions in genes (e.g. for heatmaps)
  - Bootstrapping mean signals at positions within genes (e.g. for
    metaplotting)
  - Straightforward integration with DESeq2, including:
      - Added support for discontinuous gene regions, so you can do
        things like filter out dREG peaks within a gene
      - Automatic use of reduced dispersion matrices, which means that
        DESeq2 won’t calculate dispersion from all samples within the
        dataset, but only those being compared – *important for
        perturbations with genome-wide consequences, whether or not
        there’s a global shift*
      - Easy, simultaneous generation of numerous sets of DESeq results
  - Random subsampling of reads from datasets
  - Convenient, easy-to-use modification of gene regions, e.g. to define
    promoters or genebodies
  - Consistent support for binning measurements throughout

For more, see the
[documentation](https://mdeber.github.io/reference/index.html), which
includes demonstrative example code. The package currently includes
example PRO-seq data<sup>\[1\]</sup>.

## Coming Soon

  - Calculation of spike-in normalization factors using several
    different methods
  - Functions to assess sufficiency of sequencing depth using meaningful
    metrics, e.g. pausing indices
  - Support for length-scaled, bootstrapped signal counting, i.e. signal
    is binned by percentage of gene length
  - Function to summarize/aggregate data in bins generated using
    `binNdimensions`
  - *Potentially:* Support of “Views” objects to avoid loading data into
    memory

## Known issues and limitations

  - Currently no support for multicore processing on Windows
  - Currently no support for multiple fields in `getStrandedCoverage`
  - Certain gene names can cause `getDESeqDataSet` to return an error
      - Systematic naming schemes work (e.g. ensembl IDs) while some
        lists of conventional gene names, i.e. “symbols”, will cause
        failure
  - In some contexts, using multicore causes errors in `getDESeqResults`
    If this occurs in your environment, set `ncores = 1`

-----

1.  Hojoong Kwak, Nicholas J. Fuda, Leighton J. Core, John T. Lis
    (2013). Precise Maps of RNA Polymerase Reveal How Promoters Direct
    Initiation and Pausing. *Science* **339**(6122): 950–953.
    <https://doi.org/10.1126/science.1229386>
