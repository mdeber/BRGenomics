
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BRGenomics

<!-- badges: start -->

<!-- badges: end -->

Efficient tools for the analysis of high-resolution genomics data in R.

Explore the documentation: <https://mdeber.github.io/>

# Installation

As of Bioconductor 3.11 (release date April 28, 2020), BRGenomics can be
installed directly from Bioconductor:

``` r
# install.packages("BiocManager")
BiocManager::install("BRGenomics")
```

Alternatively, the latest development version can be installed from  
[GitHub](https://github.com/mdeber/BRGenomics):

``` r
# install.packages("remotes")
remotes::install_github("mdeber/BRGenomics@R3")
```

*BRGenomics (and Bioconductor 3.11) require R version 4.0 (release April
24, 2020). Installing the `R3` branch, as shown above, is required to
install under R \>= 3.5*

If you install the development version from Github and you’re using
Windows, [Rtools for
Windows](https://cran.rstudio.com/bin/windows/Rtools/) is required.

## Features

See the [documentation website](https://mdeber.github.io/), which
includes an introductory vignette, as well as the documentation for
currently implemented functions, complete with demonstrative example
code. The package currently includes example PRO-seq
data<sup>\[1\]</sup>.

## Limitations for Windows users

  - No support for parallel/multicore processing
  - No support for import bigWig files

## To Do

  - Convert method dispatch to S4 generics
  - Formalize support for `GRangesList` objects (progress underway)
  - Write methods for `BigWigFile`/`BigWigFileList` objects (to avoid
    loading data into memory)
  - (Possibly) use `GPos` objects by default

-----

1.  Hojoong Kwak, Nicholas J. Fuda, Leighton J. Core, John T. Lis
    (2013). Precise Maps of RNA Polymerase Reveal How Promoters Direct
    Initiation and Pausing. *Science* **339**(6122): 950–953.
    <https://doi.org/10.1126/science.1229386>
