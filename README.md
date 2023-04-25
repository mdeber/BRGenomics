
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BRGenomics

<!-- badges: start (versions currently manual...) -->

[![platforms](https://img.shields.io/badge/platforms-linux%20%7C%20osx%20%7C%20win-yellow.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/BRGenomics/)
[![release](https://img.shields.io/badge/release%20version-1.12.0-navy.svg)](https://www.bioconductor.org/packages/BRGenomics)
[![updated](http://www.bioconductor.org/shields/lastcommit/release/bioc/BRGenomics.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/BRGenomics/)
[![devel](https://img.shields.io/badge/devel%20version-1.13.0-orange.svg)](https://github.com/mdeber/BRGenomics)
[![last
commit](https://img.shields.io/github/last-commit/mdeber/BRGenomics.svg)](https://github.com/mdeber/BRGenomics/commits/master)
[![BioC
status](http://www.bioconductor.org/shields/build/release/bioc/BRGenomics.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/BRGenomics/)
[![test
coverage](https://codecov.io/gh/mdeber/BRGenomics/branch/master/graph/badge.svg)](https://codecov.io/gh/mdeber/BRGenomics)

<!-- badges: end -->

Efficient tools for the analysis of high-resolution genomics data in R.

Explore the documentation: <https://mdeber.github.io/>

# Installation

BRGenomics can be installed directly from Bioconductor:

``` r
if (!"BiocManager" %in% installed.packages())
    install.packages("BiocManager")
BiocManager::install("BRGenomics")
```

Alternatively, the latest development version can be installed from
[GitHub](https://github.com/mdeber/BRGenomics):

``` r
if (!"remotes" %in% installed.packages())
    install.packages("remotes")
remotes::install_github("mdeber/BRGenomics")
```

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

-   No support for parallel/multicore processing
-   No support for import bigWig files

------------------------------------------------------------------------

1.  Hojoong Kwak, Nicholas J. Fuda, Leighton J. Core, John T. Lis
    (2013). Precise Maps of RNA Polymerase Reveal How Promoters Direct
    Initiation and Pausing. *Science* **339**(6122): 950–953.
    <https://doi.org/10.1126/science.1229386>
