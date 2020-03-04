(Changelog from 0.99.x temporary; will be collapsed)

## BRGenomics 0.99.12-13

* Revert rtracklayer::import link to point to function (rtracklayer:io not working...)

## BRGenomics 0.99.11

* Change tidyChromosomes test
* Remove indirect links in doc pages (use only exact names of man pages)

## BRGenomics 0.99.10

* (Fixed missed bigWig import in makeGRangesBRG example)

## BRGenomics 0.99.9

* Avoid importing bigWigs in test/examples if on Windows

## BRGenomics 0.99.8

* Made all examples and tests use single core (for Windows test build)

## BRGenomics 0.99.7

* Redefine mcMap to avoid errors on Windows, and allow Windows again

## BRGenomics 0.99.6

* Formally remove support for Windows

## BRGenomics 0.99.4

* For Bioconductor, made new branch `R3` to allow users to install under R version >=3.5
    + Branch `master` now requires R version 4.0 (development version)
* Updated readme to have users install from branch `R3`

## BRGenomics 0.99.1

* Remove .Rproj file from git repository

## BRGenomics 0.99.0 (Bioconductor submission)

* Update version for Bioconductor submission
* Various minor formatting changes to codebase

## BRGenomics 0.8.1

* Fixed bug in `import_bam()` that produced warnings when `shift` argument used to shift both 5' and 3' ends of reads (i.e. when `length(shift) == 2`)
* Updated included external datasets to be much smaller
    + Included bam file now <200 reads
    + Included bigWig and bedGraph files derived from the bam file
* Minor update to vignette to reflect the change in the bam file
* Updated the included `data()` objects (`PROseq` and `PROseq_paired`)
    + Shifted PRO-seq 3' bases to remove the run-on base, and updated associated package tests
    + Added xz compression to the files
* Streamlined some method dispatch in `genebodies()` and `getDESeqDataSet()` functions

## BRGenomics 0.8.0

* Add support for lists in data import functions
* Add the convenience function `applyNFsGRanges()`
* Significant internal changes to `import_bam()`
    + New test for paired-end reads (deprecated use of `Rsamtools::testPairedEndBam()`)
    + Avoids any internal use of `bpiterate()`
    + Dropped dependency on `GenomicFiles` package

## BRGenomics 0.7.10

* Add `intersectByGene()` and `reduceByGene()` functions
* Minor vignette updates

## BRGenomics 0.7.8

* Substantially updated vignette
* Fully load `rtracklayer` (so completely exported to users)
* Add `isBRG()` function
* Fixed bug in `spikeInNormGRanges()` that failed to remove spike-in reads (aside from maintaining those reads, normalization was otherwise correct)
* Fixed minor bug in which `metaSubsample()` automatically added rownames with list input

## BRGenomics 0.7.7

* Bug fix in `aggregateByNdimensionalBins()` affecting simultaneous aggregation of multiple data
* Minor updates to documentation, including an error in `getDESeqResults()`
* Slightly expanded vignette

## BRGenomics 0.7.5

* Update bootstrapping functions
    + Add blacklisting support for `metaSubsample()`
    + Related to blacklisting, NA values now ignored in bootstrapping
* Add additional `melt` options for signal counting functions
* Further expanded support for list inputs (lists of GRanges datasets), including in `getStrandedCoverage()`
* Add explicit support for blacklisting in `getDESeqDataSet()`
* Rewrite n-dimensional binning functions, and add function for aggregating data within n-dimensional bins
    + Changed arguments in `binNdimensions()` to only accept dataframe inputs
    + Add `densityInNdimensionalBins()` function to count points in each bin
    + Add `aggregateByNdimensionalBins()` function to aggregate data within bins using arbitrary functions
* Added arguments for setting sample names in spike-in/normalization functions
* Various improvements and streamlining for method dispatch and flexibility

## BRGenomics 0.7.0

* Added functions for counting and filtering spike-in reads
* Added functions for generating spike-in normalization factors
* Added support for lists of GRanges datasets throughout, including all signal counting functions
* Updating signal counting functions with a `blacklist` argument, for ignoring reads from user-supplied regions
* Added wrappers for `import_bam()` for several common use cases
* Update `getCountsByPositions()`: 
    + Added a `melt` option for returning melted dataframes
    + Now returns an error by default if multi-width regions are given (must be explicit) 
* Changed argument order in `getMaxPositionsBySignal()`

## BRGenomics 0.5.6

* Update `import_bam()` function
    + Added support for paired-end bam files
    + Added the `shift` argument
* Made `metaSubsample()` functions robust to unevaluated inputs 
* Small performance improvement to `genebodies()` function
* Multicore usage is again the default for `getDESeqResults()`

## BRGenomics 0.5.3

* Substantial performance improvement for `mergeGRangesData()`
* Make single-core the default for `getDESeqResults()`
* Fixed errant warning message in `binNdimensions()` for integer inputs
* Update namespace to fully import GenomicRanges and S4Vectors
* Changed R dependency, evidently required by the updated data objects
* Various documentation updates
* Added package documentation page
* Minor changes to examples

## BRGenomics 0.5.0

* Increased support for multiplexed GRanges objects across all functions
* Increased performance for bootstrapping over multiple fields/multiplexed GRanges objects
* Removed requirement to set 'field' argument for gettings counts over multiplexed GRanges

## BRGenomics 0.4.7

* Rewrote `mergeGRangesData()` to support the creation of multiplexed GRanges objects
* Made `getCountsByRegions()` and `getCountsByPositions()` to return integers if input signal is integer

## BRGenomics 0.4.5

* Added and modified numerous arguments
* Increased support for normalization factors across signal counting functions
* Modified behavior of bootstrapping over multiple fields by removing explicit access to .Random.seed

## BRGenomics 0.4.1

* Various documentation updates

## BRGenomics 0.4.0

* Added a `NEWS.md` file to track changes
* Added support for bam files
