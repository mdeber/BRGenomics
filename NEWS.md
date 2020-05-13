## BRGenomics 1.0.1 (release) & 1.1.1 (devel)

* Merge changes from 0.99.33 into release branch
* Bug workaround for `getStrandedCoverage()` and dependent methods for getting
coverage of normalized data (apparent bug in `IRanges::coverage()` for 
weighting by normalized values)
* Bug fix in `getCountsByPositions()` for getting counts over an unstranded
region with `expand_ranges=TRUE`

## BRGenomics 0.99.33

* Add support for log fold-change shrinkage in `getDESeqResults()`
* Minor bugfix generating a warning under some conditions when non-integer 
binsizes are used 

## BRGenomics 0.99.30

* Change all `ncores` options to default to `getOption("mc.cores", 2L)`
* Rename n-dimensional binning functions `aggregateByNdimBins()` and 
`densityInNdimBins()`
* Redone documentation website, splitting up the vignette
* Updated package vignettes/guide pages in several places
* Other minor internal changes

## BRGenomics 0.99.25

* Bug fixes with `expand_ranges` arguments affecting `getCountsByRegions()`, 
`subsampleGRanges()`, and `getSpikeInNFs()`
* Expanded testing, particularly for `expand_ranges` arguments and 
`import_bam()`
* Added options in `mergeReplicates()`
* Various small doc updates and minor internal changes

## BRGenomics 0.99.23

* Added support for non-basepair-resolution GRanges throughout, via the 
`expand_ranges` argument
    + Substantial performance benefits for less-sparse datasets (e.g. whole 
    read coverage)
    + Supported everywhere, including counting functions, subsampling, merging, 
    normalization, etc.
* Rewrite of `mergeGRangesData()`: 
    + Substantial performance improvements for most datasets
    + No longer requires basepair-resolution GRanges objects
    + Added options and flexibility for merging reads as well as coverage data
* Add a `mergeReplicates()` function
* Rewrite of `makeGRangesBRG()` that significantly improves performance for 
sparser datasets (the datasets for which using the function makes the most 
sense)
* `subsampleGRanges()` no longer returns a basepair-resolution GRanges by 
default
* When `field=NULL`, `applyNFsGRanges()` no longer returns a 
basepair-resolution GRanges by default
* Add `use_bin_numbers` option to n-dimension binning functions; setting to 
false allows returning of bin values (the bin center) instead of the ordinal 
bin numbers (indexes)
* Quietly adding support for `GRangesList` objects (not fully tested)

## BRGenomics 0.99.15

* Add pre-filtering to counting functions for performance
* Some additional clarification of readcounts vs. coverage signal in counting 
and import functions
* Change tidyChromosomes test
* Remove indirect links in doc pages (use only exact names of man pages)

## BRGenomics 0.99.10

* Code modifications to pass Bioconductor test builds on Windows:
    + Make all examples and tests single core
    + Internally (not exported) redefine mcMap (current implementation in 
    package parallel needs to be modified)
    + In tests and examples, test if on Windows before attempting any bigWig 
    file import

## BRGenomics 0.99.0

* Changes for Bioconductor submission:
    + Move to package versioning 0.99.x
    + Update R requirement to version 4.0
    + Add new branch `R3` to allow users to install under R version >=3.5
* Various minor formatting changes to codebase

## BRGenomics 0.8.1

* Fixed bug in `import_bam()` that produced warnings when `shift` argument used 
to shift both 5' and 3' ends of reads (i.e. when `length(shift) == 2`)
* Updated included external datasets to be much smaller
    + Included bam file now <200 reads
    + Included bigWig and bedGraph files derived from the bam file
* Minor update to vignette to reflect the change in the bam file
* Updated the included `data()` objects (`PROseq` and `PROseq_paired`)
    + Shifted PRO-seq 3' bases to remove the run-on base, and updated 
    associated package tests
    + Added xz compression to the files
* Streamlined some method dispatch in `genebodies()` and `getDESeqDataSet()` 
functions

## BRGenomics 0.8.0

* Add support for lists in data import functions
* Add the convenience function `applyNFsGRanges()`
* Significant internal changes to `import_bam()`
    + New test for paired-end reads (deprecated use of 
    `Rsamtools::testPairedEndBam()`)
    + Avoids any internal use of `bpiterate()`
    + Dropped dependency on `GenomicFiles` package

## BRGenomics 0.7.10

* Add `intersectByGene()` and `reduceByGene()` functions
* Minor vignette updates

## BRGenomics 0.7.8

* Substantially updated vignette
* Fully load `rtracklayer` (so completely exported to users)
* Add `isBRG()` function
* Fixed bug in `spikeInNormGRanges()` that failed to remove spike-in reads 
(aside from maintaining those reads, normalization was otherwise correct)
* Fixed minor bug in which `metaSubsample()` automatically added rownames with 
list input

## BRGenomics 0.7.7

* Bug fix in `aggregateByNdimensionalBins()` affecting simultaneous aggregation 
of multiple data
* Minor updates to documentation, including an error in `getDESeqResults()`
* Slightly expanded vignette

## BRGenomics 0.7.5

* Update bootstrapping functions
    + Add blacklisting support for `metaSubsample()`
    + Related to blacklisting, NA values now ignored in bootstrapping
* Add additional `melt` options for signal counting functions
* Further expanded support for list inputs (lists of GRanges datasets), 
including in `getStrandedCoverage()`
* Add explicit support for blacklisting in `getDESeqDataSet()`
* Rewrite n-dimensional binning functions, and add function for aggregating 
data within n-dimensional bins
    + Changed arguments in `binNdimensions()` to only accept dataframe inputs
    + Add `densityInNdimensionalBins()` function to count points in each bin
    + Add `aggregateByNdimensionalBins()` function to aggregate data within 
    bins using arbitrary functions
* Added arguments for setting sample names in spike-in/normalization functions
* Various improvements and streamlining for method dispatch and flexibility

## BRGenomics 0.7.0

* Added functions for counting and filtering spike-in reads
* Added functions for generating spike-in normalization factors
* Added support for lists of GRanges datasets throughout, including all signal 
counting functions
* Updating signal counting functions with a `blacklist` argument, for ignoring 
reads from user-supplied regions
* Added wrappers for `import_bam()` for several common use cases
* Update `getCountsByPositions()`: 
    + Added a `melt` option for returning melted dataframes
    + Now returns an error by default if multi-width regions are given (must be 
    explicit) 
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
* Increased performance for bootstrapping over multiple fields/multiplexed 
GRanges objects
* Removed requirement to set 'field' argument for gettings counts over 
multiplexed GRanges

## BRGenomics 0.4.7

* Rewrote `mergeGRangesData()` to support the creation of multiplexed GRanges 
objects
* Made `getCountsByRegions()` and `getCountsByPositions()` to return integers 
if input signal is integer

## BRGenomics 0.4.5

* Added and modified numerous arguments
* Increased support for normalization factors across signal counting functions
* Modified behavior of bootstrapping over multiple fields by removing explicit 
access to .Random.seed

## BRGenomics 0.4.1

* Various documentation updates

## BRGenomics 0.4.0

* Added a `NEWS.md` file to track changes
* Added support for bam files
