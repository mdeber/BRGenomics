## BRGenomics 0.7.7

* Bug fix in `aggregateByNdimensionalBins` affecting simultaneous aggregation of multiple data
* Minor updates to documentation, including an error in `getDESeqResults`
* Slightly expanded vignette

## BRGenomics 0.7.5

* Update bootstrapping functions
    + Add blacklisting support for `metaSubsample`
    + Related to blacklisting, NA values now ignored in bootstrapping
* Add additional `melt` options for signal counting functions
* Further expanded support for list inputs (lists of GRanges datasets), including in `getStrandedCoverage`
* Add explicit support for blacklisting in `getDESeqDataSet`
* Rewrite n-dimensional binning functions, and add function for aggregating data within n-dimensional bins
    + Changed arguments in `binNdimensions` to only accept dataframe inputs
    + Add `densityInNdimensionalBins` function to count points in each bin
    + Add `aggregateByNdimensionalBins` function to aggregate data within bins using arbitrary functions
* Added arguments for setting sample names in spike-in/normalization functions
* Various improvements and streamlining for method dispatch and flexibility

## BRGenomics 0.7.0

* Added functions for counting and filtering spike-in reads
* Added functions for generating spike-in normalization factors
* Added support for lists of GRanges datasets throughout, including all signal counting functions
* Updating signal counting functions with a `blacklist` argument, for ignoring reads from user-supplied regions
* Added wrappers for `import_bam` for several common use cases
* Update `getCountsByPositions`: 
    + Added a `melt` option for returning melted dataframes
    + Now returns an error by default if multi-width regions are given (must be explicit) 
* Changed argument order in `getMaxPositionsBySignal`

## BRGenomics 0.5.6

* Update `import_bam` function
    + Added support for paired-end bam files
    + Added the `shift` argument
* Made `metaSubsample` functions robust to unevaluated inputs 
* Small performance improvement to `genebodies` function
* Multicore usage is again the default for `getDESeqResults`

## BRGenomics 0.5.3

* Substantial performance improvement for `mergeGRangesData`
* Make single-core the default for `getDESeqResults`
* Fixed errant warning message in `binNdimensions` for integer inputs
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

* Rewrote `mergeGRangesData` to support the creation of multiplexed GRanges objects
* Made `getCountsByRegions` and `getCountsByPositions` to return integers if input signal is integer

## BRGenomics 0.4.5

* Added and modified numerous arguments
* Increased support for normalization factors across signal counting functions
* Modified behavior of bootstrapping over multiple fields by removing explicit access to .Random.seed

## BRGenomics 0.4.1

* Various documentation updates

## BRGenomics 0.4.0

* Added a `NEWS.md` file to track changes
* Added support for bam files
