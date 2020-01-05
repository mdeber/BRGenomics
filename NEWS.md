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
