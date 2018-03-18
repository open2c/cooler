### 0.7.8 (2018-03-18) ###

* New `cooler cload pairs` command provides index-free loading of pairs.
* Changed name of `create_from_unsorted` to more correct `create_from_unordered`.

Bug fixes

* Fixed broken use of single-file temporary store in `create_from_unordered`.
* Added heuristic in pairix cload to prevent excessively large chunks. #92
* Added extra checks in `cload pairix` and `cload tabix`. #62, #75


### 0.7.7 (2018-03-16) ###

New features

* Implementation of unsorted (index-free) loading
    * `cooler.io.create_from_unsorted` takes an iterable of pixel dataframe chunks that need not be properly sorted.
    * Use input sanitization procedures for pairs `sanitize_records` and binned data `sanitize_pixels` to feed data to `create_from_unsorted`. #87 #108 #109
    * The `cooler load` command is now index-free: unsorted `COO` and `BG2` input data can be streamed in. #90. This will soon be implemented as an option for loading pairs as well.
* Prevent `cooler balance` command from exiting with non-zero status upon failed convergence using convergence error policies. #93
* Improve the `create` API to support pandas read_csv-style `columns` and `dtype` kwargs to add extra value columns or override default dtypes. #108
* Experimental implementation of trans-only balancing. #56

Bug fixes

* Fix argmax deprecation. #99


### 0.7.6 (2017-10-31) ###
New features

* Cooler zoomify with explicit resolutions
* Towards standardization of multicooler structure
* Support for loading 1-based COO triplet input files

Bug fixes

* Fixed issue of exceeding header limit with too many scaffolds. If header size is exceeded, chrom IDs are stored as raw integers instead of HDF5 enums. There should be no effect at the API level.
* Fixed issue of single-column chromosomes files not working in `cload`.
* Fixed edge case in performing joins when using both `as_pixels` and `join` options in the matrix selector.

Happy Halloween!

### 0.7.5 (2017-07-13) ###
* Fix pandas issue affecting cases when loading single chromosomes
* Add transform options to higlass API

### 0.7.4 (2017-05-25) ###
* Fix regression in automatic --balance option in cooler zoomify
* Fix special cases where cooler.io.create and append would not work with certain inputs

### 0.7.3 (2017-05-22) ###
* Added function to print higlass zoom resolutions for a given genome and base resolution.

### 0.7.2 (2017-05-09) ###

* Improve chunking and fix pickling issue with aggregating very large text datasets
* Restore zoom binsize metadata to higlass files

### 0.7.1 (2017-04-29) ###

* `cooler load` command can now accept supplemental pixel fields and custom field numbers
* Fix parsing errors with unused pixel fields
* Eliminate hard dependence on dask to make pip installs simpler. Conda package will retain dask as a run time requirement.

### 0.7.0 (2017-04-27) ###

New features

* New Cooler URIs: Full support for Cooler objects anywhere in the data hierarchy of a .cool file
* Experimental dask support via `cooler.contrib.dask`
* New explicit bin blacklist option for `cooler balance`
* Various new CLI tools:
    * `cooler list`
    * `cooler copy`
    * `cooler merge`
* `cooler csort` now produces Pairix files by default
* `cooler load` now accepts two types of matrix text input formats
    * 3-column sparse matrix
    * 7-column bg2.gz (2D bedGraph) indexed with Pairix (e.g. using csort)
* `cooler coarsegrain` renamed `cooler coarsen`
* Multi-resolution HiGlass input files can now be generated with the `cooler zoomify` command
* More flexible API functions to create and append columns to Coolers in `cooler.io`

Backwards-incompatible changes

* `cooler.io.create` signature changed; `chromsizes` argument is deprecated.
* `cooler csort` argument order changed

Bug fixes

* Chromosome name length restriction removed
* `Cooler.open` function now correctly opens the specific root group of the Cooler and behaves like a proper context manager in all cases

### 0.6.6 (2017-03-21) ###

* Chromosome names longer than 32 chars are forbidden for now
* Improved pairix and tabix iterators, dropped need for slow first pass over contacts

### 0.6.5 (2017-03-18) ###

* Fixed pairix aggregator to properly deal with autoflipping of pairs

### 0.6.4 (2017-03-17) ###

* Migrated higlass multires aggregator to `cooler coarsegrain` command
* Fixed pairix aggregator to properly deal with autoflipping of pairs

### 0.6.3 (2017-02-22) ###

* Merge PairixAggregator patch from Soo.
* Update repr string
* Return matrix scale factor in balance stats rather than the bias scale factor: #35.

### 0.6.2 (2017-02-12) ###

Fixed regressions in

* cooler cload tabix/pairix failed on non-fixed sized bins
* cooler show

### 0.6.1 (2017-02-06) ###

* This fixes stale build used in bdist_wheel packaging that broke 0.6.0. #29

### 0.6.0 (2017-02-03) ###

General

* Dropped Python 3.3 support. Added 3.6 support.
* Added `contrib` subpackage containing utilities for higlass, including multires aggregation.
* Fixed various issues with synchronizing read/write multiprocessing with HDF5.
* Replacing prints with logging.
* Added sandboxed `tools` module to develop utilities for out-of-core algorithms using Coolers.

New features

* Cooler objects have additional convenience properties `chromsizes`, `chromnames`.
* New file introspection functions `ls` and `is_cooler` to support nested Cooler groups.
* Cooler initializer can accept a file path and path to Cooler group.
* `cload` accepts contact lists in hiclib-style HDF5 format, the legacy tabix-indexed format, and new pairix-indexed format.

Backwards-incompatible changes

* `create` only accepts a file path and optional group path instead of an open file object.
* `Cooler.matrix` selector now returns a balanced dense 2D NumPy array by default. Explicitly set `balance` to False to get raw counts and set `sparse` to True to get a `coo_matrix` as per old behavior.
* Command line parameters of `cload` changed significantly

Bug fixes

* Fixed bug in `csort` that led to incorrect triangularity of trans read pairs.


### 0.5.3 (2016-09-10) ###

* Check for existence of required external tools in CLI
* Fixed `cooler show` incompatibility with older versions of matplotlib
* Fixed `cooler.annotate` to work on empty dataframe input
* Fixed broken pipe signals not getting suppressed on Python 2
* `cooler cload` raises a warning when bin file lists a contig missing from the contact list


### 0.5.2 (2016-08-26) ###

* Fix bug in `cooler csort` parsing of chromsizes file.
* Workaround for two locale-related issues on Python 3. Only affects cases where a machine's locale is set to ASCII or Unices which use the ambiguous C or POSIX locales.
* Fix typo in setup.py and add pysam to dependencies.


### 0.5.1 (2016-08-24) ###

* Bug fix in input parser to `cooler csort`
* Update triu reording awk template in `cooler csort`
* Rename `cooler binnify` to `cooler makebins`. Binnify sounds like "aggregate" which is what `cload` does.


### 0.5.0 (2016-08-24) ###

* Most scripts ported over to a new command line interface using the Click framework with many updates.
* New `show` and `info` scripts.
* Updated Readme.
* Minor bug fixes.


### 0.4.0 (2016-08-18) ###

Schema

* Updated file schema: v2
* `/bins/chroms` is now an enum instead of string column

API changes

* Table views are a bit more intuitive: selecting field names on table view objects returns a new view on the subset of columns.
* New API function: `cooler.annotate` for doing joins

Features

* Support for nested Cooler "trees" at any depth in an HDF5 hierarchy
* Refactored `cooler.io` to provide "contact readers" that process different kinds of input (aggregate from a contact list, load from an existing matrix, etc.)
* Added new scripts for contact aggregation, loading, dumping and balancing


### 0.3.0 (2016-02-18) ###

* 2D range selector `matrix()` now provides either rectangular data as coo_matrix or triangular data as a pixel table dataframe.
* Added binning support for any genome segmentation (i.e., fixed or variable bin width).
* Fixed issues with binning data from mapped read files.
* Genomic locus string parser now accepts ENSEMBL-style number-only chromosome names and FASTA-style sequence names containing pipes.


### 0.2.1 (2016-02-07) ###

* Fixed bintable region fetcher


### 0.2 (2016-01-17) ###

* First beta release


### 0.1 (2015-11-22) ###

* Working initial prototype.
