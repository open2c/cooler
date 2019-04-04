# Release notes #

## v0.8.4 ##

Date: 2019-04-04

### Enhancements
* When creating coolers from unordered input, change the default temporary dir to be the same as the output file instead of the system tmp (pass '-' to use the system one). #150
* `cooler ls` and `list_coolers()` now output paths in natural order. #153
* New option in `cooler.matrix()` to handle divisive balancing weight vectors.

### Bug fixes
* Restore function of `--count-as-float` option to `cooler load`
* Fixed partitioning issue sometimes causing some bins to get split during coarsen
* `rename_chroms()` will refresh cached chromosome names #147
* `Cooler.bins()` selector will always properly convert bins/chrom integer IDs to categorical chromosome names when the number of contigs is very large and therefore the HDF5 ENUM header is missing. Before this would only happen when explicitly requesting `convert_enum=True`.

## v0.8.3 #

Date: 2019-02-11

### Bug fixes
* Fixed import bug in `rename_chroms`
* `create_cooler` no longer requires a "count" column when specifying custom value columns

## v0.8.2 ##

Date: 2019-01-20

### Enhancements

New options for `cooler dump` pixel output:
* `--matrix` option: Applies to symmetric-upper coolers; no-op for square coolers. Generates all lower triangular pixels necessary to fill the requested genomic query window. Without this option, `cooler dump` will only return the data explicity stored in the pixel table (i.e. upper triangle).
* `-one-based-ids` and `--one-based-starts` convenience options.

### Bug fixes

* A bug was introduced into the matrix-as-pixels selector in 0.8.0 that also affected `cooler dump`. The behavior has been restored to that in 0.7.


## v0.8.1 ##

Date: 2019-01-02

### Enhancements

* `cooler zoomify` command can take additional base resolutions as input.

### Bug fixes

* Fixed regression that slowed down pre-processing during coarsen.
* Fixed missing import on handling bad URIs.
* Restore but deprecate `cooler.io.ls` for backwards compatibility.

## v0.8.0 ##

Date: 2018-12-31

This is a major release from 0.7 and includes an updated format version, and several API changes and deprecations.

### Schema

* New schema version: v3
* Adds required `storage-mode` metadata attribute. Two possible values: `"symmetric-upper"` indicates a symmetric matrix encoded as upper triangle (previously the only storage mode); `"square"` indicates no special encoding (e.g. for non-symmetric matrices).

### New features

* Support for **non-symmetric** matrices, e.g. RNA-DNA maps.
    * Create function accepts a boolean `symmetric_upper` option to set the storage mode. Default is `True`.
    * Creation commands also use `symmetric_upper` by default, which can be overridden with a flag.
* All main functionality exposed through top-level functions (create, merge, coarsen, zoomify, balance)
* New commands for generic file operations and file inspection.

### API changes

* `cooler.annotate()` option `replace` now defaults to `False`.

* Submodule renaming. Old names are preserved as aliases but are deprecated.
    * `cooler.io` -> `cooler.create`.
    * `cooler.ice` -> `cooler.balance`.

* New top level public functions:
    * `cooler.create_cooler()`. Use instead of `cooler.io.create` and `cooler.io.create_from_unordered`.
    * `cooler.merge_coolers()`
    * `cooler.coarsen_cooler()`
    * `cooler.zoomify_cooler()`
    * `cooler.balance_cooler()`. Alias: `cooler.balance.iterative_correction()`.

* Refactored file operations available in `cooler.fileops`. See the API reference.

### CLI changes

* Various output options added to `cooler info`, `cooler dump`, `cooler makebins` and `cooler digest`.
* Generic data and attribute hierarchy viewers `cooler tree` and `cooler attrs`.
* Generic `cp`, `mv` and `ln` convenience commands.
* New verbosity and process info options.

### Maintenance

* Unit tests refactored and re-written for pytest.

## v0.7.11 ##

Date: 2018-08-17

* Genomic range parser supports humanized units (k/K(b), m/M(b), g/G(b))
* Experimental support for arbitrary aggregation operations in `cooler csort` (e.g. mean, median, max, min)
* Documentation updates

Bug fixes
* Fix newline handling for csort when p1 or p2 is last column.
* Fix `--count-as-float` regression in load/cload.

## v0.7.10 ##

Date: 2018-05-07

* Fix a shallow copy bug in validate pixels causing records to sometimes flip twice.
* Add ignore distance (bp) filter to cooler balance
* Start using shuffle filter by default

## v0.7.9 ##

Date: 2018-03-30

* Indexed pairs loading commands now provide option for 0- or 1-based positions (1-based by default). #115
* Fixed error introduced into cload pairix in last release.

## v0.7.8 ##

Date: 2018-03-18

### Enhancements

* New `cooler cload pairs` command provides index-free loading of pairs.
* Changed name of `create_from_unsorted` to more correct `create_from_unordered`.

### Bug fixes

* Fixed broken use of single-file temporary store in `create_from_unordered`.
* Added heuristic in pairix cload to prevent excessively large chunks. #92
* Added extra checks in `cload pairix` and `cload tabix`. #62, #75


## v0.7.7 ##

Date: 2018-03-16

### Enhancements

* Implementation of unsorted (index-free) loading
    * `cooler.io.create_from_unsorted` takes an iterable of pixel dataframe chunks that need not be properly sorted.
    * Use input sanitization procedures for pairs `sanitize_records` and binned data `sanitize_pixels` to feed data to `create_from_unsorted`. #87 #108 #109
    * The `cooler load` command is now index-free: unsorted `COO` and `BG2` input data can be streamed in. #90. This will soon be implemented as an option for loading pairs as well.
* Prevent `cooler balance` command from exiting with non-zero status upon failed convergence using convergence error policies. #93
* Improve the `create` API to support pandas read_csv-style `columns` and `dtype` kwargs to add extra value columns or override default dtypes. #108
* Experimental implementation of trans-only balancing. #56

### Bug fixes

* Fix argmax deprecation. #99


## v0.7.6 ##

Date: 2017-10-31

### Enhancements

* Cooler zoomify with explicit resolutions
* Towards standardization of multicooler structure
* Support for loading 1-based COO triplet input files

### Bug fixes

* Fixed issue of exceeding header limit with too many scaffolds. If header size is exceeded, chrom IDs are stored as raw integers instead of HDF5 enums. There should be no effect at the API level.
* Fixed issue of single-column chromosomes files not working in `cload`.
* Fixed edge case in performing joins when using both `as_pixels` and `join` options in the matrix selector.

Happy Halloween!

## v0.7.5 ##

Date: 2017-07-13

* Fix pandas issue affecting cases when loading single chromosomes
* Add transform options to higlass API

## v0.7.4 ##

Date: 2017-05-25

* Fix regression in automatic --balance option in cooler zoomify
* Fix special cases where cooler.io.create and append would not work with certain inputs

## v0.7.3 ##

Date: 2017-05-22

* Added function to print higlass zoom resolutions for a given genome and base resolution.

## v0.7.2 ##

Date: 2017-05-09

* Improve chunking and fix pickling issue with aggregating very large text datasets
* Restore zoom binsize metadata to higlass files

## v0.7.1 ##

Date: 2017-04-29

* `cooler load` command can now accept supplemental pixel fields and custom field numbers
* Fix parsing errors with unused pixel fields
* Eliminate hard dependence on dask to make pip installs simpler. Conda package will retain dask as a run time requirement.

## v0.7.0 ##

Date: 2017-04-27

### New features

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

#### API/CLI changes

* `cooler.io.create` signature changed; `chromsizes` argument is deprecated.
* `cooler csort` argument order changed

### Bug fixes

* Chromosome name length restriction removed
* `Cooler.open` function now correctly opens the specific root group of the Cooler and behaves like a proper context manager in all cases

## v0.6.6 ##

Date: 2017-03-21

* Chromosome names longer than 32 chars are forbidden for now
* Improved pairix and tabix iterators, dropped need for slow first pass over contacts

## v0.6.5 ##

Date: 2017-03-18

* Fixed pairix aggregator to properly deal with autoflipping of pairs

## v0.6.4 ##

Date: 2017-03-17

* Migrated higlass multires aggregator to `cooler coarsegrain` command
* Fixed pairix aggregator to properly deal with autoflipping of pairs

## v0.6.3 ##

Date: 2017-02-22

* Merge PairixAggregator patch from Soo.
* Update repr string
* Return matrix scale factor in balance stats rather than the bias scale factor: #35.

## v0.6.2 ##

Date: 2017-02-12

Fixed regressions in

* cooler cload tabix/pairix failed on non-fixed sized bins
* cooler show

## v0.6.1 ##

Date: 2017-02-06

* This fixes stale build used in bdist_wheel packaging that broke 0.6.0. #29

## v0.6.0 ##

Date: 2017-02-03

### Enhancements

* Dropped Python 3.3 support. Added 3.6 support.
* Added `contrib` subpackage containing utilities for higlass, including multires aggregation.
* Fixed various issues with synchronizing read/write multiprocessing with HDF5.
* Replacing prints with logging.
* Added sandboxed `tools` module to develop utilities for out-of-core algorithms using Coolers.

### New features

* Cooler objects have additional convenience properties `chromsizes`, `chromnames`.
* New file introspection functions `ls` and `is_cooler` to support nested Cooler groups.
* Cooler initializer can accept a file path and path to Cooler group.
* `cload` accepts contact lists in hiclib-style HDF5 format, the legacy tabix-indexed format, and new pairix-indexed format.

### API/CLI changes

* `create` only accepts a file path and optional group path instead of an open file object.
* `Cooler.matrix` selector now returns a balanced dense 2D NumPy array by default. Explicitly set `balance` to False to get raw counts and set `sparse` to True to get a `coo_matrix` as per old behavior.
* Command line parameters of `cload` changed significantly

### Bug fixes

* Fixed bug in `csort` that led to incorrect triangularity of trans read pairs.


## v0.5.3 ##

Date: 2016-09-10

* Check for existence of required external tools in CLI
* Fixed `cooler show` incompatibility with older versions of matplotlib
* Fixed `cooler.annotate` to work on empty dataframe input
* Fixed broken pipe signals not getting suppressed on Python 2
* `cooler cload` raises a warning when bin file lists a contig missing from the contact list


## v0.5.2 ##

Date: 2016-08-26

* Fix bug in `cooler csort` parsing of chromsizes file.
* Workaround for two locale-related issues on Python 3. Only affects cases where a machine's locale is set to ASCII or Unices which use the ambiguous C or POSIX locales.
* Fix typo in setup.py and add pysam to dependencies.


## v0.5.1 ##

Date: 2016-08-24

* Bug fix in input parser to `cooler csort`
* Update triu reording awk template in `cooler csort`
* Rename `cooler binnify` to `cooler makebins`. Binnify sounds like "aggregate" which is what `cload` does.


## v0.5.0 ##

Date: 2016-08-24

* Most scripts ported over to a new command line interface using the Click framework with many updates.
* New `show` and `info` scripts.
* Updated Readme.
* Minor bug fixes.


## v0.4.0 ##

Date: 2016-08-18

### Schema

* Updated file schema: v2
* `/bins/chroms` is now an enum instead of string column

### API changes

* Table views are a bit more intuitive: selecting field names on table view objects returns a new view on the subset of columns.
* New API function: `cooler.annotate` for doing joins

### New Features

* Support for nested Cooler "trees" at any depth in an HDF5 hierarchy
* Refactored `cooler.io` to provide "contact readers" that process different kinds of input (aggregate from a contact list, load from an existing matrix, etc.)
* Added new scripts for contact aggregation, loading, dumping and balancing


## v0.3.0 ##

Date: 2016-02-18

* 2D range selector `matrix()` now provides either rectangular data as coo_matrix or triangular data as a pixel table dataframe.
* Added binning support for any genome segmentation (i.e., fixed or variable bin width).
* Fixed issues with binning data from mapped read files.
* Genomic locus string parser now accepts ENSEMBL-style number-only chromosome names and FASTA-style sequence names containing pipes.


## v0.2.1 ##

Date: 2016-02-07

* Fixed bintable region fetcher


## v0.2 ##

Date: 2016-01-17

* First beta release


## v0.1 ##

Date: 2015-11-22

* Working initial prototype.
