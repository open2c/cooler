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
