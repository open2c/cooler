# Release notes #

## Unreleased

### New features
* `cooler.log_to_stderr()` prints cooler's log messages to stderr, for scripts and notebooks that don't configure logging.
* The 2D range query engines (`cooler.core.DirectRangeQuery2D`, `FillLowerRangeQuery2D` and `CSRReader`) accept a list of value columns, and return them all in frame output. Matrix output still needs a single column.
* `cooler cload pairs --field <name>=<col>:agg=mean`, and `create_cooler(..., ordered=False, agg={...: "mean"})`, now work for inputs larger than one chunk. Means are computed over all records of a pixel (ignoring missing values) by carrying sums and counts through the external sort. Other aggregations that can't be combined across chunks (e.g. `median`) still raise an error.
* `cooler balance --rescale-marginals/--no-rescale-marginals` (#448). With `--no-rescale-marginals`, the balanced matrix's rows sum to the mean marginal sum of the unbalanced matrix, stored as the `scale` attribute.
* `balance_cooler(..., ignore_dist=...)` ignores intra-chromosomal pixels whose bin midpoints are less than a given distance apart, which also works with variable-size bins.

### Performance
* Range queries filter pixels in a single vectorized pass instead of looping over matrix rows in Python. The pixel filtering step is about 7-16x faster, and matrix queries on HDF5 files about 1.5x faster end to end, with peak memory about the same as before. Queries spanning many rows benefit most.
* Cis-only balancing (`cooler balance --cis-only`) balances every chromosome in the same passes over the pixels, instead of separate passes for each chromosome, each of which also read the whole bin table. For assemblies with many contigs, it's as fast as genome-wide balancing: 20 iterations on a cooler with 2,010 chromosomes took 4.6 s instead of 316 s. The weights are unchanged. Chromosomes that don't converge are named in a single warning, and progress is logged once per iteration.

### Deprecations
* `cooler.set_verbosity_level` and `cooler.get_verbosity_level` are deprecated and will be removed in cooler 0.13. Use `cooler.log_to_stderr(level)` or `logging.getLogger("cooler").setLevel(level)`. Verbosity levels remain a CLI concept (`-v`).
* `cooler cload pairs --input-copy-status` is deprecated and will be removed in cooler 0.13; pairs should list each contact once. Passing it logs a warning. `cooler load --input-copy-status` (matrix input) is unchanged, with a clearer description of `duplex`.

### Removed
* `cooler.sandbox`, including the experimental `cooler.sandbox.dask.read_table`, which no longer worked with current cooler or dask. For dask dataframes of pixels, the range query engines have `to_dask_frame()`.
* `cooler zoomify --legacy` and `legacy_zoomify`, which wrote multi-resolution files in the layout used before cooler 0.8, with integer-labeled zoom levels (`::0`, `::1`, ...). Use the default `resolutions/<binsize>` layout; files in the old layout can still be read.

### Bug fixes
* fix: `Cooler.info` and `cooler attrs` no longer parse text attributes as JSON. For example, an assembly named `"123"` was returned as the integer 123. Attributes that cooler defines as text (`genome-assembly`, `generated-by`, `creation-date`, `format`, `format-url`, `bin-type`, `storage-mode`) are returned as stored. Other string attributes are still parsed as JSON, including `metadata` and numbers stored as strings by older versions.
* fix: `list_coolers`, `is_cooler` and `cooler ls` recognize coolers in format version 2 or earlier that have no `format` attribute.
* fix: A missing space in the error message for missing balancing weights.
* fix: Importing cooler no longer configures logging. Previously, it attached a handler to the `cooler` logger that printed warnings to stdout and stopped records from propagating, so they bypassed the user's logging configuration and mixed into piped output (including from worker processes started with `spawn` or `forkserver`). The `cooler` logger now only has a `NullHandler`, as recommended for libraries: records propagate to the root logger, and are silent unless logging is configured. The warning about records on chromosomes missing from the bin table is now issued with `warnings.warn`, so it is shown either way.
* fix: The CLI sets up logging for the duration of each command and restores the previous state afterwards. Running commands repeatedly in the same process (e.g. with click's `CliRunner`) no longer logs to a stale stream, and `logging.captureWarnings` is no longer left on. Warnings are logged as one line each (`WARNING  cooler: ConvergenceWarning: ...`), instead of Python's multi-line warning output. The process information printed with `-vv` is now logged when the command finishes.
* fix: `cooler dump` writes every column of the pixel table by default, and `-c/--columns` selects pixel columns. Previously, only `count` could be dumped, and `-c` was ignored for pixels. Unknown column names raise an error.
* fix: `cooler dump --one-based-ids` works on its own. It only took effect with `--join`, `--balanced` or `--annotate`.
* fix: `cooler dump` closes its output file and the cooler file when finished.
* fix: `cooler ls` accepts a URI and lists only the coolers under that group (e.g. `cooler ls data.mcool::resolutions`).
* fix: Selecting a bin table column whose name contains "chrom" (e.g. `bins()["chrom_arm"]`) no longer decodes it as chromosome names.
* fix: Region strings with decimal units are parsed exactly: `"4.1Mb"` gave 4,099,999 because of float rounding. Values that aren't a whole number of base pairs (`"1.0001kb"`) raise an error, and so do region strings with trailing text (`"chr1:100-200-300"`) or more than one `:`. For chromosome names containing `:`, use a `(chrom, start, end)` tuple.
* fix: Bin tables must tile each chromosome: the bins of a chromosome are consecutive rows, start at 0, and have no gaps or overlaps (#254). Previously, such tables were accepted, and queries with variable-size bins could return bins from the previous chromosome. Queries on existing files with such tables no longer do.
* fix: Out-of-range indices on table and matrix selectors raise `IndexError` instead of returning the wrong rows (e.g. `bins()[-n - 5]`). Slices are clamped like Python sequences, so a slice past the end is truncated and `m[10:5]` is empty. Float indices and lists of integers raise `TypeError`; a list selects columns by name.
* fix: Single-column table selectors (e.g. `bins()["start"]`) support `columns`, `dtypes`, `keys()` and `in`.
* fix: Bin and chromosome tables with extra string columns can be written; this crashed. Appending to a table with `cooler.core.put` maps categorical values to the existing column's categories instead of writing each chunk's own codes, and raises instead of silently truncating strings longer than the existing column allows.
* fix: Files whose string attributes or datasets are stored as bytes or variable-length strings (written by other tools) are read correctly. A bytes `storage-mode` made symmetric files be queried as square, without the lower triangle.
* fix: `Cooler`, `open_hdf5` and URI parsing accept `pathlib.Path` and other path-like objects.
* fix: `annotate` accepts unsigned 64-bit bin IDs.
* fix: Pixel pipelines (`cooler.parallel.split`) can be pickled and copied.
* fix: Coarse zoom levels whose chromosomes are each a single bin are stored with a fixed bin size; they were stored as variable-size bins. `cooler.create.create` gains a `binsize` parameter for when the bin size can't be inferred from the bin table, which coarsening and merging now pass.
* fix: The `Cooler.matrix` docs give the right `join` default and divisive weight names (`VC_SQRT`, not `SQRT_VC`), and describe NaN values with balancing: dense output is NaN in every row and column of a bin with a NaN weight, while sparse output only has NaN values at stored pixels.
* fix: `cooler load` rejects duplicate pixels wherever they occur in the input. Previously, duplicates in the same chunk raised an error but duplicates in different chunks were silently summed, so the result depended on `--chunksize`. `create_cooler(..., ordered=False)` gains a `duplicates` option: `"combine"` (default) combines repeated pixels with `agg`, now also within a chunk; `"error"` raises.
* fix: Creating a cooler into a new or truncated file (mode "w") writes to a temporary file that replaces the destination only when complete, so a failure leaves an existing file intact and no partial file behind.
* fix: The `fillvalue` and `scaleoffset` storage options work (they were applied to string datasets and failed), and `maxshape`, which could never work, is no longer accepted.
* fix: `create_cooler(..., symmetric_upper=False)` no longer warns about `triucheck` unless it was explicitly set to True.
* fix: `cooler makebins` and `cooler digest` close their output files.
* fix: `coarsen_cooler` (and `cooler coarsen`/`cooler zoomify --field`) kept only the `count` column, silently dropping other value columns, and failed if `count` wasn't selected. All requested columns are now coarsened.
* fix: `zoomify_cooler` with several base coolers (`cooler zoomify --base-uri`) kept only one of them, and could coarsen a base resolution from a finer one, overwriting its data. Base resolutions are now all kept as given; base coolers with the same resolution are rejected.
* fix: `cooler zoomify -r <N>b` (binary progression from N) crashed.
* fix: `cooler zoomify --balance` balances every level in the output, including base resolutions, and recognizes a column name set with `--balance-args "--name ..."`.
* fix: Coarsening and zoomifying keep the input's `genome-assembly` and `metadata`.
* fix: Merges no longer buffer one bin more than `mergebuf` per epoch (which also triggered spurious buffering warnings), and the merger no longer leaves input files open.
* fix: `cooler coarsen -k` requires a factor of at least 2 (a factor of 1 failed with an `AssertionError`); `coarsen_cooler` no longer modifies the caller's `dtypes`.
* fix: The `--field` help of `cooler coarsen` and `cooler merge` showed the wrong syntax (`<name>:dtype=...`, not `<name>,dtype=...`).
* fix: Balancing raises an error for coolers stored as a square matrix (`symmetric_upper=False`), which may hold asymmetric matrices that balancing with a single weight vector doesn't apply to. Previously, it gave wrong weights (off by a factor of √2 even for symmetric data) and halved the `min_nnz`/`min_count` thresholds.
* fix: `--ignore-diags` in balancing only ignores intra-chromosomal pixels. Previously, it also dropped trans pixels between the last bins of one chromosome and the first bins of the next (#211), so weights of bins at chromosome ends change slightly.
* fix: `cooler balance --ignore-dist` works with variable-size bins (#464).
* fix: Bins with no data left after filtering get NaN weights. Previously they could get finite weights with `--min-nnz 0`, or when a whole chromosome was empty, since the MAD filter skipped chromosomes without data.
* fix: Balancing reports non-convergence (and a NaN `scale`) when every bin is filtered out. Previously it reported convergence. In cis-only balancing, an empty chromosome still counts as converged, unless all of them are empty.
* fix: `balance_cooler(..., max_iters=0)` no longer raises an `UnboundLocalError`. `x0` is no longer modified and returned as the weights, integer `x0` works, and an `x0` of the wrong length raises an error.
* fix: `cooler balance --blacklist` ignores regions on chromosomes that aren't in the cooler (e.g. chrY or chrM), with a warning, instead of crashing.
* fix: Storing balancing weights no longer fails when the stats attributes exceed HDF5's 64 KiB limit for attributes, as with cis-only balancing of assemblies with many thousands of contigs (#226). New columns are created with an object header that stores large attributes separately.
* fix: Balancing stats record whether `trans_only` was used. The `--mad-max` help describes the filter correctly: marginals are normalized by their chromosome's median, then compared with a genome-wide cutoff.
* fix: `cooler balance --force` and `balance_cooler(..., store=True)` replace an existing weight column only once the new weights are written. Previously the old weights were deleted before balancing, so a discarded result (`--convergence-policy discard`/`error`) or an error lost them.
* fix: `create_scool` honors `ordered=False` (and `mergebuf`, `max_merge`, `temp_dir`, `delete_temp`) for cells given as iterables of unsorted chunks, checks that the cell names of a per-cell bins dict match the pixels dict, and checks that per-cell bin tables match the shared one. A failure while writing a cell leaves the cells already written intact.
* fix: Creating a cooler from sorted input (`create_cooler(..., ordered=True)` and the functions built on it) checks that pixels are sorted and unique across all chunks, not just within each one. Previously, out-of-order chunks silently produced a cooler with a corrupt index.
* fix: A failure while creating a cooler no longer leaves a partial cooler at the destination, or removes an existing cooler being replaced: coolers are written into a temporary group and moved into place when complete.
* fix: Coolers with no pixels no longer contain placeholder rows in their pixel datasets.
* fix: `cooler.create.create` raises `TypeError` for unknown keyword arguments instead of silently ignoring them. This exposed options that `create_scool` passed but that had no effect (`ordered`, `mergebuf`, `delete_temp`, `temp_dir`, `max_merge`).
* fix: `cooler cload` rejects positions at or beyond the end of a chromosome. Previously, a position equal to the chromosome length (one past the end, after conversion to zero-based) was assigned to the first bin of the next chromosome, turning a cis contact into a trans contact, or silently dropped on the last chromosome. `cload pairs` now also keeps pixel bounds checking on, and `cload tabix`/`pairix` check the second position.
* fix: `cooler cload pairs --field` and `cooler load --field` no longer swap the values of fields given out of column order (e.g. `--field b=9 --field a=8`).
* fix: `cooler cload pairs --field <name>=<col>:agg=<agg>` combined the per-chunk results by summing them, so `max`, `min` and other aggregations were wrong for inputs larger than `--chunksize`. Pixels are now combined with the requested aggregation; aggregations that can't be combined across chunks (e.g. `mean`) raise an error when the input spans more than one chunk. `create_cooler` gains an `agg` parameter for this.
* fix: `create_cooler(..., ordered=False)` sorts each input chunk, so unsorted chunks no longer produce unsorted coolers with duplicate pixels.
* fix: `create_cooler(..., ordered=False)` and `cload pairs` handle empty input (e.g. no records on the requested chromosomes) by writing an empty cooler instead of crashing. `--max-merge` values of 1 or 2 no longer crash.
* fix: Chromosomes longer than 2,147,483,647 bp (the int32 maximum) are stored correctly. Chromosome lengths and bin coordinates were always written as int32, so they silently overflowed to negative values. They are now written as int64 when a chromosome needs it, and as int32 otherwise, as before (#480).
* fix: Pixel values that don't fit in the column's integer type (e.g. counts above 2^31 - 1 with the default int32) raise an error suggesting int64, instead of being silently clipped.
* fix: `create_scool` sorts DataFrame pixels, rejects cell names that would be stored under the same name, and works with a non-root group URI.
* fix: `cooler cload pairs --comment-char` is now used, and `.zip` input works.
* fix: `cooler load --field` no longer drops the `count` column.
* fix: Malformed bin (BED) files raise a parse error instead of an `AttributeError`.
* fix: Records on chromosomes not in the bin table no longer rely on a deprecated pandas behavior that will raise in pandas 4.
* fix: `cooler cp`/`mv`/`ln` (and `cooler.fileops.cp/mv/ln`) no longer modify or truncate the destination file when the operation fails, e.g. on a mistyped source path or a hard link between different files. Previously, `-w/--overwrite` could empty the destination file.
* fix: `-w/--overwrite` replaces only the destination group or dataset; the whole file is replaced only when the destination is the root of another file. It now also works within a single file.
* fix: `mv` between two different files raises an error instead of silently copying. Copying or linking a group into itself, or moving a group onto one of its ancestors, is rejected.
* fix: Same-file operations are detected by the file's real path, so `cp ./a.cool a.cool::/x` works.
* fix: External links (`ln -s` across files) store the target path relative to the file containing the link, so they resolve from any working directory.
* fix: `list_coolers`, `cooler ls`, `cooler tree` and `cooler attrs` report objects under the link paths used to reach them, skip dangling links, and no longer recurse forever on link cycles.
* fix: `cooler attrs` no longer crashes on array-valued attributes.
* fix: `is_cooler` returns False for a missing group, and `is_scool_file` for a missing or non-HDF5 file, as documented. `list_scool_cells` lists only the cells under `/cells`.
* fix: `to_dask_frame()` on range queries works with current dask versions; it relied on dask's legacy `DataFrame` constructor, which newer dask releases removed. By @nvictus in https://github.com/open2c/cooler/pull/501

### Maintenance
* CLI log messages have timestamps and a more readable format: `2026-10-26 14:03:12 INFO     cooler.cli.balance: Balancing ...` (previously `INFO:cooler.cli.balance:Balancing ...`).
* The balancing and reduction loggers are named `cooler.balance` and `cooler.reduce` (previously `cooler._balance` and `cooler._reduce`).
* Log messages are clearer: balancing logs the iteration number with the variance, merging and zoomifying summarize what they do, and internal details (chunk ranges, merge epochs, input fields, shell commands in `csort`) moved to DEBUG. `cooler balance` logs non-convergence once, as a warning unless `--convergence-policy error` is used. Messages show the destination of a cooler being written, not its temporary file.
* Removed the internal `set_logging_context`, `get_logging_context` and `configure` logging helpers.
* Fixed deprecation warnings from NumPy (generic timedelta unit), Click (`CliRunner.isolated_filesystem`) and pandas in the code and tests.
* `cooler cload` logs a warning naming chromosomes whose records are skipped because they are not in the bin table (once per chromosome).
* Removed the unused, untested `SparseBlockLoader` and `ArrayBlockLoader` classes, and compatibility code for pandas < 1.2.
* The minimum h5py version is raised to 3.0 (from 2.5). h5py 3.0 is the oldest release with wheels for every supported Python version.
* Removed the undocumented `cooler.create.append` function. Columns are now added to existing tables through an internal column writer (`cooler.create._writer._ColumnWriter`), which writes several columns at once in batches of rows (e.g. from polars `LazyFrame.sink_batches`), and adds or replaces them only once every row is written. A public interface will follow with the new API.
* `create_scool` is built on an internal scool writer (`cooler.create._writer._ScoolWriter`) that writes the shared tables once and adds cells one at a time (`add_cell`), replacing the scool-specific options of `cooler.create.create` (`append_scool`, `scool_root_uri`).
* `create_from_unordered` is built on an internal external-sort writer (`cooler.create._writer._UnorderedCoolerWriter`), which validates its options up front and removes its temporary files on error.
* Cooler creation is built on an internal chunked writer (`cooler.create._writer._CoolerWriter`) that validates, writes and finalizes a cooler one chunk at a time, keeping the output file open between chunks unless a lock is shared with readers.
* The `to_sparse_array()` method of range query engines returns a SciPy `coo_array` instead of a pydata/sparse `COO`, so it no longer needs the optional `sparse` package. `to_sparse_matrix()` now converts the same `coo_array` to a `coo_matrix`. The minimum SciPy version is raised to 1.11, the oldest compatible with NumPy 1.26.
* Property-based tests of 2D range queries against dense NumPy references. By @nvictus in https://github.com/open2c/cooler/pull/501

### Documentation
* The `agg` parameter of `merge_coolers`, `coarsen_cooler` and `zoomify_cooler`, and the `--field` help of `cooler merge`, `coarsen` and `zoomify`, describe how aggregations treat pixels that aren't stored: they don't count as zeros, so `mean` averages over stored pixels only (or, when merging, over the coolers that store the pixel). Zoomify coarsens each resolution from the previous one, so `mean` there is an unweighted mean of means.

## [v0.11.0](https://github.com/open2c/cooler/compare/v0.10.4...v0.11.0)

### New features
* `cooler show` and `cooler dump` accept `--weight-name` to select a balancing weight column other than `weight` (implies `--balanced`). By @dmalzl and @nvictus in https://github.com/open2c/cooler/pull/472
* `cooler.parallel.get_lock()` returns the global multiprocessing lock, creating it on first use. By @nvictus in https://github.com/open2c/cooler/pull/497

### Bug fixes
* fix: Count main diagonal pixels once when computing balancing marginals. This only affects balancing with `ignore_diags=0` (`--ignore-diags 0`), whose weights will change. By @jpaggi in https://github.com/open2c/cooler/pull/471
* fix: `cooler balance --blacklist` no longer fails on single-line or empty BED files, and skips header, track and browser lines. By @ShigrafS in https://github.com/open2c/cooler/pull/462
* fix: `read_chromsizes` raises a clear error for malformed rows (e.g. spaces instead of tabs) instead of failing later with a cryptic error. By @ShigrafS in https://github.com/open2c/cooler/pull/458
* fix: Appending cells to an scool file with `create_scool(..., mode="a")` reuses the existing bin table instead of writing a new copy each time, raises an error if the bin table differs, and keeps the `ncells` count correct. By @Snowman-cpu and @nvictus in https://github.com/open2c/cooler/pull/457
* fix: Wait for worker processes to exit after parallel `balance`, `coarsen`/`zoomify` and `cload` runs, so no lingering workers can hold HDF5 file locks. By @nvictus in https://github.com/open2c/cooler/pull/495
* fix: `cooler` can be imported in environments without process-shared semaphores, such as AWS Lambda (#329). The global lock in `cooler.parallel` is now created on first use; `cooler.parallel.lock` still works, and `cooler.parallel.get_lock()` is the preferred accessor. By @nvictus in https://github.com/open2c/cooler/pull/497
* fix: `cooler show` no longer blurs heatmaps by interpolating between pixels.
* fix: The process info logged at exit with `-vv` is limited to a small set of attributes that are stable across psutil versions and platforms. It no longer includes environment variables, which could leak secrets into logs. By @nvictus in https://github.com/open2c/cooler/pull/498

### Behavior changes
* Balancing weights computed with `ignore_diags=0` will differ from previous versions (#471).
* `read_chromsizes` raises an error for malformed rows that were previously dropped or misparsed (#458).
* Appending cells with a different bin table to an existing scool file now raises an error (#457).

### Maintenance
* Removed the `simplejson` dependency in favor of the standard library `json` module. By @a-detiste in https://github.com/open2c/cooler/pull/454
* Development tooling moved to uv and PEP 735 dependency groups, with a committed `uv.lock`. The `cooler[test]`, `cooler[dev]` and `cooler[docs]` extras are removed; use `uv sync` or `pip install -e .[all] --group dev` instead. https://github.com/open2c/cooler/pull/489, https://github.com/open2c/cooler/pull/491
* Upgraded ruff to v0.17 and fixed new lint findings, including a missing `raise` in `preferred_sequence`. https://github.com/open2c/cooler/pull/489

### Documentation
* Documented that genomic coordinates and range strings are zero-based and half-open. https://github.com/open2c/cooler/pull/496
* Clarified which compression formats `cooler cload pairs` accepts. By @sapunyangkut in https://github.com/open2c/cooler/pull/483

## [v0.10.4](https://github.com/open2c/cooler/compare/v0.10.3...v0.10.4)

### Bug fixes
* fix: Fix indexing in `merge_breakpoints` by @conradbzura in https://github.com/open2c/cooler/pull/468

## [v0.10.3](https://github.com/open2c/cooler/compare/v0.10.2...v0.10.3)

### Maintenance
* Update to support numpy>=1.26 and numpy 2.x by @lrauschning in https://github.com/open2c/cooler/pull/446
* Dropped support for numpy < 1.26 and Python 3.8.
* Dropping pypairix dependency from `cooler[all]` until compilation issues are resolved.

## [v0.10.2](https://github.com/open2c/cooler/compare/v0.10.1...v0.10.2)

### Maintenance
* NumPy 2.0 was released. Pin `numpy < 2` until we achieve compatibility.

## [v0.10.1](https://github.com/open2c/cooler/compare/v0.10.0...v0.10.1)

### Bug fixes
* fix: Pixel annotation with partial bin table in https://github.com/open2c/cooler/pull/426

## [v0.10.0](https://github.com/open2c/cooler/compare/v0.9.3...v0.10.0)

### New features
* Support custom contig separator characters in pairix index by @nvictus in https://github.com/open2c/cooler/pull/398

### Performance improvements and bug fixes
* Make annotate faster by using indexing instead of merge by @Phlya in https://github.com/open2c/cooler/pull/353
* Fix input text file handling on certain pandas versions by @thomas-reimonn in https://github.com/open2c/cooler/pull/384
* Reduce memory consumption of very high-resolution merges by @nvictus in https://github.com/open2c/cooler/pull/408

### Maintenance
* Add windows and macos platforms to CI by @nvictus in https://github.com/open2c/cooler/pull/383
* Add type annotations by @GarrettNg in https://github.com/open2c/cooler/pull/387
* docs: Replace ftp with OSN storage bucket for example data. by @nvictus in https://github.com/open2c/cooler/pull/399

### New Contributors
* @a-detiste made their first contribution in https://github.com/open2c/cooler/pull/346
* @manzt made their first contribution in https://github.com/open2c/cooler/pull/363
* @thomas-reimonn made their first contribution in https://github.com/open2c/cooler/pull/380
* @GarrettNg made their first contribution in https://github.com/open2c/cooler/pull/387

## [v0.9.3](https://github.com/open2c/cooler/compare/v0.9.2...v0.9.3)

### Bug fixes

* Fix estimation of mean bin size when using variable-sized bins #332.
* Fix regression bug to handle multiple convergence statuses in balance CLI #334.

## [v0.9.2](https://github.com/open2c/cooler/compare/v0.9.1...v0.9.2)

### Bug fixes

Several fixes from @robomics
* Improve handling of uint attributes
* Fix incorrect variance stored by _balance_cisonly
* cload.py: fix detection of pandas v2.*.* causing IOHandles error

## [v0.9.1](https://github.com/open2c/cooler/compare/v0.9.0...v0.9.1)

Date 2023-01-23

### Maintenance

* Export `_IndexingMixin` from `cooler.core` to keep private API used by cooltools.

## [v0.9.0](https://github.com/open2c/cooler/compare/v0.8.11...v0.9.0)

Date 2023-01-18

### New features

* New query engine implementation (no user-facing API changes).
* New logging context and verbosity system.

### API changes

* `cooler.balance_cooler` now uses the same default argument values as the CLI.

### Schema

* `cooler balance` now inserts `divisive_weights: False` metadata to balancing weights.

### Maintenance

* Dropped Python 2 support.
* Removed deprecated `io` module.
* Deprecated `tools` module; renamed `parallel`.
* Various dependency maintenance upgrades.
* Modernized CI and dev tools.
* Migrated to pyproject.toml and modernized packaging infra.

## [v0.8.11](https://github.com/open2c/cooler/compare/v0.8.10...v0.8.11)

Date : 2021-04-01

### Bug fixes

* Hotfix `cooler cload pairs` file header parsing to work with the latest version of pandas (>=1.2).
* Update CLI help for `--balance-args` to clarify usage. By @gfudenberg.

## [v0.8.10](https://github.com/open2c/cooler/compare/v0.8.9...v0.8.10)

Date : 2020-09-25

### Bug fixes
* Fixed the new header parsing in `cooler cload pairs` to handle esoteric file stream implementations. Specifically `GzipFile` had stopped working. By @golobor

## [v0.8.9](https://github.com/open2c/cooler/compare/v0.8.8...v0.8.9)

Date : 2020-07-17

### Enhancements
* Added single-cell cooler file flavor (.scool) (#201)


## [v0.8.8](https://github.com/open2c/cooler/compare/v0.8.7...v0.8.8)

Date : 2020-06-23

### Maintenance
* Improved code coverage
* Added missing autodoc for cooler balance
* Dropped pysam and biopython as hard dependencies
* Officially sunsetting Python 2.7 support

### Enhancements
* Added zoom progressions (#203)

### Bug fixes
* Allow hashes in read IDs in cload pairs (#193)


## [v0.8.7](https://github.com/open2c/cooler/compare/v0.8.6...v0.8.7)

Date: 2020-01-12

### Maintenance
* Code styling with black
* Add coverage reporting

### Bug fixes
* Replace `json` with `simplejson` to deal with attrs stored as bytes


## [v0.8.6](https://github.com/open2c/cooler/compare/v0.8.5...v0.8.6)

Date: 2019-08-12

### Maintenance
* Added contributing guidelines

### Bug fixes
* Fixed a related regression that affected selection of the `chrom` column.

Post-release `v0.8.6.post0`: requirements files added to MANIFEST.in


## [v0.8.5](https://github.com/open2c/cooler/compare/v0.8.4...v0.8.5)

Date: 2019-04-08

### Bug fixes
* Fixed a regression that prevented selection of bins excluding the `chrom` column.


## [v0.8.4](https://github.com/open2c/cooler/compare/v0.8.3...v0.8.4)

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


## [v0.8.3](https://github.com/open2c/cooler/compare/v0.8.2...v0.8.3)

Date: 2019-02-11

### Bug fixes
* Fixed import bug in `rename_chroms`
* `create_cooler` no longer requires a "count" column when specifying custom value columns


## [v0.8.2](https://github.com/open2c/cooler/compare/v0.8.1...v0.8.2)

Date: 2019-01-20

### Enhancements

New options for `cooler dump` pixel output:
* `--matrix` option: Applies to symmetric-upper coolers; no-op for square coolers. Generates all lower triangular pixels necessary to fill the requested genomic query window. Without this option, `cooler dump` will only return the data explicity stored in the pixel table (i.e. upper triangle).
* `-one-based-ids` and `--one-based-starts` convenience options.

### Bug fixes

* A bug was introduced into the matrix-as-pixels selector in 0.8.0 that also affected `cooler dump`. The behavior has been restored to that in 0.7.


## [v0.8.1](https://github.com/open2c/cooler/compare/v0.8.0...v0.8.1)

Date: 2019-01-02

### Enhancements

* `cooler zoomify` command can take additional base resolutions as input.

### Bug fixes

* Fixed regression that slowed down pre-processing during coarsen.
* Fixed missing import on handling bad URIs.
* Restore but deprecate `cooler.io.ls` for backwards compatibility.


## [v0.8.0](https://github.com/open2c/cooler/compare/v0.7.11...v0.8.0)

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


## [v0.7.11](https://github.com/open2c/cooler/compare/v0.7.10...v0.7.11)

Date: 2018-08-17

* Genomic range parser supports humanized units (k/K(b), m/M(b), g/G(b))
* Experimental support for arbitrary aggregation operations in `cooler csort` (e.g. mean, median, max, min)
* Documentation updates

Bug fixes
* Fix newline handling for csort when p1 or p2 is last column.
* Fix `--count-as-float` regression in load/cload.


## [v0.7.10](https://github.com/open2c/cooler/compare/v0.7.9...v0.7.10)

Date: 2018-05-07

* Fix a shallow copy bug in validate pixels causing records to sometimes flip twice.
* Add ignore distance (bp) filter to cooler balance
* Start using shuffle filter by default


## [v0.7.9](https://github.com/open2c/cooler/compare/v0.7.8...v0.7.9)

Date: 2018-03-30

* Indexed pairs loading commands now provide option for 0- or 1-based positions (1-based by default). #115
* Fixed error introduced into cload pairix in last release.


## [v0.7.8](https://github.com/open2c/cooler/compare/v0.7.7...v0.7.8)

Date: 2018-03-18

### Enhancements

* New `cooler cload pairs` command provides index-free loading of pairs.
* Changed name of `create_from_unsorted` to more correct `create_from_unordered`.

### Bug fixes

* Fixed broken use of single-file temporary store in `create_from_unordered`.
* Added heuristic in pairix cload to prevent excessively large chunks. #92
* Added extra checks in `cload pairix` and `cload tabix`. #62, #75


## [v0.7.7](https://github.com/open2c/cooler/compare/v0.7.6...v0.7.7)

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


## [v0.7.6](https://github.com/open2c/cooler/compare/v0.7.5...v0.7.6)

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


## [v0.7.5](https://github.com/open2c/cooler/compare/v0.7.4...v0.7.5)

Date: 2017-07-13

* Fix pandas issue affecting cases when loading single chromosomes
* Add transform options to higlass API


## [v0.7.4](https://github.com/open2c/cooler/compare/v0.7.3...v0.7.4)

Date: 2017-05-25

* Fix regression in automatic --balance option in cooler zoomify
* Fix special cases where cooler.io.create and append would not work with certain inputs


## [v0.7.3](https://github.com/open2c/cooler/compare/v0.7.2...v0.7.3)

Date: 2017-05-22

* Added function to print higlass zoom resolutions for a given genome and base resolution.


## [v0.7.2](https://github.com/open2c/cooler/compare/v0.7.1...v0.7.2)

Date: 2017-05-09

* Improve chunking and fix pickling issue with aggregating very large text datasets
* Restore zoom binsize metadata to higlass files


## [v0.7.1](https://github.com/open2c/cooler/compare/v0.7.0...v0.7.1)

Date: 2017-04-29

* `cooler load` command can now accept supplemental pixel fields and custom field numbers
* Fix parsing errors with unused pixel fields
* Eliminate hard dependence on dask to make pip installs simpler. Conda package will retain dask as a run time requirement.


## [v0.7.0](https://github.com/open2c/cooler/compare/v0.6.6...v0.7.0)

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


## [v0.6.6](https://github.com/open2c/cooler/compare/v0.6.5...v0.6.6)

Date: 2017-03-21

* Chromosome names longer than 32 chars are forbidden for now
* Improved pairix and tabix iterators, dropped need for slow first pass over contacts


## [v0.6.5](https://github.com/open2c/cooler/compare/v0.6.4...v0.6.5)

Date: 2017-03-18

* Fixed pairix aggregator to properly deal with autoflipping of pairs


## [v0.6.4](https://github.com/open2c/cooler/compare/v0.6.3...v0.6.4)

Date: 2017-03-17

* Migrated higlass multires aggregator to `cooler coarsegrain` command
* Fixed pairix aggregator to properly deal with autoflipping of pairs


## [v0.6.3](https://github.com/open2c/cooler/compare/v0.6.2...v0.6.3)

Date: 2017-02-22

* Merge PairixAggregator patch from Soo.
* Update repr string
* Return matrix scale factor in balance stats rather than the bias scale factor: #35.


## [v0.6.2](https://github.com/open2c/cooler/compare/v0.6.1...v0.6.2)

Date: 2017-02-12

Fixed regressions in

* cooler cload tabix/pairix failed on non-fixed sized bins
* cooler show


## [v0.6.1](https://github.com/open2c/cooler/compare/v0.6.0...v0.6.1)

Date: 2017-02-06

* This fixes stale build used in bdist_wheel packaging that broke 0.6.0. #29


## [v0.6.0](https://github.com/open2c/cooler/compare/v0.5.3...v0.6.0)

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


## [v0.5.3](https://github.com/open2c/cooler/compare/v0.5.2...v0.5.3)

Date: 2016-09-10

* Check for existence of required external tools in CLI
* Fixed `cooler show` incompatibility with older versions of matplotlib
* Fixed `cooler.annotate` to work on empty dataframe input
* Fixed broken pipe signals not getting suppressed on Python 2
* `cooler cload` raises a warning when bin file lists a contig missing from the contact list


## [v0.5.2](https://github.com/open2c/cooler/compare/v0.5.1...v0.5.2)

Date: 2016-08-26

* Fix bug in `cooler csort` parsing of chromsizes file.
* Workaround for two locale-related issues on Python 3. Only affects cases where a machine's locale is set to ASCII or Unices which use the ambiguous C or POSIX locales.
* Fix typo in setup.py and add pysam to dependencies.


## [v0.5.1](https://github.com/open2c/cooler/compare/v0.5.0...v0.5.1)

Date: 2016-08-24

* Bug fix in input parser to `cooler csort`
* Update triu reording awk template in `cooler csort`
* Rename `cooler binnify` to `cooler makebins`. Binnify sounds like "aggregate" which is what `cload` does.


## [v0.5.0](https://github.com/open2c/cooler/compare/v0.4.0...v0.5.0)

Date: 2016-08-24

* Most scripts ported over to a new command line interface using the Click framework with many updates.
* New `show` and `info` scripts.
* Updated Readme.
* Minor bug fixes.


## [v0.4.0](https://github.com/open2c/cooler/compare/v0.3.0...v0.4.0)

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


## [v0.3.0](https://github.com/open2c/cooler/compare/v0.2.1...v0.3.0)

Date: 2016-02-18

* 2D range selector `matrix()` now provides either rectangular data as coo_matrix or triangular data as a pixel table dataframe.
* Added binning support for any genome segmentation (i.e., fixed or variable bin width).
* Fixed issues with binning data from mapped read files.
* Genomic locus string parser now accepts ENSEMBL-style number-only chromosome names and FASTA-style sequence names containing pipes.


## [v0.2.1](https://github.com/open2c/cooler/compare/v0.2...v0.2.1)

Date: 2016-02-07

* Fixed bintable region fetcher


## [v0.2](https://github.com/open2c/cooler/compare/v0.1...v0.2)

Date: 2016-01-17

* First beta release


## [v0.1](https://github.com/open2c/cooler/releases/tag/v0.1)

Date: 2015-11-22

* Working initial prototype.
