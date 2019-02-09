.. _cli-reference:

CLI Reference
=============

.. toctree::
   :maxdepth: 1


Quick reference
---------------

.. program:: cooler
.. code-block:: shell

    cooler [OPTIONS] COMMAND [ARGS]...

See the cooler_cli.ipynb Jupyter Notebook for specific examples on usage: (https://github.com/mirnylab/cooler-binder).



.. rubric:: Commands

+------------------------+
| Data ingest            |
+========================+
| `cooler cload`_        |
+------------------------+
| `cooler load`_         |
+------------------------+

+------------------------+
| Reduction              |
+========================+
| `cooler merge`_        |
+------------------------+
| `cooler coarsen`_      |
+------------------------+
| `cooler zoomify`_      |
+------------------------+

+------------------------+
| Normalization          |
+========================+
| `cooler balance`_      |
+------------------------+

+------------------------+
| Export/visualization   |
+========================+
| `cooler info`_         |
+------------------------+
| `cooler dump`_         |
+------------------------+
| `cooler show`_         |
+------------------------+

+------------------------+
| File manipulation/info |
+========================+
| `cooler tree`_         |
+------------------------+
| `cooler attrs`_        |
+------------------------+
| `cooler ls`_           |
+------------------------+
| `cooler cp`_           |
+------------------------+
| `cooler mv`_           |
+------------------------+
| `cooler ln`_           |
+------------------------+

+------------------------+
| Helper commands        |
+========================+
| `cooler makebins`_     |
+------------------------+
| `cooler digest`_       |
+------------------------+
| `cooler csort`_        |
+------------------------+

.. rubric:: Options

.. option:: -v, --verbose

    Verbose logging.

.. option:: -d, --debug

    On error, drop into the post-mortem debugger shell.

.. option:: -V, --version

    Show the version and exit.

----

cooler cload
----------------

::

Create a cooler from genomic pairs and bins.

Choose a subcommand based on the format of the input contact list.

.. program:: cooler cload
.. code-block:: shell

    cooler cload [OPTIONS] COMMAND [ARGS]...

.. rubric:: Commands

.. hlist::
  :columns: 4

  * .. object:: hiclib
  * .. object:: pairix
  * .. object:: pairs
  * .. object:: tabix


----

cooler cload pairs
----------------

::

Bin any text file or stream of pairs.

Pairs data need not be sorted. Accepts compressed files.
To pipe input from stdin, set PAIRS_PATH to '-'.

BINS : One of the following

    <TEXT:INTEGER> : 1. Path to a chromsizes file, 2. Bin size in bp

    <TEXT> : Path to BED file defining the genomic bin segmentation.

PAIRS_PATH : Path to contacts (i.e. read pairs) file.

COOL_PATH : Output COOL file path or URI.

.. program:: cooler cload pairs
.. code-block:: shell

    cooler cload pairs [OPTIONS] BINS PAIRS_PATH COOL_PATH

.. rubric:: Arguments

.. option:: BINS

    Required argument

.. option:: PAIRS_PATH

    Required argument

.. option:: COOL_PATH

    Required argument

.. rubric:: Options

.. option:: --metadata <metadata>

    Path to JSON file containing user metadata.

.. option:: --assembly <assembly>

    Name of genome assembly (e.g. hg19, mm10)

.. option:: -c1, --chrom1 <chrom1>

    chrom1 field number (one-based)  [required]

.. option:: -p1, --pos1 <pos1>

    pos1 field number (one-based)  [required]

.. option:: -c2, --chrom2 <chrom2>

    chrom2 field number (one-based)  [required]

.. option:: -p2, --pos2 <pos2>

    pos2 field number (one-based)  [required]

.. option:: --chunksize <chunksize>

    Number of input lines to load at a time

.. option:: -0, --zero-based

    Positions are zero-based  [default: False]

.. option:: --comment-char <comment_char>

    Comment character that indicates lines to ignore.  [default: #]

.. option:: -N, --no-symmetric-upper

    Create a complete square matrix without implicit symmetry. This allows for distinct upper- and lower-triangle values

.. option:: --input-copy-status <input_copy_status>

    Copy status of input data when using symmetric-upper storage. | `unique`: Incoming data comes from a unique half of a symmetric map, regardless of how the coordinates of a pair are ordered. `duplex`: Incoming data contains upper- and lower-triangle duplicates. All input records that map to the lower triangle will be discarded! | If you wish to treat lower- and upper-triangle input data as distinct, use the ``--no-symmetric-upper`` option.   [default: unique]

.. option:: --field <field>

    Specify quantitative input fields to aggregate into value columns using the syntax ``--field <field-name>=<field-number>``. Optionally, append ``:`` followed by ``dtype=<dtype>`` to specify the data type (e.g. float), and/or ``agg=<agg>`` to specify an aggregation function different from sum (e.g. mean). Field numbers are 1-based. Passing 'count' as the target name will override the default behavior of storing pair counts. Repeat the ``--field`` option for each additional field.

.. option:: --temp-dir <temp_dir>

    Create temporary files in specified directory.

.. option:: --no-delete-temp

    Do not delete temporary files when finished.

.. option:: --max-merge <max_merge>

    Maximum number of chunks to merge before invoking recursive merging  [default: 200]

.. option:: --storage-options <storage_options>

    Options to modify the data filter pipeline. Provide as a comma-separated list of key-value pairs of the form 'k1=v1,k2=v2,...'. See http://docs.h5py.org/en/stable/high/dataset.html#filter-pipeline for more details.


----

cooler cload pairix
----------------

::

Bin a pairix-indexed contact list file.

BINS : One of the following

    <TEXT:INTEGER> : 1. Path to a chromsizes file, 2. Bin size in bp

    <TEXT> : Path to BED file defining the genomic bin segmentation.

PAIRS_PATH : Path to contacts (i.e. read pairs) file.

COOL_PATH : Output COOL file path or URI.

See also: 'cooler csort' to sort and index a contact list file

Pairix on GitHub: <https://github.com/4dn-dcic/pairix>.

.. program:: cooler cload pairix
.. code-block:: shell

    cooler cload pairix [OPTIONS] BINS PAIRS_PATH COOL_PATH

.. rubric:: Arguments

.. option:: BINS

    Required argument

.. option:: PAIRS_PATH

    Required argument

.. option:: COOL_PATH

    Required argument

.. rubric:: Options

.. option:: --metadata <metadata>

    Path to JSON file containing user metadata.

.. option:: --assembly <assembly>

    Name of genome assembly (e.g. hg19, mm10)

.. option:: -p, --nproc <nproc>

    Number of processes to split the work between.  [default: 8]

.. option:: -0, --zero-based

    Positions are zero-based  [default: False]

.. option:: -s, --max-split <max_split>

    Divide the pairs from each chromosome into at most this many chunks. Smaller chromosomes will be split less frequently or not at all. Increase ths value if large chromosomes dominate the workload on multiple processors.  [default: 2]


----

cooler cload tabix
----------------

::

Bin a tabix-indexed contact list file.

BINS : One of the following

    <TEXT:INTEGER> : 1. Path to a chromsizes file, 2. Bin size in bp

    <TEXT> : Path to BED file defining the genomic bin segmentation.

PAIRS_PATH : Path to contacts (i.e. read pairs) file.

COOL_PATH : Output COOL file path or URI.

See also: 'cooler csort' to sort and index a contact list file

Tabix manpage: <http://www.htslib.org/doc/tabix.html>.

.. program:: cooler cload tabix
.. code-block:: shell

    cooler cload tabix [OPTIONS] BINS PAIRS_PATH COOL_PATH

.. rubric:: Arguments

.. option:: BINS

    Required argument

.. option:: PAIRS_PATH

    Required argument

.. option:: COOL_PATH

    Required argument

.. rubric:: Options

.. option:: --metadata <metadata>

    Path to JSON file containing user metadata.

.. option:: --assembly <assembly>

    Name of genome assembly (e.g. hg19, mm10)

.. option:: -p, --nproc <nproc>

    Number of processes to split the work between.  [default: 8]

.. option:: -c2, --chrom2 <chrom2>

    chrom2 field number (one-based)

.. option:: -p2, --pos2 <pos2>

    pos2 field number (one-based)

.. option:: -0, --zero-based

    Positions are zero-based  [default: False]

.. option:: -s, --max-split <max_split>

    Divide the pairs from each chromosome into at most this many chunks. Smaller chromosomes will be split less frequently or not at all. Increase ths value if large chromosomes dominate the workload on multiple processors.  [default: 2]


----

cooler cload hiclib
----------------

::

Bin a hiclib HDF5 contact list (frag) file.

BINS : One of the following

    <TEXT:INTEGER> : 1. Path to a chromsizes file, 2. Bin size in bp

    <TEXT> : Path to BED file defining the genomic bin segmentation.

PAIRS_PATH : Path to contacts (i.e. read pairs) file.

COOL_PATH : Output COOL file path or URI.

hiclib on BitBucket: <https://bitbucket.org/mirnylab/hiclib>.

.. program:: cooler cload hiclib
.. code-block:: shell

    cooler cload hiclib [OPTIONS] BINS PAIRS_PATH COOL_PATH

.. rubric:: Arguments

.. option:: BINS

    Required argument

.. option:: PAIRS_PATH

    Required argument

.. option:: COOL_PATH

    Required argument

.. rubric:: Options

.. option:: --metadata <metadata>

    Path to JSON file containing user metadata.

.. option:: --assembly <assembly>

    Name of genome assembly (e.g. hg19, mm10)

.. option:: -c, --chunksize <chunksize>

    Control the number of pixels handled by each worker process at a time.  [default: 100000000]


----

cooler load
----------------

::

Create a cooler from a pre-binned matrix.

BINS_PATH : One of the following

    <TEXT:INTEGER> : 1. Path to a chromsizes file, 2. Bin size in bp

    <TEXT> : Path to BED file defining the genomic bin segmentation.

PIXELS_PATH : Text file containing nonzero pixel values. May be gzipped.
Pass '-' to use stdin.

COOL_PATH : Output COOL file path or URI.

**Notes**

Two input format options (tab-delimited).
Input pixel file may be compressed.

COO: COO-rdinate sparse matrix format (a.k.a. ijv triple).
3 columns: "bin1_id, bin2_id, count",

BG2: 2D version of the bedGraph format.
7 columns: "chrom1, start1, end1, chrom2, start2, end2, count"

**Examples**

cooler load -f bg2 <chrom.sizes>:<binsize> in.bg2.gz out.cool

.. program:: cooler load
.. code-block:: shell

    cooler load [OPTIONS] BINS_PATH PIXELS_PATH COOL_PATH

.. rubric:: Arguments

.. option:: BINS_PATH

    Required argument

.. option:: PIXELS_PATH

    Required argument

.. option:: COOL_PATH

    Required argument

.. rubric:: Options

.. option:: -f, --format <format>

    'coo' refers to a tab-delimited sparse triplet file (bin1, bin2, count). 'bg2' refers to a 2D bedGraph-like file (chrom1, start1, end1, chrom2, start2, end2, count).  [required]

.. option:: --metadata <metadata>

    Path to JSON file containing user metadata.

.. option:: --assembly <assembly>

    Name of genome assembly (e.g. hg19, mm10)

.. option:: --field <field>

    Add supplemental value fields or override default field numbers for the specified format. Specify quantitative input fields to aggregate into value columns using the syntax ``--field <field-name>=<field-number>``. Optionally, append ``:`` followed by ``dtype=<dtype>`` to specify the data type (e.g. float). Field numbers are 1-based. Repeat the ``--field`` option for each additional field.

.. option:: -c, --chunksize <chunksize>

    Size (in number of lines/records) of data chunks to read and process from the input file at a time. These chunks will be saved as temporary partial Coolers and merged at the end. Also specifies the size of the buffer during the merge step.

.. option:: --count-as-float

    Store the 'count' column as floating point values instead of as integers. Can also be specified using the `--field` option.

.. option:: --one-based

    Pass this flag if the bin IDs listed in a COO file are one-based instead of zero-based.

.. option:: --comment-char <comment_char>

    Comment character that indicates lines to ignore.  [default: #]

.. option:: -N, --no-symmetric-upper

    Create a complete square matrix without implicit symmetry. This allows for distinct upper- and lower-triangle values

.. option:: --input-copy-status <input_copy_status>

    Copy status of input data when using symmetric-upper storage. | `unique`: Incoming data comes from a unique half of a symmetric matrix, regardless of how element coordinates are ordered. Execution will be aborted if duplicates are detected. `duplex`: Incoming data contains upper- and lower-triangle duplicates. All lower-triangle input elements will be discarded! | If you wish to treat lower- and upper-triangle input data as distinct, use the ``--no-symmetric-upper`` option instead.   [default: unique]

.. option:: --storage-options <storage_options>

    Options to modify the data filter pipeline. Provide as a comma-separated list of key-value pairs of the form 'k1=v1,k2=v2,...'. See http://docs.h5py.org/en/stable/high/dataset.html#filter-pipeline for more details.


----

cooler merge
----------------

::

Merge multiple coolers with identical axes.

OUT_PATH : Output file path or URI.

IN_PATHS : Input file paths or URIs of coolers to merge.

**Notes**

Data columns merged:

    pixels/bin1_id, pixels/bin2_id, pixels/<value columns>

Data columns preserved:

    chroms/name, chroms/length
    bins/chrom, bins/start, bins/end

Additional columns in the the input files are not transferred to the output.

.. program:: cooler merge
.. code-block:: shell

    cooler merge [OPTIONS] OUT_PATH [IN_PATHS]...

.. rubric:: Arguments

.. option:: OUT_PATH

    Required argument

.. option:: IN_PATHS

    Optional argument(s)

.. rubric:: Options

.. option:: -c, --chunksize <chunksize>

    Size of the merge buffer in number of pixel table rows.  [default: 20000000]

.. option:: --field <field>

    Specify the names of value columns to merge as '<name>'. Repeat the `--field` option for each one. Use '<name>,dtype=<dtype>' to specify the dtype. Include ',agg=<agg>' to specify an aggregation function different from 'sum'.


----

cooler coarsen
----------------

::

Coarsen a cooler to a lower resolution.

Works by pooling *k*-by-*k* neighborhoods of pixels and aggregating.
Each chromosomal block is coarsened individually.

COOL_PATH : Path to a COOL file or Cooler URI.

.. program:: cooler coarsen
.. code-block:: shell

    cooler coarsen [OPTIONS] COOL_PATH

.. rubric:: Arguments

.. option:: COOL_PATH

    Required argument

.. rubric:: Options

.. option:: -k, --factor <factor>

    Gridding factor. The contact matrix is coarsegrained by grouping each chromosomal contact block into k-by-k element tiles  [default: 2]

.. option:: -n, -p, --nproc <nproc>

    Number of processes to use for batch processing chunks of pixels [default: 1, i.e. no process pool]

.. option:: -c, --chunksize <chunksize>

    Number of pixels allocated to each process  [default: 10000000]

.. option:: --field <field>

    Specify the names of value columns to merge as '<name>'. Repeat the `--field` option for each one. Use '<name>,dtype=<dtype>' to specify the dtype. Include ',agg=<agg>' to specify an aggregation function different from 'sum'.

.. option:: -o, --out <out>

    Output file or URI  [required]


----

cooler zoomify
----------------

::

Generate a multi-resolution cooler file by coarsening.

COOL_PATH : Path to a COOL file or Cooler URI.

.. program:: cooler zoomify
.. code-block:: shell

    cooler zoomify [OPTIONS] COOL_PATH

.. rubric:: Arguments

.. option:: COOL_PATH

    Required argument

.. rubric:: Options

.. option:: -n, -p, --nproc <nproc>

    Number of processes to use for batch processing chunks of pixels [default: 1, i.e. no process pool]

.. option:: -c, --chunksize <chunksize>

    Number of pixels allocated to each process  [default: 10000000]

.. option:: -r, --resolutions <resolutions>

    Comma-separated list of target resolutions.

.. option:: --balance

    Apply balancing to each zoom level. Off by default.

.. option:: --balance-args <balance_args>

    Additional arguments to pass to cooler balance

.. option:: -i, --base-uri <base_uri>

    One or more additional base coolers to aggregate from, if needed.

.. option:: -o, --out <out>

    Output file or URI

.. option:: --field <field>

    Specify the names of value columns to merge as '<name>'. Repeat the `--field` option for each one. Use '<name>,dtype=<dtype>' to specify the dtype. Include ',agg=<agg>' to specify an aggregation function different from 'sum'.

.. option:: --legacy

    Use the legacy layout of integer-labeled zoom levels.


----

cooler balance
----------------

::

Out-of-core matrix balancing.

Matrix must be symmetric. See the help for various filtering options to
mask out poorly mapped bins.

COOL_PATH : Path to a COOL file.

.. program:: cooler balance
.. code-block:: shell

    cooler balance [OPTIONS] COOL_PATH

.. rubric:: Arguments

.. option:: COOL_PATH

    Required argument

.. rubric:: Options

.. option:: -p, --nproc <nproc>

    Number of processes to split the work between.  [default: 8]

.. option:: -c, --chunksize <chunksize>

    Control the number of pixels handled by each worker process at a time.  [default: 10000000]

.. option:: --mad-max <mad_max>

    Ignore bins from the contact matrix using the 'MAD-max' filter: bins whose log marginal sum is less than ``mad-max`` median absolute deviations below the median log marginal sum of all the bins in the same chromosome.  [default: 5]

.. option:: --min-nnz <min_nnz>

    Ignore bins from the contact matrix whose marginal number of nonzeros is less than this number.  [default: 10]

.. option:: --min-count <min_count>

    Ignore bins from the contact matrix whose marginal count is less than this number.  [default: 0]

.. option:: --blacklist <blacklist>

    Path to a 3-column BED file containing genomic regions to mask out during the balancing procedure, e.g. sequence gaps or regions of poor mappability.

.. option:: --ignore-diags <ignore_diags>

    Number of diagonals of the contact matrix to ignore, including the main diagonal. Examples: 0 ignores nothing, 1 ignores the main diagonal, 2 ignores diagonals (-1, 0, 1), etc.  [default: 2]

.. option:: --ignore-dist <ignore_dist>

    Distance in bp to ignore.

.. option:: --tol <tol>

    Threshold value of variance of the marginals for the algorithm to converge.  [default: 1e-05]

.. option:: --max-iters <max_iters>

    Maximum number of iterations to perform if convergence is not achieved.  [default: 200]

.. option:: --cis-only

    Calculate weights against intra-chromosomal data only instead of genome-wide.

.. option:: --trans-only

    Calculate weights against inter-chromosomal data only instead of genome-wide.

.. option:: --name <name>

    Name of column to write to.  [default: weight]

.. option:: -f, --force

    Overwrite the target dataset, 'weight', if it already exists.

.. option:: --check

    Check whether a data column 'weight' already exists.

.. option:: --stdout

    Print weight column to stdout instead of saving to file.

.. option:: --convergence-policy <convergence_policy>

    What to do with weights when balancing doesn't converge in max_iters.  [default: store_final]


----

cooler info
----------------

::

Display a cooler's info and metadata.

COOL_PATH : Path to a COOL file or cooler URI.

.. program:: cooler info
.. code-block:: shell

    cooler info [OPTIONS] COOL_PATH

.. rubric:: Arguments

.. option:: COOL_PATH

    Required argument

.. rubric:: Options

.. option:: -f, --field <field>

    Print the value of a specific info field.

.. option:: -m, --metadata

    Print the user metadata in JSON format.

.. option:: -o, --out <out>

    Output file (defaults to stdout)


----

cooler dump
----------------

::

Dump a cooler's data to a text stream.

COOL_PATH : Path to COOL file or cooler URI.

.. program:: cooler dump
.. code-block:: shell

    cooler dump [OPTIONS] COOL_PATH

.. rubric:: Arguments

.. option:: COOL_PATH

    Required argument

.. rubric:: Options

.. option:: -t, --table <table>

    Which table to dump. Choosing 'chroms' or 'bins' will cause all pixel-related options to be ignored. Note that for coolers stored in symmetric-upper mode, 'pixels' only holds the upper triangle values of the matrix.  [default: pixels]

.. option:: -c, --columns <columns>

    Restrict output to a subset of columns, provided as a comma-separated list.

.. option:: -H, --header

    Print the header of column names as the first row.  [default: False]

.. option:: --na-rep <na_rep>

    Missing data representation. Default is empty ''.

.. option:: --float-format <float_format>

    Format string for floating point numbers (e.g. '.12g', '03.2f').  [default: g]

.. option:: -r, --range <range>

    The coordinates of a genomic region shown along the row dimension, in UCSC-style notation. (Example: chr1:10,000,000-11,000,000). If omitted, the entire contact matrix is printed.

.. option:: -r2, --range2 <range2>

    The coordinates of a genomic region shown along the column dimension. If omitted, the column range is the same as the row range.

.. option:: -m, --matrix

    For coolers stored in symmetric-upper mode, ensure any empty areas of the genomic query window are populated by generating the lower-triangular pixels.  [default: False]

.. option:: -b, --balanced, --no-balance

    Apply balancing weights to data. This will print an extra column called `balanced`  [default: False]

.. option:: --join

    Print the full chromosome bin coordinates instead of bin IDs. This will replace the `bin1_id` column with `chrom1`, `start1`, and `end1`, and the `bin2_id` column with `chrom2`, `start2` and `end2`.  [default: False]

.. option:: --annotate <annotate>

    Join additional columns from the bin table against the pixels. Provide a comma separated list of column names (no spaces). The merged columns will be suffixed by '1' and '2' accordingly.

.. option:: --one-based-ids

    Print bin IDs as one-based rather than zero-based.

.. option:: --one-based-starts

    Print start coordinates as one-based rather than zero-based.

.. option:: -k, --chunksize <chunksize>

    Sets the amount of pixel data loaded from disk at one time. Can affect the performance of joins on high resolution datasets. Default is to load as many rows as there are bins.

.. option:: -o, --out <out>

    Output text file If .gz extension is detected, file is written using zlib. Default behavior is to stream to stdout.


----

cooler show
----------------

::

Display and browse a cooler in matplotlib.

COOL_PATH : Path to a COOL file or Cooler URI.

RANGE : The coordinates of the genomic region to display, in UCSC notation.
Example: chr1:10,000,000-11,000,000

.. program:: cooler show
.. code-block:: shell

    cooler show [OPTIONS] COOL_PATH RANGE

.. rubric:: Arguments

.. option:: COOL_PATH

    Required argument

.. option:: RANGE

    Required argument

.. rubric:: Options

.. option:: -r2, --range2 <range2>

    The coordinates of a genomic region shown along the column dimension. If omitted, the column range is the same as the row range. Use to display asymmetric matrices or trans interactions.

.. option:: -b, --balanced

    Show the balanced contact matrix. If not provided, display the unbalanced counts.

.. option:: -o, --out <out>

    Save the image of the contact matrix to a file. If not specified, the matrix is displayed in an interactive window. The figure format is deduced from the extension of the file, the supported formats are png, jpg, svg, pdf, ps and eps.

.. option:: --dpi <dpi>

    The DPI of the figure, if saving to a file

.. option:: -s, --scale <scale>

    Scale transformation of the colormap: linear, log2 or log10. Default is log10.

.. option:: -f, --force

    Force display very large matrices (>=10^8 pixels). Use at your own risk as it may cause performance issues.

.. option:: --zmin <zmin>

    The minimal value of the color scale. Units must match those of the colormap scale. To provide a negative value use a equal sign and quotes, e.g. -zmin='-0.5'

.. option:: --zmax <zmax>

    The maximal value of the color scale. Units must match those of the colormap scale. To provide a negative value use a equal sign and quotes, e.g. -zmax='-0.5'

.. option:: --cmap <cmap>

    The colormap used to display the contact matrix. See the full list at http://matplotlib.org/examples/color/colormaps_reference.html

.. option:: --field <field>

    Pixel values to display.  [default: count]


----

cooler tree
----------------

::

Display a file's data hierarchy.

.. program:: cooler tree
.. code-block:: shell

    cooler tree [OPTIONS] URI

.. rubric:: Arguments

.. option:: URI

    Required argument

.. rubric:: Options

.. option:: -L, --level <level>


----

cooler attrs
----------------

::

Display a file's attribute hierarchy.

.. program:: cooler attrs
.. code-block:: shell

    cooler attrs [OPTIONS] URI

.. rubric:: Arguments

.. option:: URI

    Required argument

.. rubric:: Options

.. option:: -L, --level <level>


----

cooler ls
----------------

::

List all coolers inside a file.

.. program:: cooler ls
.. code-block:: shell

    cooler ls [OPTIONS] COOL_PATH

.. rubric:: Arguments

.. option:: COOL_PATH

    Required argument

.. rubric:: Options

.. option:: -l, --long

    Long listing format


----

cooler cp
----------------

::

Copy a cooler from one file to another or within the same file.

See also: h5copy, h5repack tools from HDF5 suite.

.. program:: cooler cp
.. code-block:: shell

    cooler cp [OPTIONS] SRC_URI DST_URI

.. rubric:: Arguments

.. option:: SRC_URI

    Required argument

.. option:: DST_URI

    Required argument

.. rubric:: Options

.. option:: -w, --overwrite

    Truncate and replace destination file if it already exists.


----

cooler mv
----------------

::

Rename a cooler within the same file.

.. program:: cooler mv
.. code-block:: shell

    cooler mv [OPTIONS] SRC_URI DST_URI

.. rubric:: Arguments

.. option:: SRC_URI

    Required argument

.. option:: DST_URI

    Required argument

.. rubric:: Options

.. option:: -w, --overwrite

    Truncate and replace destination file if it already exists.


----

cooler ln
----------------

::

Create a hard link to a cooler (rather than a true copy) in the same file.
Also supports soft links (in the same file) or external links (different
files).

.. program:: cooler ln
.. code-block:: shell

    cooler ln [OPTIONS] SRC_URI DST_URI

.. rubric:: Arguments

.. option:: SRC_URI

    Required argument

.. option:: DST_URI

    Required argument

.. rubric:: Options

.. option:: -w, --overwrite

    Truncate and replace destination file if it already exists.

.. option:: -s, --soft

    Creates a soft link rather than a hard link if the source and destination file are the same. Otherwise, creates an external link. This type of link uses a path rather than a pointer.


----

cooler makebins
----------------

::

Generate fixed-width genomic bins.

Output a genome segmentation at a fixed resolution as a BED file.

CHROMSIZES_PATH : UCSC-like chromsizes file, with chromosomes in desired
order.

BINSIZE : Resolution (bin size) in base pairs <int>.

.. program:: cooler makebins
.. code-block:: shell

    cooler makebins [OPTIONS] CHROMSIZES_PATH BINSIZE

.. rubric:: Arguments

.. option:: CHROMSIZES_PATH

    Required argument

.. option:: BINSIZE

    Required argument

.. rubric:: Options

.. option:: -o, --out <out>

    Output file (defaults to stdout)

.. option:: -H, --header

    Print the header of column names as the first row.  [default: False]

.. option:: -i, --rel-ids <rel_ids>

    Include a column of relative bin IDs for each chromosome. Choose whether to report them as 0- or 1-based.


----

cooler digest
----------------

::

Generate fragment-delimited genomic bins.

Output a genome segmentation of restriction fragments as a BED file.

CHROMSIZES_PATH : UCSC-like chromsizes file, with chromosomes in desired
order.

FASTA_PATH : Genome assembly FASTA file or folder containing FASTA files
(uncompressed).

ENZYME : Name of restriction enzyme

.. program:: cooler digest
.. code-block:: shell

    cooler digest [OPTIONS] CHROMSIZES_PATH FASTA_PATH ENZYME

.. rubric:: Arguments

.. option:: CHROMSIZES_PATH

    Required argument

.. option:: FASTA_PATH

    Required argument

.. option:: ENZYME

    Required argument

.. rubric:: Options

.. option:: -o, --out <out>

    Output file (defaults to stdout)

.. option:: -H, --header

    Print the header of column names as the first row.  [default: False]

.. option:: -i, --rel-ids <rel_ids>

    Include a column of relative bin IDs for each chromosome. Choose whether to report them as 0- or 1-based.


----

cooler csort
----------------

::

Sort and index a contact list.

Order the mates of each pair record so that all contacts are upper
triangular with respect to the chromosome ordering given by the chromosomes
file, sort contacts by genomic location, and index the resulting file.

PAIRS_PATH : Contacts (i.e. read pairs) text file, optionally compressed.

CHROMOSOMES_PATH : File listing desired chromosomes in the desired order.
May be tab-delimited, e.g. a UCSC-style chromsizes file. Contacts mapping to
other chromosomes will be discarded.

**Notes**

| - csort can also be used to sort and index a text representation of
|   a contact *matrix* in bedGraph-like format. In this case, substitute
|   `pos1` and `pos2` with `start1` and `start2`, respectively.
| - Requires Unix tools: sort, bgzip + tabix or pairix.

If indexing with Tabix, the output file will have the following properties:

| - Upper triangular: the read pairs on each row are assigned to side 1 or 2
|   in such a way that (chrom1, pos1) is always "less than" (chrom2, pos2)
| - Rows are lexicographically sorted by chrom1, pos1, chrom2, pos2;
|   i.e. "positionally sorted"
| - Compressed with bgzip [*]
| - Indexed using Tabix [*] on chrom1 and pos1.

If indexing with Pairix, the output file will have the following properties:

| - Upper triangular: the read pairs on each row are assigned to side 1 or 2
|   in such a way that (chrom1, pos1) is always "less than" (chrom2, pos2)
| - Rows are lexicographically sorted by chrom1, chrom2, pos1, pos2; i.e.
|   "block sorted"
| - Compressed with bgzip [*]
| - Indexed using Pairix [+] on chrom1, chrom2 and pos1.

| [*] Tabix manpage: <http://www.htslib.org/doc/tabix.html>.
| [+] Pairix on Github: <https://github.com/4dn-dcic/pairix>

.. program:: cooler csort
.. code-block:: shell

    cooler csort [OPTIONS] PAIRS_PATH CHROMOSOMES_PATH

.. rubric:: Arguments

.. option:: PAIRS_PATH

    Required argument

.. option:: CHROMOSOMES_PATH

    Required argument

.. rubric:: Options

.. option:: -c1, --chrom1 <chrom1>

    chrom1 field number in the input file (starting from 1)  [required]

.. option:: -c2, --chrom2 <chrom2>

    chrom2 field number  [required]

.. option:: -p1, --pos1 <pos1>

    pos1 field number  [required]

.. option:: -p2, --pos2 <pos2>

    pos2 field number  [required]

.. option:: -i, --index <index>

    Select the preset sort and indexing options  [default: pairix]

.. option:: --flip-only

    Only flip mates; no sorting or indexing. Write to stdout.  [default: False]

.. option:: -p, --nproc <nproc>

    Number of processors  [default: 8]

.. option:: -0, --zero-based

    Read positions are zero-based  [default: False]

.. option:: --sep <sep>

    Data delimiter in the input file  [default: \t]

.. option:: --comment-char <comment_char>

    Comment character to skip header  [default: #]

.. option:: --sort-options <sort_options>

    Quoted list of additional options to `sort` command

.. option:: -o, --out <out>

    Output gzip file

.. option:: -s1, --strand1 <strand1>

    strand1 field number (deprecated)

.. option:: -s2, --strand2 <strand2>

    strand2 field number (deprecated)


----

