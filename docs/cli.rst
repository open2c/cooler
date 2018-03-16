.. _cli-reference:

CLI Reference
=============

.. toctree::
   :maxdepth: 1


Quick reference
---------------

+------------------------+
| Helper commands        |
+========================+
| `cooler makebins`_     |
+------------------------+
| `cooler digest`_       |
+------------------------+
| `cooler csort`_        |
+------------------------+

+------------------------+
| Ingesting data         |
+========================+
| `cooler cload`_        |
+------------------------+
| `cooler load`_         |
+------------------------+

+------------------------+
| File manipulation/info |
+========================+
| `cooler list`_         |
+------------------------+
| `cooler info`_         |
+------------------------+
| `cooler copy`_         |
+------------------------+

+------------------------+
| Export/visualization   |
+========================+
| `cooler dump`_         |
+------------------------+
| `cooler show`_         |
+------------------------+

+------------------------+
| Operations             |
+========================+
| `cooler balance`_      |
+------------------------+
| `cooler merge`_        |
+------------------------+
| `cooler coarsen`_      |
+------------------------+
| `cooler zoomify`_      |
+------------------------+


cooler
------

::

    Usage: cooler [OPTIONS] COMMAND [ARGS]...
    
      Type -h or --help after any subcommand for more information.
    
    Options:
      --debug / --no-debug  Verbose logging
      -pm, --post-mortem    Post mortem debugging
      --version             Show the version and exit.
      --help                Show this message and exit.
    
    Commands:
      balance      Out-of-core contact matrix balancing.
      cload        Create a Cooler from a sorted list of...
      coarsegrain  Deprecated in favor of separate "coarsen" and...
      coarsen      Coarsen a contact matrix by uniformly...
      copy         Copy a Cooler from one file to another or...
      csort        Sort and index a contact list.
      digest       Generate fragment-delimited genomic bins.
      dump         Dump a contact matrix.
      info         Display file info and metadata.
      list         List all Coolers inside a COOL file.
      load         Load a pre-binned contact matrix into a COOL...
      makebins     Generate fixed-width genomic bins.
      merge        Merge multiple contact matrices with...
      show         Display a contact matrix.
      zoomify      Generate zoom levels for HiGlass by...


cooler makebins
----------------

::

    Usage: cooler makebins [OPTIONS] CHROMSIZES_PATH BINSIZE
    
      Generate fixed-width genomic bins. Output a genome segmentation at a fixed
      resolution as a BED file.
    
      CHROMSIZES_PATH : UCSC-like chromsizes file, with chromosomes in desired
      order.
    
      BINSIZE : Resolution (bin size) in base pairs <int>.
    
    Options:
      -o, --out TEXT  Output file (defaults to stdout)
      --help          Show this message and exit.


cooler digest
----------------

::

    Usage: cooler digest [OPTIONS] CHROMSIZES_PATH FASTA_PATH ENZYME
    
      Generate fragment-delimited genomic bins. Output a genome segmentation of
      restriction fragments as a BED file.
    
      CHROMSIZES_PATH : UCSC-like chromsizes file, with chromosomes in desired
      order.
    
      FASTA_PATH : Genome assembly FASTA file or folder containing FASTA files
      (uncompressed).
    
      ENZYME : Name of restriction enzyme
    
    Options:
      -o, --out TEXT  Output file (defaults to stdout)
      --help          Show this message and exit.


cooler csort
----------------

::

    Usage: cooler csort [OPTIONS] PAIRS_PATH CHROMOSOMES_PATH
    
      Sort and index a contact list.
    
      Order the mates of each pair record so that all contacts are upper
      triangular with respect to the chromosome ordering given by the
      chromosomes file, sort contacts by genomic location, and index the
      resulting file.
    
      Notes:
    
      - csort can also be used to sort and index a text representation of 
        a contact _matrix_ in bedGraph-like format. In this case, substitute 
        `pos1` and `pos2` with `start1` and `start2`, respectively.
      - Requires Unix tools: sort, bgzip + tabix or pairix.
    
      If indexing with Tabix, the output file will have the following
      properties:
    
      - Upper triangular: the read pairs on each row are assigned to side 1 or 2
        in such a way that (chrom1, pos1) is always "less than" (chrom2, pos2)
      - Rows are lexicographically sorted by chrom1, pos1, chrom2, pos2; 
        i.e. "positionally sorted"
      - Compressed with bgzip [*]
      - Indexed using Tabix [*] on chrom1 and pos1.
    
      If indexing with Pairix, the output file will have the following
      properties:
    
      - Upper triangular: the read pairs on each row are assigned to side 1 or 2
        in such a way that (chrom1, pos1) is always "less than" (chrom2, pos2)
      - Rows are lexicographically sorted by chrom1, chrom2, pos1, pos2; i.e. 
        "block sorted"
      - Compressed with bgzip [*]
      - Indexed using Pairix [+] on chrom1, chrom2 and pos1.
    
      [*] Tabix manpage: <http://www.htslib.org/doc/tabix.html>.
      [+] Pairix on Github: <https://github.com/4dn-dcic/pairix>
    
      Arguments:
    
      PAIRS_PATH : Contacts (i.e. read pairs) text file, optionally compressed.
    
      CHROMOSOMES_PATH : File listing desired chromosomes in the desired order.
      May be tab-delimited, e.g. a UCSC-like chromsizes file. Contacts mapping
      to other chromosomes will be discarded.
    
    Options:
      -c1, --chrom1 INTEGER       chrom1 field number in the input file (starting
                                  from 1)  [required]
      -c2, --chrom2 INTEGER       chrom2 field number  [required]
      -p1, --pos1 INTEGER         pos1 field number  [required]
      -p2, --pos2 INTEGER         pos2 field number  [required]
      -i, --index [tabix|pairix]  Select the preset sort and indexing options
                                  [default: pairix]
      --flip-only                 Only flip mates; no sorting or indexing. Write
                                  to stdout.  [default: False]
      -p, --nproc INTEGER         Number of processors  [default: 8]
      -0, --zero-based            Read positions are zero-based  [default: False]
      --sep TEXT                  Data delimiter in the input file  [default: \t]
      --comment-char TEXT         Comment character to skip header  [default: #]
      --sort-options TEXT         Quoted list of additional options to `sort`
                                  command
      -o, --out TEXT              Output gzip file
      -s1, --strand1 INTEGER      strand1 field number (deprecated)
      -s2, --strand2 INTEGER      strand2 field number (deprecated)
      --help                      Show this message and exit.


cooler cload
----------------

::

    Usage: cooler cload [OPTIONS] COMMAND [ARGS]...
    
      Create a Cooler from a sorted list of contacts and a list of genomic bins.
      Choose a subcommand based on the format of the input contact list.
    
    Options:
      --help  Show this message and exit.
    
    Commands:
      hiclib  Bin a hiclib HDF5 contact list (frag) file.
      pairix  Bin a pairix-indexed contact list file.
      tabix   Bin a tabix-indexed contact list file.        
        
        cooler cload hiclib
        ~~~~~~~~~~~~~~~~~~~
        Usage: cooler cload hiclib [OPTIONS] BINS PAIRS_PATH COOL_PATH
        
          Bin a hiclib HDF5 contact list (frag) file.
        
          BINS : One of the following
        
              <TEXT:INTEGER> : 1. Path to a chromsizes file, 2. Bin size in bp
              <TEXT> : Path to BED file defining the genomic bin segmentation.
        
          PAIRS_PATH : Path to contacts (i.e. read pairs) file.
        
          COOL_PATH : Output COOL file path.
        
          hiclib on BitBucket: <https://bitbucket.org/mirnylab/hiclib>.
        
        Options:
          --metadata TEXT          Path to JSON file containing user metadata.
          --assembly TEXT          Name of genome assembly (e.g. hg19, mm10)
          -c, --chunksize INTEGER  Control the number of pixels handled by each worker
                                   process at a time.  [default: 100000000]
          --help                   Show this message and exit.
                
        
        cooler cload pairix
        ~~~~~~~~~~~~~~~~~~~
        Usage: cooler cload pairix [OPTIONS] BINS PAIRS_PATH COOL_PATH
        
          Bin a pairix-indexed contact list file.
        
          BINS : One of the following
        
              <TEXT:INTEGER> : 1. Path to a chromsizes file, 2. Bin size in bp
              <TEXT> : Path to BED file defining the genomic bin segmentation.
        
          PAIRS_PATH : Path to contacts (i.e. read pairs) file.
        
          COOL_PATH : Output COOL file path.
        
          See also: 'cooler csort' to sort and index a contact list file
        
          Pairix on GitHub: <https://github.com/4dn-dcic/pairix>.
        
        Options:
          --metadata TEXT          Path to JSON file containing user metadata.
          --assembly TEXT          Name of genome assembly (e.g. hg19, mm10)
          -p, --nproc INTEGER      Number of processes to split the work between.
                                   [default: 8]
          -s, --max-split INTEGER  Divide the pairs from each chromosome into at most
                                   this many chunks. Smaller chromosomes will be split
                                   less frequently or not at all. Increase ths value
                                   if large chromosomes dominate the workload on
                                   multiple processors.  [default: 2]
          --help                   Show this message and exit.
                
        
        cooler cload tabix
        ~~~~~~~~~~~~~~~~~~
        Usage: cooler cload tabix [OPTIONS] BINS PAIRS_PATH COOL_PATH
        
          Bin a tabix-indexed contact list file.
        
          BINS : One of the following
        
              <TEXT:INTEGER> : 1. Path to a chromsizes file, 2. Bin size in bp
              <TEXT> : Path to BED file defining the genomic bin segmentation.
        
          PAIRS_PATH : Path to contacts (i.e. read pairs) file.
        
          COOL_PATH : Output COOL file path.
        
          See also: 'cooler csort' to sort and index a contact list file
        
          Tabix manpage: <http://www.htslib.org/doc/tabix.html>.
        
        Options:
          --metadata TEXT          Path to JSON file containing user metadata.
          --assembly TEXT          Name of genome assembly (e.g. hg19, mm10)
          -p, --nproc INTEGER      Number of processes to split the work between.
                                   [default: 8]
          -c2, --chrom2 INTEGER    chrom2 field number (one-based)
          -p2, --pos2 INTEGER      pos2 field number (one-based)
          -s, --max-split INTEGER  Divide the pairs from each chromosome into at most
                                   this many chunks. Smaller chromosomes will be split
                                   less frequently or not at all. Increase ths value
                                   if large chromosomes dominate the workload on
                                   multiple processors.  [default: 2]
          --help                   Show this message and exit.
        


cooler load
----------------

::

    Usage: cooler load [OPTIONS] BINS_PATH PIXELS_PATH COOL_PATH
    
      Load a pre-binned contact matrix into a COOL file.
    
      Two input format options (tab-delimited):
    
      * COO: COO-rdinate sparse matrix format (a.k.a. ijv triple). 3 columns.
    
      - columns: "bin1_id, bin2_id, count",
    
      * BG2: 2D version of the bedGraph format. 7 columns.
    
      - columns: "chrom1, start1, end1, chrom2, start2, end2, count"
    
      Input pixel file may be compressed.
    
      **New in v0.7.7: Input files no longer need to be sorted or indexed!**
    
      Example:
    
      cooler load -f bg2 <chrom.sizes>:<binsize> in.bg2.gz out.cool
    
      Arguments:
    
      BINS_PATH : One of the following
    
          <TEXT:INTEGER> : 1. Path to a chromsizes file, 2. Bin size in bp
          <TEXT> : Path to BED file defining the genomic bin segmentation.
    
      PIXELS_PATH : Text file containing nonzero pixel values. May be gzipped.
      Pass '-' to use stdin.
    
      COOL_PATH : Output COOL file path
    
    Options:
      -f, --format [coo|bg2]        'coo' refers to a tab-delimited sparse triplet
                                    file (bin1, bin2, count). 'bg2' refers to a 2D
                                    bedGraph-like file (chrom1, start1, end1,
                                    chrom2, start2, end2, count).  [required]
      --metadata TEXT               Path to JSON file containing user metadata.
      --assembly TEXT               Name of genome assembly (e.g. hg19, mm10)
      --field TEXT                  Add supplemental value fields or override
                                    default field numbers for the specified
                                    format. Specify as '<name>,<number>' or as
                                    '<name>,<number>,<dtype>' to enforce a dtype
                                    other than `float` or the default for a
                                    standard column. Field numbers are 1-based.
                                    Repeat the `--field` option for each
                                    additional field. [Changed in v0.7.7: use a
                                    comma separator, rather than a space.]
      -c, --chunksize INTEGER       Size (in number of lines/records) of data
                                    chunks to read and process from the input file
                                    at a time. These chunks will be saved as
                                    temporary partial Coolers and merged at the
                                    end. Also specifies the size of the buffer
                                    during the merge step.
      --count-as-float              Store the 'count' column as floating point
                                    values instead of as integers. Can also be
                                    specified using the `--field` option.
      --one-based                   Pass this flag if the bin IDs listed in a COO
                                    file are one-based instead of zero-based.
      --comment-char TEXT           Comment character that indicates lines to
                                    ignore.  [default: #]
      --tril-action [reflect|drop]  How to handle lower triangle pixels.
                                    'reflect': make lower triangle pixels upper
                                    triangular. Use this if your input data comes
                                    only from a unique half of a symmetric matrix
                                    (but may not respect the specified chromosome
                                    order).'drop': discard all lower triangle
                                    pixels. Use this if your input data is derived
                                    from a complete symmetric matrix.  [default:
                                    reflect]
      --help                        Show this message and exit.


cooler list
----------------

::

    Usage: cooler list [OPTIONS] COOL_PATH
    
      List all Coolers inside a COOL file.
    
    Options:
      -l, --long  Long listing format
      --help      Show this message and exit.


cooler info
----------------

::

    Usage: cooler info [OPTIONS] COOL_PATH
    
      Display file info and metadata.
    
      COOL_PATH : Path to a COOL file or Cooler URI.
    
    Options:
      -f, --field TEXT  Print the value of a specific info field.
      -m, --metadata    Print the user metadata in JSON format.
      -o, --out TEXT    Output file (defaults to stdout)
      --help            Show this message and exit.


cooler copy
----------------

::

    Usage: cooler copy [OPTIONS] SRC_URI DST_URI
    
      Copy a Cooler from one file to another or within the same file.
    
      See also: h5copy, h5repack tools from HDF5 suite
    
      Arguments:
    
      SRC_URI : Path to source file or URI to source Cooler group
    
      DST_URI : Path to destination file or URI to destination Cooler group
    
    Options:
      -w, --overwrite  Truncate and replace destination file if it already exists.
      -l, --link       If the source and destination file are the same, create a
                       hard link to the source group instead of a true copy.
      -m, --rename     If the source and destination file are the same, create a
                       hard link to the source group and remove the original
                       reference.
      -s, --soft-link  If the source and destination file are the same, create a
                       soft link. If the destination file is different, create an
                       external link. This type of link uses a path rather than a
                       pointer.
      --help           Show this message and exit.


cooler dump
----------------

::

    Usage: cooler dump [OPTIONS] COOL_PATH
    
      Dump a contact matrix. Print the contents of a COOL file to tab-delimited
      text.
    
      COOL_PATH : Path to COOL file or Cooler URI.
    
    Options:
      -t, --table [chroms|bins|pixels]
                                      Which table to dump. Choosing 'chroms' or
                                      'bins' will cause all pixel-related options
                                      to be ignored. Note that dumping 'pixels'
                                      will only provide data for the upper
                                      triangle of the contact matrix.   [default:
                                      pixels]
      --header                        Print the header of column names as the
                                      first row.  [default: False]
      -k, --chunksize INTEGER         Sets the amount of pixel data loaded from
                                      disk at one time. Can affect the performance
                                      of joins on high resolution datasets.
                                      Default is to load as many rows as there are
                                      bins.
      -r, --range TEXT                The coordinates of a genomic region shown
                                      along the row dimension, in UCSC notation.
                                      (Example: chr1:10,000,000-11,000,000). If
                                      omitted, the entire contact matrix is
                                      printed.
      -r2, --range2 TEXT              The coordinates of a genomic region shown
                                      along the column dimension. If omitted, the
                                      column range is the same as the row range.
      -b, --balanced / --no-balance   Apply balancing weights to data. This will
                                      print an extra column called `balanced`
                                      [default: False]
      --join                          Print the full chromosome bin coordinates
                                      instead of bin IDs. This will replace the
                                      `bin1_id` column with `chrom1`, `start1`,
                                      and `end1`, and the `bin2_id` column with
                                      `chrom2`, `start2` and `end2`.  [default:
                                      False]
      --annotate TEXT                 Join additional columns from the bin table
                                      against the pixels. Provide a comma
                                      separated list of column names (no spaces).
                                      The merged columns will be suffixed by '1'
                                      and '2' accordingly.
      -o, --out TEXT                  Output text file If .gz extension is
                                      detected, file is written using zlib.
                                      Default behavior is to stream to stdout.
      --help                          Show this message and exit.


cooler show
----------------

::

    Usage: cooler show [OPTIONS] COOL_PATH RANGE
    
      Display a contact matrix. Display a region of a contact matrix stored in a
      COOL file.
    
      Arguments:
    
      COOL_PATH : Path to a COOL file or Cooler URI.
    
      RANGE : The coordinates of the genomic region to display, in UCSC
      notation. Example: chr1:10,000,000-11,000,000
    
    Options:
      -r2, --range2 TEXT              The coordinates of a genomic region shown
                                      along the column dimension. If omitted, the
                                      column range is the same as the row range.
                                      Use to display asymmetric matrices or trans
                                      interactions.
      -b, --balanced                  Show the balanced contact matrix. If not
                                      provided, display the unbalanced counts.
      -o, --out TEXT                  Save the image of the contact matrix to a
                                      file. If not specified, the matrix is
                                      displayed in an interactive window. The
                                      figure format is deduced from the extension
                                      of the file, the supported formats are png,
                                      jpg, svg, pdf, ps and eps.
      --dpi INTEGER                   The DPI of the figure, if saving to a file
      -s, --scale [linear|log2|log10]
                                      Scale transformation of the colormap:
                                      linear, log2 or log10. Default is log10.
      -f, --force                     Force display very large matrices (>=10^8
                                      pixels). Use at your own risk as it may
                                      cause performance issues.
      --zmin FLOAT                    The minimal value of the color scale. Units
                                      must match those of the colormap scale. To
                                      provide a negative value use a equal sign
                                      and quotes, e.g. -zmin='-0.5'
      --zmax FLOAT                    The maximal value of the color scale. Units
                                      must match those of the colormap scale. To
                                      provide a negative value use a equal sign
                                      and quotes, e.g. -zmax='-0.5'
      --cmap TEXT                     The colormap used to display the contact
                                      matrix. See the full list at http://matplotl
                                      ib.org/examples/color/colormaps_reference.ht
                                      ml
      --field TEXT                    Pixel values to display.  [default: count]
      --help                          Show this message and exit.


cooler balance
----------------

::

    Usage: cooler balance [OPTIONS] COOL_PATH
    
      Out-of-core contact matrix balancing.
    
      Assumes uniform binning. See the help for various filtering options to
      ignore poorly mapped bins.
    
      COOL_PATH : Path to a COOL file.
    
    Options:
      -p, --nproc INTEGER             Number of processes to split the work
                                      between.  [default: 8]
      -c, --chunksize INTEGER         Control the number of pixels handled by each
                                      worker process at a time.  [default:
                                      10000000]
      --mad-max INTEGER               Ignore bins from the contact matrix using
                                      the 'MAD-max' filter: bins whose log
                                      marginal sum is less than ``mad-max`` median
                                      absolute deviations below the median log
                                      marginal sum of all the bins in the same
                                      chromosome.  [default: 5]
      --min-nnz INTEGER               Ignore bins from the contact matrix whose
                                      marginal number of nonzeros is less than
                                      this number.  [default: 10]
      --min-count INTEGER             Ignore bins from the contact matrix whose
                                      marginal count is less than this number.
                                      [default: 0]
      --blacklist PATH                Path to a 3-column BED file containing
                                      genomic regions to mask out during the
                                      balancing procedure, e.g. sequence gaps or
                                      regions of poor mappability.
      --ignore-diags INTEGER          Number of diagonals of the contact matrix to
                                      ignore, including the main diagonal.
                                      Examples: 0 ignores nothing, 1 ignores the
                                      main diagonal, 2 ignores diagonals (-1, 0,
                                      1), etc.  [default: 2]
      --tol FLOAT                     Threshold value of variance of the marginals
                                      for the algorithm to converge.  [default:
                                      1e-05]
      --max-iters INTEGER             Maximum number of iterations to perform if
                                      convergence is not achieved.  [default: 200]
      --cis-only                      Calculate weights against intra-chromosomal
                                      data only instead of genome-wide.
      --trans-only                    Calculate weights against inter-chromosomal
                                      data only instead of genome-wide.
      --name TEXT                     Name of column to write to.  [default:
                                      weight]
      -f, --force                     Overwrite the target dataset, 'weight', if
                                      it already exists.
      --check                         Check whether a data column 'weight' already
                                      exists.
      --stdout                        Print weight column to stdout instead of
                                      saving to file.
      --convergence-policy [store_final|store_nan|discard|error]
                                      What to do with weights when balancing
                                      doesn't converge in max_iters.  [default:
                                      store_final]
      --help                          Show this message and exit.


cooler merge
----------------

::

    Usage: cooler merge [OPTIONS] OUT_PATH [IN_PATHS]...
    
      Merge multiple contact matrices with identical axes.
    
      Data columns merged:
    
          pixels/bin1_id, pixels/bin2_id, pixels/count
    
      Data columns preserved:
    
          chroms/name, chroms/length     bins/chrom, bins/start, bins/end
    
      Additional columns in the the input files are not preserved in the output.
    
    Options:
      -c, --chunksize INTEGER  [default: 20000000]
      --help                   Show this message and exit.


cooler coarsen
----------------

::

    Usage: cooler coarsen [OPTIONS] COOL_PATH
    
      Coarsen a contact matrix by uniformly gridding the elements of each
      chromosomal block and summing the elements inside the grid tiles, i.e. a
      2-D histogram.
    
      Arguments:
    
      COOL_PATH : Path to a COOL file or Cooler URI.
    
    Options:
      -k, --factor INTEGER     Gridding factor. The contact matrix is
                               coarsegrained by grouping each chromosomal contact
                               block into k-by-k element tiles  [default: 2]
      -n, -p, --nproc INTEGER  Number of processes to use for batch processing
                               chunks of pixels [default: 1, i.e. no process pool]
      -c, --chunksize INTEGER  Number of pixels allocated to each process
                               [default: 10000000]
      -o, --out TEXT           Output file or URI  [required]
      --help                   Show this message and exit.


cooler zoomify
----------------

::

    Usage: cooler zoomify [OPTIONS] COOL_PATH
    
      Generate zoom levels for HiGlass by recursively generating 2-by-2 element
      tiled aggregations of the contact matrix until reaching a minimum
      dimension. The aggregations are stored in a multi-resolution file.
    
      Arguments:
    
      COOL_PATH : Path to a COOL file or Cooler URI.
    
    Options:
      -n, -p, --nproc INTEGER   Number of processes to use for batch processing
                                chunks of pixels [default: 1, i.e. no process
                                pool]
      -c, --chunksize INTEGER   Number of pixels allocated to each process
                                [default: 10000000]
      --balance / --no-balance  Apply balancing to each zoom level  [default:
                                False]
      --balance-args TEXT       Additional arguments to pass to cooler balance
      -o, --out TEXT            Output file or URI
      -r, --resolutions TEXT    Comma-separated list of target resolutions
      --help                    Show this message and exit.


