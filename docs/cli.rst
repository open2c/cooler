.. _cli-reference:

CLI Reference
=============

.. toctree::
   :maxdepth: 1


cooler
------

::

   Usage: cooler [OPTIONS] COMMAND [ARGS]...

   Options:
     -h, --help  Show this message and exit.

   Commands:
     balance   Out-of-core contact matrix balancing.
     cload     Aggregate and load a sorted contact list.
     csort     Sort and index a contact list.
     digest    Make fragment-delimited genomic bins.
     dump      Dump a contact matrix.
     info      Display file info and metadata.
     load      Load a contact matrix.
     makebins  Make fixed-width genomic bins.
     show      Display a contact matrix.



cooler makebins
---------------

::

    Usage: cooler makebins [OPTIONS] CHROMSIZES_PATH BINSIZE

      Make fixed-width genomic bins. Output a genome segmentation at a fixed
      resolution as a BED file.

      CHROMSIZES_PATH : UCSC-like chromsizes file, with chromosomes in desired
      order.

      BINSIZE : Resolution (bin size) in base pairs <int>.

    Options:
      -o, --out TEXT  Output file (defaults to stdout)
      -h, --help      Show this message and exit.


cooler digest
-------------

::

    Usage: cooler digest [OPTIONS] CHROMSIZES_PATH FASTA_PATH ENZYME

      Make fragment-delimited genomic bins. Output a genome segmentation of
      restriction fragments as a BED file.

      CHROMSIZES_PATH : UCSC-like chromsizes file, with chromosomes in desired
      order.

      FASTA_PATH : Genome assembly FASTA file or folder containing FASTA files
      (uncompressed).

      ENZYME : Name of restriction enzyme

    Options:
      -o, --out TEXT  Output file (defaults to stdout)
      -h, --help      Show this message and exit.


cooler csort
------------

::

    Usage: cooler csort [OPTIONS] CHROMSIZES_PATH PAIRS_PATH

      Sort and index a contact list. Arrange the reads of each pair so that all
      contacts are upper triangular with respect to the chromosome ordering
      given by the chromsizes file, and sort contacts by position.

      Requires Unix tools: pigz, sort, bgzip, tabix

      CHROMSIZES_PATH : UCSC-like chromsizes file, with chromosomes in desired
      order.

      PAIRS_PATH : Contacts (i.e. read pairs) text file.

      The output file will have the following properties:

      - Upper triangular: the read pairs on each row are assigned to side 1 or 2
        in such a way that (chrom1, pos1) is always "less than" (chrom2, pos2),
        according to the desired chromosome order as given by the chromsizes
        file.
      - Rows are lexically sorted by chrom1, pos1, chrom2, pos2. Here, the way
        chromosomes are sorted does not need to respect the desired order.
      - Compressed with bgzip [*]
      - Indexed using Tabix [*] on chrom1 and pos1: `tabix -0 -s1 -b2 -e2`

      [*] Tabix manpage: <http://www.htslib.org/doc/tabix.html>.

    Options:
      -c1, --chrom1 INTEGER   chrom1 field number in the input file (starting from
                              1)
      -p1, --pos1 INTEGER     pos1 field number
      -s1, --strand1 INTEGER  strand1 field number
      -c2, --chrom2 INTEGER   chrom2 field number
      -p2, --pos2 INTEGER     pos2 field number
      -s2, --strand2 INTEGER  strand2 field number
      -p, --nproc INTEGER     number of processors
      -o, --out TEXT          Output gzip file
      -h, --help              Show this message and exit.



cooler cload
------------

::

    Usage: cooler cload [OPTIONS] BINS_PATH PAIRS_PATH COOL_PATH

      Aggregate and load a sorted contact list. Create a COOL file from a list
      of contacts and a list of bins.

      BINS_PATH : Path to BED file defining the genomic bin segmentation.

      PAIRS_PATH : Path to contacts (i.e. read pairs) file whose first six
      columns are `chrom1`, `pos1`, `strand1`, `chrom2`, `pos2`, `strand2`. The
      contacts file must be:

      - Tab-delimited
      - Upper triangular: reads on each row are oriented such that
        (chrom1, pos1) is "less than" (chrom2, pos2) according to the
        desired chromosome ordering
      - Lexically sorted by chrom1, pos1, chrom2, pos2. Here, the way
        chromosomes are ordered is not crucial because of indexing (below).
      - Compressed with bgzip [*]
      - Indexed using Tabix [*] on chrom1 and pos1: `tabix -0 -s1 -b2 -e2`

      COOL_PATH : Output COOL file path.

      See also: 'cooler csort' to sort and index a contact list file

      [*] Tabix manpage: <http://www.htslib.org/doc/tabix.html>.

    Options:
      --metadata TEXT  Path to JSON file containing user metadata.
      --assembly TEXT  Name of genome assembly (e.g. hg19, mm10)
      -h, --help       Show this message and exit.



cooler load
-----------

::

    Usage: cooler load [OPTIONS] BINS_PATH PIXELS_PATH COOL_PATH

      Load a contact matrix. Load a sparse-formatted text dump of a contact
      matrix into a COOL file.

      BINS_PATH : BED-like file containing genomic bin segmentation

      PIXELS_PATH : Three-column, sorted sparse matrix text file in ijv-triple
      format, a.k.a. COO format. May be gzipped. Must be lexically sorted by
      bin1_id and bin2_id.

      COOL_PATH : Output COOL file path

    Options:
      --metadata TEXT  Path to JSON file containing user metadata.
      --assembly TEXT  Name of genome assembly (e.g. hg19, mm10)
      -h, --help       Show this message and exit.



cooler balance
---------------

::

    Usage: cooler balance [OPTIONS] COOL_PATH

      Out-of-core contact matrix balancing.

      Assumes uniform binning. See the help for various filtering options to
      ignore poorly mapped bins.

      COOL_PATH : Path to a COOL file.

    Options:
      -p, --nproc INTEGER      Number of processes to split the work between.
                               [default: 8]
      -c, --chunksize INTEGER  Control the number of pixels handled by each worker
                               process at a time.  [default: 100000000]
      --mad-max INTEGER        Ignore bins from the contact matrix using the 'MAD-
                               max' filter: bins whose log marginal sum is less
                               than ``mad-max`` mean absolute deviations below the
                               median log marginal sum of all the bins in the same
                               chromosome.  [default: 3]
      --min-nnz INTEGER        Ignore bins from the contact matrix whose marginal
                               number of nonzeros is less than this number.
                               [default: 10]
      --min-count INTEGER      Ignore bins from the contact matrix whose marginal
                               count is less than this number.  [default: 0]
      --ignore-diags INTEGER   Number of diagonals of the contact matrix to
                               ignore, including the main diagonal. Examples: 0
                               ignores nothing, 1 ignores the main diagonal, 2
                               ignores diagonals (-1, 0, 1), etc.  [default: 3]
      --tol FLOAT              Threshold value of variance of the marginals for
                               the algorithm to converge.  [default: 1e-05]
      --max-iters INTEGER      Maximum number of iterations to perform if
                               convergence is not achieved.  [default: 200]
      --cis-only               Calculate weights against intra-chromosomal data
                               only instead of genome-wide.
      -f, --force              Overwrite the target dataset, 'weight', if it
                               already exists.
      -h, --help               Show this message and exit.


cooler info
-----------

::

    Usage: cooler info [OPTIONS] COOL_PATH

      Display file info and metadata.

      COOL_PATH : Path to a COOL file.

    Options:
      -f, --field TEXT  Print the value of a specific info field.
      -m, --metadata    Print the user metadata in JSON format.
      -o, --out TEXT    Output file (defaults to stdout)
      -h, --help        Show this message and exit.


cooler dump
-----------

::

    Usage: cooler dump [OPTIONS] COOL_PATH

      Dump a contact matrix. Print the contents of a COOL file to tab-delimited
      text.

      COOL_PATH : Path to COOL file containing a contact matrix.

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
      -b, --balanced                  Apply balancing weights to data. This will
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
      -h, --help                      Show this message and exit.




cooler show
-----------

::

    Usage: cooler show [OPTIONS] COOLER_PATH RANGE

      Display a contact matrix. Display a region of a contact matrix stored in a
      COOL file.

      COOLER_PATH : Path to a COOL file

      RANGE : The coordinates of the genomic region to display, in UCSC
      notation. Example: chr1:10,000,000-11,000,000

    Options:
      -r2, --range2 TEXT              The coordinates of a genomic region shown
                                      along the column dimension. If omitted, the
                                      column range is the same as the row range.
                                      Use to display asymmetric matrices or trans
                                      interactions.
      -o, --out TEXT                  Save the image of the contact matrix to a
                                      file. If not specified, the matrix is
                                      displayed in an interactive window. The
                                      figure format is deduced from the extension
                                      of the file, the supported formats are png,
                                      jpg, svg, pdf, ps and eps.
      --dpi INTEGER                   The DPI of the figure, if saving to a file
      --raw                           Show the raw (i.e. unbalanced) contact
                                      matrix. If not provided, display the
                                      balanced contact matrix.
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
      -h, --help                      Show this message and exit.


