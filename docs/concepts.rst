Concepts
========

Resource String
---------------
The default location for a single-cooler .cool file is the root group ``/`` of the HDF5 file. It does not need to be explicitly specified.

.. code-block:: python

    >>> import cooler
    >>> c = cooler.Cooler('data/WT.DpnII.10kb.cool')
    >>> c = cooler.Cooler('data/WT.DpnII.10kb.cool::/')  # same as above

However, coolers can be stored at any level of the HDF5 hierarchy and qualified using a URI string of the form ``/path/to/cool/file::/path/to/cooler/group``.

.. code-block:: python

    >>> c1 = cooler.Cooler('data/WT.DpnII.mcool::resolutions/10000')
    >>> c2 = cooler.Cooler('data/WT.DpnII.mcool::resolutions/1000')

The current standard for Hi-C coolers is to name multi-resolution coolers under ``.mcool`` extension,
and store differrent resolutions in an HDF5 group ``resolutions``, as shown above.

Genomic coordinates
-------------------
Cooler uses **zero-based, half-open** genomic intervals throughout, the same convention as BED files and libraries such as bioframe and pysam. An interval ``(chrom, start, end)`` covers positions ``start`` through ``end - 1``, counting from 0, and has length ``end - start``. This applies both to the ``start`` and ``end`` columns of the bin table and to any genomic range used in a query.

Genomic range strings of the form ``{chrom}:{start}-{end}`` follow the same convention: ``"chr1:10,000,000-11,000,000"`` means exactly the same interval as the tuple ``("chr1", 10_000_000, 11_000_000)``. Note that this differs from the position box of the UCSC Genome Browser, which displays one-based, fully closed coordinates: the UCSC position ``chr1:1-1,000`` corresponds to the cooler range ``chr1:0-1,000``.

For example, with 10 kb bins:

- ``"chr1:0-10k"`` selects only the first bin of ``chr1``, [0, 10000).
- ``"chr1:10k-20k"`` selects only the second bin, [10000, 20000).
- ``"chr1:9,999-10,001"`` selects the first two bins, since the interval overlaps both.

Data selection
--------------
Several :class:`cooler.Cooler` methods return data selectors. Those include selecting tables and matrices (see below). Data selectors don't retrieve any data from disk until queried. There are several ways to query using selectors. Genomic ranges may be provided as 3-tuples ``(chrom: str, start: int, end: int)`` or as range strings of the form ``{chrom}:{start}-{end}``, both interpreted as zero-based, half-open intervals (see `Genomic coordinates`_). Unit prefixes ``k, M, G`` and thousands separators are supported in range strings. For regions with start and end that are not multiples of the resolution, selectors return the shortest range of bins that fully contains the half-open interval [start, end).


Table selectors (chroms, bins, pixels)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are data selectors for the three tables: :meth:`cooler.Cooler.chroms`, :meth:`cooler.Cooler.bins`, :meth:`cooler.Cooler.pixels`.
They support the following:

- lazily select columns or lists of columns, returning new selectors
- query table rows using integer/slice indexing syntax

.. code-block:: python

    >>> c.bins()
    <cooler.core.RangeSelector1D at 0x7fdb2e4f0710>

    >>> c.bins()[:10]
    chrom    start       end    weight
    0  chr1        0   1000000       NaN
    1  chr1  1000000   2000000  1.243141
    2  chr1  2000000   3000000  1.313995
    3  chr1  3000000   4000000  1.291705
    4  chr1  4000000   5000000  1.413288
    5  chr1  5000000   6000000  1.165382
    6  chr1  6000000   7000000  0.811824
    7  chr1  7000000   8000000  1.056107
    8  chr1  8000000   9000000  1.058915
    9  chr1  9000000  10000000  1.035910

    >>> c.pixels()[:10]
       bin1_id  bin2_id  count
    0        0        0  18578
    1        0        1  11582
    2        0        2    446
    3        0        3    196
    4        0        4     83
    5        0        5    112
    6        0        6    341
    7        0        7    255
    8        0        8    387
    9        0        9    354

    >>> c.bins()['weight']
     <cooler.core.RangeSelector1D at 0x7fdb2e509240>

    >>> weights = c.bins()['weight'].fetch('chr3')
    >>> weights.head()
    494    1.144698
    495    1.549848
    496    1.212580
    497    1.097539
    498    0.871931
    Name: weight, dtype: float64

    >>> mybins1 = c.bins().fetch('chr3:10,000,000-20,000,000')
    >>> mybins2 = c.bins().fetch( ('chr3', 10000000, 20000000) )
    >>> mybins2.head()
        chrom     start       end    weight
    504  chr3  10000000  11000000  0.783160
    505  chr3  11000000  12000000  0.783806
    506  chr3  12000000  13000000  0.791204
    507  chr3  13000000  14000000  0.821171
    508  chr3  14000000  15000000  0.813079

Indexing the pixel table selects rows of the table, not bins. The pixel table stores each nonzero pixel once, sorted by ``bin1_id`` then ``bin2_id``, and with the default symmetric-upper storage only the upper triangle is stored (``bin1_id <= bin2_id``). So ``c.pixels().fetch(region)`` returns the pixels whose **first** bin is in the region as they occur in storage: the region's rows of the matrix, as stored. This includes pixels whose second bin is outside the region, even on other chromosomes. With symmetric-upper storage, these rows only hold their upper triangular part: the region's contacts with bins before it are stored in earlier rows, with the two bins swapped, so they aren't returned. With square storage, the rows are complete. To select pixels by both bins, use the matrix selector (see below):

.. code-block:: python

    >>> # Stored pixels whose first bin is on chr3
    >>> c.pixels().fetch('chr3')

    >>> # Pixels with both bins on chr3
    >>> c.matrix(balance=False, as_pixels=True).fetch('chr3')

    >>> # All contacts of chr3 with the rest of the genome, as full rows of the
    >>> # symmetric matrix (the lower triangle is filled in)
    >>> lo, hi = c.extent('chr3')
    >>> c.matrix(balance=False, sparse=True)[lo:hi, :]



Matrix selector
~~~~~~~~~~~~~~~

The :meth:`cooler.Cooler.matrix` selector supports two types of queries:

- 2D bin range queries using slice indexing syntax
- 2D genomic range range queries using the ``fetch`` method

The matrix selector’s fetch method is intended to represent a **2D range query** (rectangular window), similar to the slice semantics of a 2D array. Given a matrix selector ``sel``, when calling ``sel.fetch(region1, region2)`` the ``region1`` and ``region2`` are single contiguous genomic ranges along the first and second axes of the contact matrix. This mirrors the global slice indexing interface of the matrix selector ``sel[a:b, c:d]``, where the only difference is that the genomic range syntax cannot cross chromosome boundaries. If ``region2`` is not provided, it is taken to be the same as ``region1``. That means that ``sel.fetch('chr2:10M-20M')`` returns the same result as ``sel.fetch('chr2:10M-20M', 'chr2:10M-20M')``. As a single rectangular window, queries like ``sel.fetch('chr2', 'chr3')`` will return *inter*-chromosomal values and not intra-chromosomal ones.


.. code-block:: python

    >>> c.matrix(balance=False)[1000:1005, 1000:1005]
    array([[120022,  34107,  17335,  14053,   4137],
           [ 34107,  73396,  47427,  16125,   3642],
           [ 17335,  47427,  80458,  25105,   5394],
           [ 14053,  16125,  25105, 104536,  27214],
           [  4137,   3642,   5394,  27214, 114135]])

    >>> matrix = c.matrix(sparse=True, balance=False)
    >>> matrix
    <cooler.core.RangeSelector2D at 0x7fdb2e245908>

    >>> matrix[:]
    <3114x3114 sparse matrix of type '<class 'numpy.int64'>'
        with 8220942 stored elements in COOrdinate format>

    >>> c.matrix(balance=False, as_pixels=True, join=True)[1000:1005, 1000:1005]
       chrom1     start1       end1 chrom2     start2       end2   count
    0    chr5  115000000  116000000   chr5  115000000  116000000  120022
    1    chr5  115000000  116000000   chr5  116000000  117000000   34107
    2    chr5  115000000  116000000   chr5  117000000  118000000   17335
    3    chr5  115000000  116000000   chr5  118000000  119000000   14053
    4    chr5  115000000  116000000   chr5  119000000  120000000    4137
    5    chr5  116000000  117000000   chr5  116000000  117000000   73396
    6    chr5  116000000  117000000   chr5  117000000  118000000   47427
    7    chr5  116000000  117000000   chr5  118000000  119000000   16125
    8    chr5  116000000  117000000   chr5  119000000  120000000    3642
    9    chr5  117000000  118000000   chr5  117000000  118000000   80458
    10   chr5  117000000  118000000   chr5  118000000  119000000   25105
    11   chr5  117000000  118000000   chr5  119000000  120000000    5394
    12   chr5  118000000  119000000   chr5  118000000  119000000  104536
    13   chr5  118000000  119000000   chr5  119000000  120000000   27214
    14   chr5  119000000  120000000   chr5  119000000  120000000  114135


    >>> A1 = c.matrix().fetch('chr1')
    >>> A2 = c.matrix().fetch('chr3:10,000,000-20,000,000')
    >>> A3 = c.matrix().fetch( ('chr3', 10000000, 20000000) )
    >>> A4 = c.matrix().fetch('chr2', 'chr3')
    >>> A5 = c.matrix().fetch('chr3:10M-20M', 'chr3:35M-40M')


Create a scool file
-------------------

The creation of a single-cell cooler file is similar to a regular cooler file. Each cell needs to have a name, bin table and a pixel table.
All cells must have the same dimensions, and the bins and pixels needs to be provided as two dicts with the cell names as keys.


.. code-block:: python

    >>> name_pixel_dict = {'cell1': pixels_cell1, 'cell2': pixels_cell2, 'cell3': pixels_cell3}
    >>> name_bins_dict = {'cell1': bins_cell1, 'cell2': bins_cell2, 'cell3': bins_cell3}
    >>> cooler.create_scool('single_cell_cool.scool', name_bins_dict, name_pixel_dict)

To read the content, each individual cell must be handled as a regular cool file.

.. code-block:: python

    >> content_of_scool = cooler.fileops.list_coolers('single_cell_cool.scool')
    ['/', '/cells/cell1', '/cells/cell2', '/cells/cell3']
    >>> c1 = cooler.Cooler('single_cell_cool.scool::cells/cell1')
    >>> c2 = cooler.Cooler('single_cell_cool.scool::cells/cell2')
    >>> c3 = cooler.Cooler('single_cell_cool.scool::cells/cell3')
