.. _api-reference:

API Reference
=============

.. toctree::
   :maxdepth: 1


Quick reference
---------------

.. autosummary:: 
    cooler.Cooler
    cooler.get
    cooler.info
    cooler.chroms
    cooler.bins
    cooler.pixels
    cooler.matrix
    cooler.annotate

.. autosummary:: 
    cooler.io.create
    cooler.io.is_cooler
    cooler.io.ls
    cooler.io.ContactReader

.. autosummary:: 
    cooler.ice.iterative_correction


Data selectors
~~~~~~~~~~~~~~

Several :class:`cooler.Cooler` methods return data selectors. They don't retrieve any data from disk until queried. There are several ways to query using selectors.


Table selectors (chroms, bins, pixels)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- lazy select columns or lists of columns, returning new selectors
- query rows using integer/slice indexing syntax
- *bins* supports fetching genomic ranges using ``fetch`` method
- *pixels* supports fetching genomic ranges along the *bin1* axis

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

    >>> mybins1 = c.bins().fetch('chr3:10,000,000-20,000,000')
    >>> mybins2 = c.bins().fetch( ('chr3', 10000000, 20000000) )



Matrix selector
^^^^^^^^^^^^^^^

- 2D bin range queries using slice indexing syntax
- 2D genomic range range queries using the ``fetch`` method


.. code-block:: python

    >>> matrix = c.matrix(sparse=True, balance=False)
    >>> matrix
    <cooler.core.RangeSelector2D at 0x7fdb2e245908>

    >>> matrix[:]
    <3114x3114 sparse matrix of type '<class 'numpy.int64'>'
        with 8220942 stored elements in COOrdinate format>

    >>> c.matrix(balance=False)[1000:1005, 1000:1005]
    array([[120022,  34107,  17335,  14053,   4137],
           [ 34107,  73396,  47427,  16125,   3642],
           [ 17335,  47427,  80458,  25105,   5394],
           [ 14053,  16125,  25105, 104536,  27214],
           [  4137,   3642,   5394,  27214, 114135]])

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


Genomic intervals can be provided using UCSC-style strings or chrom-start-end triples ``(str, int, int)``.


cooler
------

.. autoclass:: cooler.Cooler
    :members:
.. autofunction:: cooler.get
.. autofunction:: cooler.info
.. autofunction:: cooler.chroms
.. autofunction:: cooler.bins
.. autofunction:: cooler.pixels
.. autofunction:: cooler.matrix
.. autofunction:: cooler.annotate


cooler.io
---------

.. automodule:: cooler.io
    :members:
    :undoc-members:
    :show-inheritance:
.. autoclass:: cooler.io.ContactReader
    :members: __iter__

cooler.ice
----------

.. autofunction:: cooler.ice.iterative_correction


cooler.util
-----------

.. automodule:: cooler.util
    :members:
    :undoc-members:
    :show-inheritance:
