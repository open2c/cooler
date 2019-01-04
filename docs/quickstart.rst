Quickstart
==========


Installation
------------

Install :py:mod:`cooler`  from `PyPI <https://pypi.org/project/cooler>`_ using :command:`pip`.

::

    $ pip install cooler

Requirements:

- Python 2.7 or 3.4 and higher
- libhdf5
- Python packages :py:mod:`numpy`, :py:mod:`scipy`, :py:mod:`pandas`, :py:mod:`h5py`. 

We highly recommend using the conda package manager to install scientific packages like these. To get :command:`conda`, you can download either the full `Anaconda <https://www.continuum.io/downloads>`_ Python distribution which comes with lots of data science software or the minimal `Miniconda <http://conda.pydata.org/miniconda.html>`_ distribution which is just the standalone package manager plus Python. In the latter case, you can install the packages as follows:

::

    $ conda install numpy scipy pandas h5py

If you are using conda, you can alternatively install cooler from the `bioconda channel <https://bioconda.github.io>`_.

::

    $ conda install -c conda-forge -c bioconda cooler


Command line interface
----------------------

See:

- Jupyter Notebook `CLI walkthrough <https://github.com/mirnylab/cooler-binder/blob/master/cooler_cli.ipynb>`_.
- The `CLI Reference <http://cooler.readthedocs.io/en/latest/cli.html>`_ for more information.


The :py:mod:`cooler` package includes command line tools for creating, querying and manipulating cooler files.

::

    $ cooler cload pairs hg19.chrom.sizes:10000 $PAIRS_FILE out.10000.cool
    $ cooler balance -p 10 out.10000.cool
    $ cooler dump -b -t pixels --header --join -r chr3:10M-12M -r2 chr17 out.10000.cool | head

Output:

::

    chrom1  start1  end1    chrom2  start2  end2    count   balanced
    chr3    10000000        10010000        chr17   0       10000   1       0.810766
    chr3    10000000        10010000        chr17   520000  530000  1       1.2055
    chr3    10000000        10010000        chr17   640000  650000  1       0.587372
    chr3    10000000        10010000        chr17   900000  910000  1       1.02558
    chr3    10000000        10010000        chr17   1030000 1040000 1       0.718195
    chr3    10000000        10010000        chr17   1320000 1330000 1       0.803212
    chr3    10000000        10010000        chr17   1500000 1510000 1       0.925146
    chr3    10000000        10010000        chr17   1750000 1760000 1       0.950326
    chr3    10000000        10010000        chr17   1800000 1810000 1       0.745982


Python API
----------

See: 

- Jupyter Notebook `API walkthrough <https://github.com/mirnylab/cooler-binder/blob/master/cooler_api.ipynb>`_.
- The :ref:`api-reference` for more information.

The :py:mod:`cooler` library provides a thin wrapper over the excellent NumPy-aware `h5py <http://docs.h5py.org/en/latest/>`_ Python interface to HDF5. It supports creation of cooler files and the following types of **range queries** on the data:

- Tabular selections are retrieved as Pandas DataFrames and Series.
- Matrix  selections are retrieved as NumPy arrays, DataFrames, or SciPy sparse matrices.
- Metadata is retrieved as a json-serializable Python dictionary.
- Range queries can be supplied using either integer bin indexes or genomic coordinate intervals.


::

    >>> import cooler
    >>> import matplotlib.pyplot as plt
    >>> c = cooler.Cooler('bigDataset.cool')
    >>> resolution = c.binsize
    >>> mat = c.matrix(balance=True).fetch('chr5:10,000,000-15,000,000')
    >>> plt.matshow(np.log10(mat), cmap='YlOrRd')

::

    >>> import multiprocessing as mp
    >>> import h5py
    >>> pool = mp.Pool(8)
    >>> c = cooler.Cooler('bigDataset.cool')
    >>> weights, stats = cooler.balance_cooler(c, map=pool.map, ignore_diags=3, min_nnz=10)


URI String
----------
The default location for a single-cooler .cool file is the root group ``/`` of the HDF5 file. It does not need to be explicitly specified.

.. code-block:: python

    >>> import cooler
    >>> c = cooler.Cooler('data/WT.DpnII.10kb.cool')
    >>> c = cooler.Cooler('data/WT.DpnII.10kb.cool::/')  # same as above

However, coolers can be stored at any level of the HDF5 hierarchy and qualified using a URI string of the form ``/path/to/cool/file::/path/to/cooler/group``.

.. code-block:: python
    
    >>> c1 = cooler.Cooler('data/WT.DpnII.multi.cool::10kb')
    >>> c2 = cooler.Cooler('data/WT.DpnII.multi.cool::1kb')


Data selection
--------------
Several :class:`cooler.Cooler` methods return data selectors. They don't retrieve any data from disk until queried. There are several ways to query using selectors. Genomic intervals can be provided using UCSC-style strings ``'{chrom}:{start}-{end}'`` or chrom-start-end triples ``(str, int, int)``. For regions with start and end that are not multiples of the resolution, selectors return the range of shortest range bins that fully contains the open interval [start, end).


Table selectors (chroms, bins, pixels)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- lazily select columns or lists of columns, returning new selectors
- query table rows using integer/slice indexing syntax
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



Matrix selector
~~~~~~~~~~~~~~~

- 2D bin range queries using slice indexing syntax
- 2D genomic range range queries using the ``fetch`` method


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


Dask
~~~~

Dask data structures provide a way to manipulate and distribute computations on larger-than-memory data using familiar APIs.
The sandboxed ``read_table`` function can be used to generate a dask dataframe backed by the pixel table of a Cooler as follows:

.. code-block:: python

    >>> from cooler.sandbox.dask import read_table
    >>> df = daskify(c.filename, 'pixels')

    >>> df
    Dask DataFrame Structure:
                    bin1_id bin2_id  count
    npartitions=223                       
    0                 int64   int64  int64
    9999999             ...     ...    ...
    ...                 ...     ...    ...
    2219999999          ...     ...    ...
    2220472929          ...     ...    ...
    Dask Name: daskify, 223 tasks

    >>> df = cooler.annotate(df, c.bins(), replace=False)
    >>> df
    Dask DataFrame Structure:
                    chrom1 start1   end1  weight1  chrom2 start2   end2  weight2 bin1_id bin2_id  count
    npartitions=31                                                                                     
    None            object  int64  int64  float64  object  int64  int64  float64   int64   int64  int64
    None               ...    ...    ...      ...     ...    ...    ...      ...     ...     ...    ...
    ...                ...    ...    ...      ...     ...    ...    ...      ...     ...     ...    ...
    None               ...    ...    ...      ...     ...    ...    ...      ...     ...     ...    ...
    None               ...    ...    ...      ...     ...    ...    ...      ...     ...     ...    ...
    Dask Name: getitem, 125 tasks

    >>> df = df[df.chrom1 == df.chrom2]
    >>> grouped = df.groupby(df.bin2_id - df.bin1_id)
    >>> x = grouped['count'].sum()
    >>> x
    Dask Series Structure:
    npartitions=1
    None    int64
    None      ...
    Name: count, dtype: int64
    Dask Name: series-groupby-sum-agg, 378 tasks

    >>> x.compute()
    0       476155231
    1       284724453
    2       139952477
    3        96520218
    4        71962080
    5        56085850
    6        45176881
    7        37274367
    8        31328555
    9        26781986
    10       23212616
    11       20366934
    12       18066135
    13       16159826
    14       14584058
    15       13249443
    16       12117854
    17       11149845
    ...

Learn more about the `Dask <https://dask.org/>`_ project.
