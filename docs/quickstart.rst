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
