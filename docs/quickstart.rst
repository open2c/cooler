Quickstart
==========


Installation
------------

Requirements:

- Python 2.7 or 3.4 and higher
- libhdf5
- Python packages ``numpy``, ``scipy``, ``pandas``, ``h5py``. 

We highly recommend using the conda package manager to install scientific packages like these. To get ``conda``, you can download either the full `Anaconda <https://www.continuum.io/downloads>`_ Python distribution which comes with lots of data science software or the minimal `Miniconda <http://conda.pydata.org/miniconda.html>`_ distribution which is just the standalone package manager plus Python. In the latter case, you can install the packages as follows:

::

    $ conda install numpy scipy pandas h5py


Install ``cooler`` from PyPI using pip.

::

    $ pip install cooler

All other Python package dependencies are automatically handled by pip.

.. Additionally, the following tools are required for building ``cool`` files from contact lists:
.. - Parallel gzip ``pigz``. Install using your system package manager.
.. - Tabix/bgzf. These come with `Samtools <http://www.htslib.org/download/>`_ but are also available on system package managers like ``brew`` (Mac OS) and ``apt`` (Ubuntu). Alternatively, if you are using ``conda``, consider adding the `bioconda <https://bioconda.github.io/>`_ channel to get access to many more bioinformatics packages.


Command line interface
----------------------

See:

- Jupyter Notebook `CLI walkthrough <https://github.com/mirnylab/cooler-binder/blob/master/cooler_cli.ipynb>`_.
- The `CLI Reference <http://cooler.readthedocs.io/en/latest/cli.html>`_ for more information.


The ``cooler`` library includes utilities for creating, querying, merging and manipulating .cool files and for performing out-of-core matrix balancing on a contact matrix of any resolution.

::

    $ cooler makebins $CHROMSIZES_FILE $BINSIZE > bins.10kb.bed
    $ cooler cload bins.10kb.bed $CONTACTS_FILE out.cool
    $ cooler balance -p 10 out.cool
    $ cooler dump -b -t pixels --header --join -r chr3:10,000,000-12,000,000 -r2 chr17 out.cool | head

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

The ``cooler`` library provides a thin wrapper over the excellent NumPy-aware `h5py <http://docs.h5py.org/en/latest/>`_ Python interface to HDF5. It supports creation of cooler files and the following types of **range queries** on the data:

- Tabular selections are retrieved as Pandas DataFrames and Series.
- Matrix  selections are retrieved as NumPy arrays, DataFrames, or SciPy sparse matrices.
- Metadata is retrieved as a json-serializable Python dictionary.
- Range queries can be supplied using either integer bin indexes or genomic coordinate intervals.


::

    >>>  import multiprocessing as mp
    >>>  import h5py
    >>>  pool = mp.Pool(8)
    >>>  f = h5py.File('bigDataset.cool', 'r')
    >>>  weights, stats = cooler.ice.iterative_correction(f, map=pool.map, ignore_diags=3, min_nnz=10)

::

    >>> import cooler
    >>> import matplotlib.pyplot as plt
    >>> c = cooler.Cooler('bigDataset.cool')
    >>> resolution = c.info['bin-size']
    >>> mat = c.matrix(balance=True).fetch('chr5:10,000,000-15,000,000')
    >>> plt.matshow(np.log10(mat), cmap='YlOrRd')
