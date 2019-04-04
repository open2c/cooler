# Cooler

[![Build Status](https://travis-ci.org/mirnylab/cooler.svg?branch=master)](https://travis-ci.org/mirnylab/cooler)
[![Documentation Status](https://readthedocs.org/projects/cooler/badge/?version=latest)](http://cooler.readthedocs.org/en/latest/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/cooler/README.html)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mirnylab/cooler-binder/master)
[![Join the chat at https://gitter.im/mirnylab/cooler](https://badges.gitter.im/mirnylab/cooler.svg)](https://gitter.im/mirnylab/cooler?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![DOI](https://zenodo.org/badge/49553222.svg)](https://zenodo.org/badge/latestdoi/49553222)

## A cool place to store your Hi-C

Cooler is a support library for a **sparse, compressed, binary** persistent storage format, called _cooler_, used to store genomic interaction data, such as Hi-C contact matrices. 

The _cooler_ file format is a reference implementation of a genomic matrix data model using [HDF5](https://en.wikipedia.org/wiki/Hierarchical_Data_Format) as the container format.

The `cooler` package aims to provide the following functionality:

- Build contact matrices at any resolution from a [list of contacts](https://github.com/4dn-dcic/pairix).
- Query a contact matrix.
- Export and visualize the data.
- Perform efficient out-of-core operations, such as aggregation and contact matrix normalization (a.k.a. balancing).
- Facilitate working with potentially larger-than-memory data.

To get started:

- Read the [documentation](http://cooler.readthedocs.org/en/latest/).
- See the Jupyter Notebook [walkthrough](https://github.com/mirnylab/cooler-binder).
- _cool_ files from published Hi-C data sets are available at `ftp://cooler.csail.mit.edu/coolers`.
- Many more multires (_mcool_) files are available on the [4DN data portal](https://data.4dnucleome.org/visualization/index).

Related projects:

- Process Hi-C data with [distiller](https://github.com/mirnylab/distiller).
- Downstream analysis with [cooltools](https://github.com/mirnylab/cooltools) (WIP).
- Visualize your Cooler data with [HiGlass](http://higlass.io)!


### Installation

Requirements:

- Python 2.7/3.4+
- libhdf5 and Python packages `numpy`, `scipy`, `pandas`, `h5py`. We highly recommend using the `conda` package manager to install scientific packages like these. To get it, you can either install the full [Anaconda](https://www.continuum.io/downloads) Python distribution or just the standalone [conda](http://conda.pydata.org/miniconda.html) package manager.

Install from PyPI using pip.
```sh
$ pip install cooler
```

If you are using `conda`, you can alternatively install `cooler` from the [bioconda](https://bioconda.github.io/index.html) channel.
```sh
$ conda install -c conda-forge -c bioconda cooler
```

See the [docs](http://cooler.readthedocs.org/en/latest/) for more information.


### Command line interface

The `cooler` package includes command line tools for creating, querying and manipulating cooler files.

```bash
$ cooler cload pairs hg19.chrom.sizes:10000 $PAIRS_FILE out.10000.cool
$ cooler balance -p 10 out.10000.cool
$ cooler dump -b -t pixels --header --join -r chr3:10M-12M -r2 chr17 out.10000.cool | head
```

```
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
```

See also:

- [CLI Reference](http://cooler.readthedocs.io/en/latest/cli.html).
- Jupyter Notebook [walkthrough](https://nbviewer.jupyter.org/github/mirnylab/cooler-binder/blob/master/cooler_cli.ipynb).

### Python API

The `cooler` library provides a thin wrapper over the excellent [h5py](http://docs.h5py.org/en/latest/) Python interface to HDF5. It supports creation of cooler files and the following types of **range queries** on the data:

- Tabular selections are retrieved as Pandas DataFrames and Series.
- Matrix  selections are retrieved as NumPy arrays or SciPy sparse matrices.
- Metadata is retrieved as a json-serializable Python dictionary.
- Range queries can be supplied using either integer bin indexes or genomic coordinate intervals. Note that queries with coordinate intervals that are not multiples of the bin size will return the range of shortest range bins that fully contains the open interval [start, end).

```python

>>> import cooler
>>> import matplotlib.pyplot as plt
>>> c = cooler.Cooler('bigDataset.cool')
>>> resolution = c.binsize
>>> mat = c.matrix(balance=True).fetch('chr5:10,000,000-15,000,000')
>>> plt.matshow(np.log10(mat), cmap='YlOrRd')
```

```python
>>> import multiprocessing as mp
>>> import h5py
>>> pool = mp.Pool(8)
>>> c = cooler.Cooler('bigDataset.cool')
>>> weights, stats = cooler.balance_cooler(c, map=pool.map, ignore_diags=3, min_nnz=10)
```

See also:

- [API Reference](http://cooler.readthedocs.io/en/latest/api.html).
- Jupyter Notebook [walkthrough](https://nbviewer.jupyter.org/github/mirnylab/cooler-binder/blob/master/cooler_api.ipynb).

### Schema

The _cool_ format implements a simple [data model](http://cooler.readthedocs.io/en/latest/schema.html) that stores a genomic matrix in a sparse representation, crucial for developing robust tools for use on increasingly high resolution Hi-C data sets, including streaming and [out-of-core](https://en.wikipedia.org/wiki/Out-of-core_algorithm) algorithms.

The data tables in a cooler file are stored in a **columnar** representation as HDF5 groups of 1D array datasets of equal length. A symmetric contact matrix is represented as a single table containing only the **nonzero upper triangle** pixels.


### Contributing

[Pull requests](https://akrabat.com/the-beginners-guide-to-contributing-to-a-github-project/) are welcome. The current requirements for testing are `pytest` and `mock`.

For development, clone and install in "editable" (i.e. development) mode with the `-e` option. This way you can also pull changes on the fly.
```sh
$ git clone https://github.com/mirnylab/cooler.git
$ cd cooler
$ pip install -e .
```

### License

BSD (New)


### Citation

Abdennur, N., and Mirny, L. (2019). Cooler: scalable storage for Hi-C data and other genomically-labeled arrays. _BioRxiv_, 557660. doi: [10.1101/557660](https://doi.org/10.1101/557660).
