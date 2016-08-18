# Cooler

[![Build Status](https://travis-ci.org/mirnylab/cooler.svg?branch=master)](https://travis-ci.org/mirnylab/cooler)
[![Documentation Status](https://readthedocs.org/projects/cooler/badge/?version=latest)](http://cooler.readthedocs.org/en/latest/)

## A cool place to store your Hi-C

Cooler is a **sparse, compressed, binary** persistent storage format for Hi-C contact maps based on [HDF5](https://en.wikipedia.org/wiki/Hierarchical_Data_Format).

- See example [Jupyter notebook](https://github.com/mirnylab/cooler-binder/blob/master/cooler_quickstart.ipynb) or launch the [Binder](https://github.com/mirnylab/cooler-binder).
- Read the [docs](http://cooler.readthedocs.org/en/latest/).
- Some published data sets are available at `ftp://cooler.csail.mit.edu/coolers`.

The cooler [format](http://cooler.readthedocs.io/en/latest/intro.html#data-model) implements a simple schema and data model that stores a high resolution contact matrix along with important auxiliary data such as scaffold information, genomic bin annotations, and basic metadata. Data tables are stored in a **columnar** representation as HDF5 Groups of 1D array datasets of equal length. The contact matrix itself is stored as a single table containing only the **nonzero upper triangle** pixels.

The `cooler` [library](https://github.com/mirnylab/cooler) provides a thin wrapper over the excellent [h5py](http://docs.h5py.org/en/latest/) Python interface to HDF5. It supports creation of cooler files and the following types of **range queries** on the data:

- Tablular selections are retrieved as Pandas DataFrames and Series.
- Matrix slice selections are retrieved as SciPy sparse matrices.
- Metadata is retrieved as a json-serializable Python dictionary.
- Range queries can be supplied using either integer bin indexes or genomic coordinate intervals.

The `cooler` library also includes utilities for performing contact **matrix balancing** on a cooler file of any resolution.


### Installation

Requirements:

- Python 2.7/3.3+
- libhdf5 and Python packages `numpy`, `scipy`, `pandas`, `h5py`. If you don't have them installed already, we recommend you use the [conda](http://conda.pydata.org/miniconda.html) package manager to manage these dependencies instead of pip.

Install from PyPI using pip.
```sh
$ pip install cooler
```

For development, clone and install in "editable" (i.e. development) mode with the `-e` option.
```sh
$ git clone https://github.com/mirnylab/cooler.git
$ cd cooler
$ pip install -e .
```

