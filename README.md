# Cooler

[![Build Status](https://travis-ci.org/mirnylab/cooler.svg?branch=master)](https://travis-ci.org/mirnylab/cooler)
[![Documentation Status](https://readthedocs.org/projects/cooler/badge/?version=latest)](http://cooler.readthedocs.org/en/latest/)
[![Binder](http://mybinder.org/badge.svg)](http://mybinder.org:/repo/mirnylab/cooler-binder)

## A cool place to store your Hi-C

Cooler is a support library for a **sparse, compressed, binary** persistent storage format for Hi-C contact matrices, called `cool`, which is based on [HDF5](https://en.wikipedia.org/wiki/Hierarchical_Data_Format).

Cooler aims to provide the following functionality:

- Generate contact matrices from contact lists at arbitrary resolutions.
- Store contact matrices efficiently in COOL format based on widely used HDF5 container.
- Perform out-of-core genome wide contact matrix normalization (a.k.a. balancing)
- Perform fast range queries on a contact matrix.
- Convert contact matrices between formats.
- Provide a clean and well-documented Python API to work with Hi-C data.


To get started:

- Documentation is available [here](http://cooler.readthedocs.org/en/latest/).
- Walkthrough with a [Jupyter notebook](https://github.com/mirnylab/cooler-binder).
- Some published data sets are available at `ftp://cooler.csail.mit.edu/coolers`.


### Installation

Requirements:

- Python 2.7/3.3+
- libhdf5 and Python packages `numpy`, `scipy`, `pandas`, `h5py`, `click`. If you don't have them installed already, we recommend you use the [conda](http://conda.pydata.org/miniconda.html) package manager to manage these dependencies instead of pip.
- See the [docs](http://cooler.readthedocs.org/en/latest/) for more information.

Install from PyPI using pip.
```sh
$ pip install cooler
```


### Command line interface

The `cooler` library includes utilities for performing out-of-core contact **matrix balancing** on a cooler file of any resolution. See the [docs](http://cooler.readthedocs.org/en/latest/) for more information.

```bash
$ cooler binnifiy $chromSizesFile $binsizes
$ cooler cload $bedfile $contactlist $outputcooler
$ cooler balance $coolerfile
$ cooler dump $outfilename
```

### Python API

The `cooler` [library](https://github.com/mirnylab/cooler) provides a thin wrapper over the excellent [h5py](http://docs.h5py.org/en/latest/) Python interface to HDF5. It supports creation of cooler files and the following types of **range queries** on the data:

- Tabular selections are retrieved as Pandas DataFrames and Series.
- Matrix  selections are retrieved as SciPy sparse matrices.
- Metadata is retrieved as a json-serializable Python dictionary.
- Range queries can be supplied using either integer bin indexes or genomic coordinate intervals.

```python

>>> import cooler
>>> import matplotlib.pyplot as plt
>>> c = cooler.Cooler('bigDataset.cool')
>>> resolution = c.info['bin-size']
>>> mat = c.matrix(balance=True).fetch('chr5:10,000,000-15,000,000')
>>> plt.matshow(np.log10(mat.toarray()), cmap='YlOrRd')
```

Also see the [Jupyter notebook](https://github.com/mirnylab/cooler-binder) walkthrough.


### Cooler Schema

The `cool` [format](http://cooler.readthedocs.io/en/latest/intro.html#data-model) implements a simple schema that stores a contact matrix in a sparse representation, crucial for developing robust tools for use on increasingly high resolution Hi-C data sets, including streaming and [out-of-core](https://en.wikipedia.org/wiki/Out-of-core_algorithm) algorithms.

The data tables in a `cool` file are stored in a **columnar** representation as HDF5 Groups of 1D array datasets of equal length. The contact matrix itself is stored as a single table containing only the **nonzero upper triangle** pixels.


### Contributing

[Pull requests](https://akrabat.com/the-beginners-guide-to-contributing-to-a-github-project/) are welcome. The current requirements for testing are `nose` and `mock`.

For development, clone and install in "editable" (i.e. development) mode with the `-e` option. This way you can also pull changes on the fly.
```sh
$ git clone https://github.com/mirnylab/cooler.git
$ cd cooler
$ pip install -e .
```




