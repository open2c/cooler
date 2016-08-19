# Cooler

[![Build Status](https://travis-ci.org/mirnylab/cooler.svg?branch=master)](https://travis-ci.org/mirnylab/cooler)
[![Documentation Status](https://readthedocs.org/projects/cooler/badge/?version=latest)](http://cooler.readthedocs.org/en/latest/)

## A cool place to store your Hi-C

Cooler is a **sparse, compressed, binary** persistent storage format for Hi-C contact maps based on [HDF5](https://en.wikipedia.org/wiki/Hierarchical_Data_Format).

- Documentation is available [here](http://cooler.readthedocs.org/en/latest/).
- See example [Jupyter notebook](https://github.com/mirnylab/cooler-binder/blob/master/cooler_quickstart.ipynb) or [try it live](http://mybinder.org/repo/mirnylab/cooler-binder).
- Some published data sets are available at `ftp://cooler.csail.mit.edu/coolers`.

As published Hi-C datasets increase in sequencing depth and resolution, a simple sparse representation lends itself better not only to storage but also to streaming and [out-of-core](https://en.wikipedia.org/wiki/Out-of-core_algorithm) algorithms for analysis. The cooler [format](http://cooler.readthedocs.io/en/latest/intro.html#data-model) implements a simple schema and data model that stores a high resolution contact matrix in a sparse representation along with important auxiliary data such as scaffold information, genomic bin annotations, and basic metadata. Data tables are stored in a **columnar** representation as HDF5 Groups of 1D array datasets of equal length. The contact matrix itself is stored as a single table containing only the **nonzero upper triangle** pixels.

The `cooler` [library](https://github.com/mirnylab/cooler) provides a thin wrapper over the excellent [h5py](http://docs.h5py.org/en/latest/) Python interface to HDF5. It supports creation of cooler files and the following types of **range queries** on the data:

- Tabular selections are retrieved as Pandas DataFrames and Series.
- Matrix  selections are retrieved as SciPy sparse matrices.
- Metadata is retrieved as a json-serializable Python dictionary.
- Range queries can be supplied using either integer bin indexes or genomic coordinate intervals.


```python

>>>  import cooler
>>>  import matplotlib.pyplot as plt
>>>  c = cooler.Cooler('bigDataset.cool')
>>>  resolution = c.info['bin-size']
>>>  mat = c.matrix(balance=True).fetch('chr5:10,000,000-15,000,000')
>>>  plt.matshow(np.log10(mat.toarray()), cmap='YlOrRd')
```

The `cooler` library also includes utilities for performing out-of-core contact **matrix balancing** on a cooler file of any resolution. See the [docs](http://cooler.readthedocs.org/en/latest/) for more information.


### Installation

Requirements:

- Python 2.7/3.3+
- libhdf5 and Python packages `numpy`, `scipy`, `pandas`, `h5py`. If you don't have them installed already, we recommend you use the [conda](http://conda.pydata.org/miniconda.html) package manager to manage these dependencies instead of pip.

Install from PyPI using pip.
```sh
$ pip install cooler
```

For the latest, unstable version, clone and install from master or install directly from the repo.
```sh
$ pip install git+git://github.com/mirnylab/cooler.git
```

For development, clone and install in "editable" (i.e. development) mode with the `-e` option. This way you can also pull changes on the fly.
```sh
$ git clone https://github.com/mirnylab/cooler.git
$ cd cooler
$ pip install -e .
```

### Contributing

[Pull requests](https://akrabat.com/the-beginners-guide-to-contributing-to-a-github-project/) are welcome. The current requirements for testing are `nose` and `mock`.
