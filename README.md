# Cooler

[![Build Status](https://travis-ci.org/mirnylab/cooler.svg?branch=master)](https://travis-ci.org/mirnylab/cooler)
[![Documentation Status](https://readthedocs.org/projects/cooler/badge/?version=latest)](http://cooler.readthedocs.org/en/latest/)

## A cool place to store your Hi-C

Cooler is a **sparse, compressed, binary** persistent storage format for Hi-C contact maps based on [HDF5](https://www.hdfgroup.org/HDF5/).

See [example Jupyter notebook](https://github.com/mirnylab/cooler-binder/blob/master/cooler_quickstart.ipynb).

Some published data sets available at `ftp://cooler.csail.mit.edu/coolers`.

The cooler [format](#schema) implements a simple schema and data model that stores a high resolution contact matrix along with important auxiliary data such as scaffold information, genomic bin annotations, and basic metadata. Data tables are stored in a **columnar** representation as HDF5 Groups of 1D array datasets of equal length. The contact matrix itself is stored as a single table containing only the **nonzero upper triangle** pixels.

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


### <a id="schema"></a>Schema

The tables and indexes can be represented in the [Datashape](http://datashape.readthedocs.org/en/latest/) layout language:
```
{
  chroms: {
    name:     typevar['Nchroms'] * string['ascii'],
    length:   typevar['Nchroms'] * int32,
  },
  bins: {
    chrom:    typevar['Nbins'] * categorical[typevar['name'], type=string, ordered=True],
    start:    typevar['Nbins'] * int32,
    end:      typevar['Nbins'] * int32,
    weight:   typevar['Nbins'] * float64
  },
  pixels: {
    bin1_id:  typevar['Nnz'] * int64,
    bin2_id:  typevar['Nnz'] * int64,
    count:    typevar['Nnz'] * int32
  },
  indexes: {
    chrom_offset:  typevar['Nchroms'] * int64,
  	bin1_offset:   typevar['Nbins'] * int64
  }
}
```

Attributes (metadata):
```
nchroms         : <int> Number of rows in scaffolds table
nbins           : <int> Number of rows in bins table
nnz             : <int> Number of rows in matrix table
bin-type        : {"fixed" or "variable"}
bin-size        : <int or null> Size of bins in base pairs if bin-type is "fixed"
genome-assembly : <string> Name of genome assembly
library-version : <string> Version of cooler library that created the file
format-version  : <string> The version of the current format
format-url      : <url> URL to page providing format details
creation-date   : <datetime> Date the file was built
metadata        : <json> custom metadata about the experiment
```

Matrix storage format:

- The `bins` table is lexicographically sorted by `chrom`, then `start`, respecting the order of the contigs in the `chroms` table.
- The `pixels` table is lexicographically sorted by `bin1_id`, then `bin2_id`.
- Offset pointers are used to facilitate matrix queries. This is effectively a [compressed sparse row](https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_.28CSR.2C_CRS_or_Yale_format.29) storage scheme for a symmetric matrix.


Notes:

- Any number of additional optional columns can be added to each table. (e.g. normalization vectors, quality masks).
- Genomic coordinates are assumed to be 0-based and intervals half-open (1-based ends).
