# Cooler

[![Build Status](https://travis-ci.org/mirnylab/cooler.svg?branch=master)](https://travis-ci.org/mirnylab/cooler)
[![Documentation Status](https://readthedocs.org/projects/cooler/badge/?version=latest)](http://cooler.readthedocs.org/en/latest/)

## A cool place to store your Hi-C

Cooler is a **sparse, compressed, binary** persistent storage format for Hi-C contact maps based on HDF5.

See [example Jupyter notebook](https://gist.github.com/nvictus/904160bca9d0e8d5aeeb).

The `cooler` library implements a simple **schema** to store a high resolution contact matrix along with important auxiliary data such as scaffold information, genomic bin annotations, and basic metadata.

Data tables are stored in a **columnar** representation as groups of 1D HDF5 array datasets of the same length. The contact matrix itself is stored as a table containing only the **nonzero upper triangle** pixels.

The library API provides a thin Python wrapper over [h5py](http://docs.h5py.org/en/latest/) for **range queries** on the data:
- Table selections are retrieved as Pandas `DataFrame`s
- Matrix slice selections are retrieved as SciPy sparse matrices or NumPy `ndarray`s
- The metadata is retrieved as a json-serializable dictionary.

Range queries can be supplied as either integer bin indexes or genomic coordinate intervals.

### Installation

Requirements:

- Python 2.7/3.3+
- libhdf5 and Python packages `numpy`, `scipy`, `pandas`, `h5py`. If you don't have them installed already, we recommend you use the [conda](http://conda.pydata.org/miniconda.html) package manager to manage these dependencies instead of pip.

Clone or download the source archive and install using pip, or pass the repo url directly to pip.
```sh
$ pip install git+https://github.com/mirnylab/cooler.git@v0.2
```

To install in "editable" (i.e. development) mode use the `-e` option.
```sh
$ git clone https://github.com/mirnylab/cooler.git
$ cd cooler
$ pip install -e .
```


### Schema

Required attributes (metadata):
```
genome-assembly : <string> Name of genome assembly
bin-type        : {"fixed" or "variable"}
bin-size        : <int or null> Size of bins in bp if bin-type is "fixed"
nchroms         : <int> Number of rows in scaffolds table
nbins           : <int> Number of rows in bins table
nnz             : <int> Number of rows in matrix table
format-url      : <url> URL to page providing format details
format-version  : <string> The version of the current format
generated-by    : <string> Agent that created the file
creation-date   : <datetime> Date the file was built
metadata        : <json> custom metadata about the experiment
```

The required tables and indexes can be represented in the [Datashape](http://datashape.readthedocs.org/en/latest/) layout language:
```
{
  scaffolds: {
    name:     typevar['Nchroms'] * string[32, 'A']
    length:   typevar['Nchroms'] * int64,
  },
  bins: {
    chrom_id: typevar['Nbins'] * int32,
    start:    typevar['Nbins'] * int64,
    end:      typevar['Nbins'] * int64
  },
  matrix: {
    bin1_id:  typevar['Nnz'] * int32,
    bin2_id:  typevar['Nnz'] * int32,
    count:    typevar['Nnz'] * int32
  },
  indexes: {
    chrom_offset: typevar['Nchroms'] * int32,
  	bin1_offset:   typevar['Nbins'] * int32
  }
}
```

Notes:
- Any number of additional optional columns can be added to each table. (e.g. quality masks, normalization vectors).
- Genomic coordinates are assumed to be 0-based and intervals half-open (1-based ends).

Matrix storage format:
- The `bins` table is lexicographically sorted by `chrom_id`, `start`, `end`.
- The `matrix` table is lexicographically sorted by `bin1_id`, then `bin2_id`.
- Offset pointers are used to facilitate matrix queries. This is effectively a [compressed sparse row](https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_.28CSR.2C_CRS_or_Yale_format.29) storage scheme for a symmetric matrix.

Rather than build on top of a more full-featured, opinionated library like PyTables (or `pandas.HDFStore` built on top of that), we provide a simple and transparent data layout on top of HDF5 that supports random access range queries and can be easily [migrated](https://github.com/blaze/odo).

See also:
- [hdf2tab](https://github.com/blajoie/hdf2tab) converts dense Hi-C matrices stored in HDF5 files to tabular text files.
- The [biom](https://github.com/biocore/biom-format) format is an HDF5-based format for metagenomic observation matrices.

