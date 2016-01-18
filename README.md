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

