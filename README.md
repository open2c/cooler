# Cooler

[![Build Status](https://travis-ci.org/mirnylab/cooler.svg?branch=master)](https://travis-ci.org/mirnylab/cooler)

## A cool place to store your Hi-C

Cooler is a **sparse, compressed, binary** persistent storage format for Hi-C contact maps based on HDF5.

The `cooler` library implements a simple **schema** to store a high resolution contact matrix along with important auxiliary data such as scaffold information, genomic bin annotations, and basic metadata.

Data tables are stored in a **columnar** representation as groups of 1D HDF5 array datasets of the same length. The contact matrix itself is stored as a table containing only the _nonzero upper triangle entries_.

The library API provides a thin Python wrapper over [h5py](http://docs.h5py.org/en/latest/) for **range queries** on the data: 
- Table sections are retrieved as Pandas `DataFrame`s
- Matrix slices are retrieved as SciPy sparse matrices or NumPy `ndarray`s
- The metadata is retrieved as a dictionary.

Rather than build on top of a more full-featured, opinionated library like PyTables (or `pandas.HDFStore` built on top of that), we provide a simple and transparent data layout on top of HDF5 that supports random access range queries and can be easily [migrated](https://github.com/blaze/odo).


### Schema

Required attributes (metadata):
```
id              : <string or null> Name or id for a file
bin-type        : {"fixed" or "variable"}
bin-size        : <int or null> Size of bins in bp if bin-type is "fixed"
format-url      : <url> Static URL to page providing format details
format-version  : <tuple> The version of the current format
creation-date   : <datetime> Date the file was built
genome-assembly : <string> Name of genome assembly
nchroms         : <int> Number of rows in scaffolds table
nbins			: <int> Number of rows in bins table
nnz				: <int> Number of rows in matrix table
```

The required tables and indexes can be represented in the [Datashape](http://datashape.readthedocs.org/en/latest/) layout language:
```
{
  scaffolds: {
    length:   typevar['Nchroms'] * int64, 
    name:     typevar['Nchroms'] * string[32, 'A']
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
- The `bins` table is lexicographically sorted by `chrom_id`, `start`, `end`.
- The `matrix` table is lexicographically sorted by `bin1_id`, then `bin2_id`.
- Simple offset pointer indexes are used to speed up matrix queries.


