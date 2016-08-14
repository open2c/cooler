Cooler
======

Cooler is a sparse data model, schema and HDF5-based file format for high resolution Hi-C contact maps. 

The cooler format implements a simple data model that stores a high resolution contact matrix along with important auxiliary data such as scaffold information, genomic bin annotations, and basic metadata. 


Data model
----------

We model a contact matrix using three tables.

chroms
~~~~~~

+ Required columns: *name*
+ Order: *enumeration*

An enumeration of the chromosomes, scaffolds or contigs of the assembly that the data is mapped to. This information can be extracted from the bin table below, but is included separately for convenience. This enumeration is the intended ordering of the chromosomes as they would appear in a global contact matrix. Additional columns can provide metadata on the chromosomes, such as their length.

bins
~~~~

+ Required columns: *chrom, start, end*
+ Order: *chrom [enum], start*

An enumeration of the concatenated genomic bins that make up a single dimension or axis of the global contact matrix. Genomic bins can be of fixed size or variable sizes (e.g. restriction fragments). A genomic bin is defined by the triple (chrom, start, end), where start is zero-based and end is 1-based. The order is significant: the bins are sorted by chrom (based on the chromosome enumeration) then by start, and each genomic bin is implicitly endowed with a 0-based bin ID from this ordering (i.e., its row number in the table). Additional columns can be added to describe other bin-associated properties such as normalization vectors and bin-level masks.

pixels
~~~~~~

+ Required columns: *bin1ID, bin2ID, count*
+ Order: *bin1ID, bin2ID*

The contact matrix is stored as a single table containing only the nonzero upper triangle elements, assuming the ordering of the bins given by the bin table. Each row defines a non-zero element of the contact matrix. Additional columns can be appended to store pixel-associated properties such as pixel-level masks or filtered and transformed versions of the data. Currently, the pixels are sorted lexicographically by the bin ID of the 1st axis (matrix row) then the bin ID of the 2nd axis (matrix column).


Why model it this way?
Balances the tradeoff between simplicity, terseness and flexibility. Flatter is better than nested.
It also models a text serialization that doesn’t require splitting a contact matrix into separate files for each contig-contig submatrix.
Separating bins (annotation of the axes labels) from pixels (the matrix data) allows for easy inclusion of bin-level properties without introducing redundancy
It is ideal for many types of serial and out-of-core processing of very large contact matrices, such as matrix balancing


Schema
------

See :ref:`current-version`.

Column-oriented vs record-oriented tables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Why is the reference schema column-oriented?

- Cheap column addition/removal
- Better compression ratios
- Blazingly fast I/O speed can be achieved with new compressors such as blosc
- Migrating to other column stores such as bcolz, Apache Parquet.
- The schema is fully interchangeable with a record-oriented representation (e.g., traditional SQL databases).
- Tradeoff between flexibility and number of read cycles to fetch all columns of a table.


Supporting a matrix “view”
~~~~~~~~~~~~~~~~~~~~~~~~~~

Indexes are stored as 1D datasets in a separate group. The current indexes can be thought of as run-length encodings of the bins/chrom and pixels/bin1_id columns, respectively.


Limitations
~~~~~~~~~~~

A full rectangular matrix “view” of the data must be modeled on top of this representation.
2D range queries must be computed using indexes, and the sort order on the pixels and types of usable indexing strategies are strongly related.


Container
---------

The reference cooler implementation uses HDF5 as a container format, which supports chunking, compression, and random access. It can be accessed from virtually any programming environment ...

Data tables are stored in a columnar representation as HDF5 Groups of 1D array datasets of equal length.


Library
-------


The excellent `h5py <http://docs.h5py.org/en/latest/>`_ Python interface to HDF5 provides direct access to the group and dataset structure of a cooler file. h5py translates HDF5 dataset queries directly into NumPy arrays.

The cooler library provides an additional thin wrapper over h5py to support creation and conversion of cooler files as well as both tabular and sparse matrix views on the data. Range queries can be made using either integer bin indexing or genomic interval strings.

Table range queries are retrieved as Pandas DataFrames and Series.
Matrix range queries are retrieved as SciPy sparse matrices.
Metadata is retrieved as a json-compatible Python dictionary.

The cooler library also includes utilities for performing contact matrix balancing on a cooler file of any resolution.

Try it out in a Jupyter notebook using `Binder <https://github.com/mirnylab/cooler-binder>`_.


HDF5 bindings in other languages
--------------------------------


- canonical C-library libhdf5
- rhdf5
- perl
- Java
- scala
- MATLAB



Glossary
--------

HDF5 is a general purpose binary container format for large scientific datasets.

h5py is a Python library providing low-level bindings to the libhdf5 C-library and a high-level, numpy-aware API to interact with HDF5 files on disk.

Cooler is a flexible binary schema for Hi-C data based on a two-table sparse data model.

Cooler [Cooler-HDF5?] is also the name of an implementation of the Cooler schema in HDF5.

Cooler [pycooler?] is a Python package providing an API to create cooler-hdf5 files and to interact with them both as data frames and sparse matrices.
