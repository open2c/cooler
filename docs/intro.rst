Cooler
======

Cooler is a sparse data model, schema and HDF5-based file format for high resolution Hi-C contact maps. 

The cooler format implements a simple data model that stores a high resolution contact matrix along with important auxiliary data such as scaffold information, genomic bin annotations, and basic metadata.

Why? As published Hi-C datasets increase in sequencing depth and resolution, a simple sparse representation lends itself better not only to the increasing storage demands but also to performing streaming and `out-of-core <https://en.wikipedia.org/wiki/Out-of-core_algorithm>`_ algorithms for analysis.


Data model
----------

We model a contact matrix using three tables.

chroms
~~~~~~

+ Required columns: ``name``
+ Order: *enumeration*

An enumeration of the chromosomes, scaffolds or contigs of the assembly that the data is mapped to. This information can be extracted from the bin table below, but is included separately for convenience. This enumeration is the intended ordering of the chromosomes as they would appear in a global contact matrix. Additional columns can provide metadata on the chromosomes, such as their length.

bins
~~~~

+ Required columns: ``chrom, start, end [, weight]``
+ Order: ``chrom`` (*enum*), ``start``

An enumeration of the concatenated genomic bins that make up a single dimension or axis of the global contact matrix. Genomic bins can be of fixed size or variable sizes (e.g. restriction fragments). A genomic bin is defined by the triple (chrom, start, end), where start is zero-based and end is 1-based. The order is significant: the bins are sorted by chromosome (based on the chromosome enumeration) then by start, and each genomic bin is implicitly endowed with a 0-based bin ID from this ordering (i.e., its row number in the table). A reserved but optional column called ``weight`` can store weights for normalization or matrix balancing. Additional columns can be added to describe other bin-associated properties such as additional normalization vectors or bin-level masks.

pixels
~~~~~~

+ Required columns: ``bin1ID, bin2ID, count``
+ Order: ``bin1ID, bin2ID``

The contact matrix is stored as a single table containing only the nonzero upper triangle elements, assuming the ordering of the bins given by the bin table. Each row defines a non-zero element of the contact matrix. Additional columns can be appended to store pixel-associated properties such as pixel-level masks or filtered and transformed versions of the data. Currently, the pixels are sorted lexicographically by the bin ID of the 1st axis (matrix row) then the bin ID of the 2nd axis (matrix column).


Why model it this way?

To balance the tradeoff between simplicity, terseness and flexibility in an attempt to stay `Zen <https://www.python.org/dev/peps/pep-0020/>`_. 

+ The schema is flexible enough to describe a whole genome contact matrix, or any subset of a contact matrix, including single contig-contig tiles.
+ Given the variety of ways we might want to read the data or add new columns, flatter is better than nested.
+ For one, it makes the data much easier to stream and process in chunks, which ideal for many types of out-of-core algorithms on very large contact matrices.
+ Separating bins (annotations of the axis labels) from pixels (the matrix data) allows for easy inclusion of bin-level properties without introducing redundancy.

+ The same schema [bin + pixel table combination] also defines a plain text format, a simple serialization of the binary format, that doesn’t require splitting a contact matrix into a separate file for every contig-contig submatrix. One of two forms are possible:
    - Two-file: The bin table and pixel table are stored as separate tab-delimited files (BED file + sparse triple file). See the output format from `Hi-C Pro <http://nservant.github.io/HiC-Pro/RESULTS.html#intra-and-inter-chromosomal-contact-maps>`_.
    - Single-file ("tidy"): The ``bin1_id`` and ``bin2_id`` columns of the pixel table are replaced with annotations from the bin table, suffixed with `1` or `2` accordingly. There is an increased storage cost from the redudundacy of using a single file. The result is a `BEDPE <http://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format>`_-like file.


Schema
------

See :ref:`current-version`.

Column-oriented vs record-oriented tables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Why is the reference schema column-oriented?

- Cheap column addition/removal.
- Better compression ratios.
- Blazingly fast I/O speed can be achieved with new compressors such as `blosc <http://www.blosc.org/>`_.
- Easy to migrate to other column stores such as `bcolz <https://github.com/Blosc/bcolz>`_, Apache `Parquet <https://parquet.apache.org/>`_, and Apache `Arrow <http://blog.cloudera.com/blog/2016/02/introducing-apache-arrow-a-fast-interoperable-in-memory-columnar-data-structure-standard/>`_.

There is a tradeoff between flexibility and number of read cycles required to fetch all columns of a table, however, a column-oriented schema is fully interchangeable with a record-oriented representation (e.g., traditional SQL databases, CSV files).


Supporting a matrix “view”
~~~~~~~~~~~~~~~~~~~~~~~~~~

Indexes are stored as 1D datasets in a separate group. The current indexes can be thought of as run-length encodings of the ``bins/chrom`` and ``pixels/bin1_id`` columns, respectively.


Limitations
~~~~~~~~~~~

A complete rectangular matrix “view” of the data must be modeled on top of this representation. 2D range queries must be computed with the help of indexes. The sort order on the pixels and types of indexing strategies that can be used are strongly related. This could be changed in future versions of the schema.


Container
---------

The reference cooler implementation uses `HDF5 <https://www.hdfgroup.org/HDF5/>`_ as a container format. HDF5 is a hierarchical data format for homongenenously typed multidimensional arrays, which supports chunking, compression, and random access. The HDF5 file specification and open source standard library is maintained by the nonprofit HDF Group.

HDF5 files consist of three fundamental entities: groups, datasets, and attibutes. The hierarchical organization of an HDF5 file is conceptually analogous to a file system: *groups* are akin to directories and *datasets* (arrays) are akin to files. Additionally, key-value metadata can be attached to groups and datasets using *attributes*. The standard library provides the ability to access and manipulate these entities. There are bindings for virtually every platform and programming environment. To learn more in detail about HDF5, I recommend the book `HDF5 and Python <https://www.safaribooksonline.com/library/view/python-and-hdf5/9781491944981/ch01.html>`_ by Andrew Collette, the author of ``h5py``.

To implement the Cooler data model in HDF5, data tables are stored in a columnar representation as HDF5 groups of 1D array datasets of equal length. Metadata is stored using top-level attributes.


Library
-------


The excellent `h5py <http://docs.h5py.org/en/latest/>`_ Python interface to HDF5 provides direct access to the group and dataset structure of a cooler file. h5py translates HDF5 dataset queries directly into NumPy arrays.

The cooler library provides an additional thin wrapper over h5py to support creation and conversion of cooler files as well as both tabular and sparse matrix views on the data. Range queries can be made using either integer bin indexing or genomic interval strings. Table range queries are retrieved as Pandas DataFrames and Series. Matrix range queries are retrieved as SciPy sparse matrices. Metadata is retrieved as a json-compatible Python dictionary. The cooler library also includes utilities for performing contact matrix balancing on a cooler file of any resolution.

Try it out in a Jupyter notebook using `Binder <https://github.com/mirnylab/cooler-binder>`_.


Scripts
-------

See the `scripts <https://github.com/mirnylab/cooler/tree/master/scripts>`_ folder in the git repository for examples of how to aggregate, load, dump and balance contact matrices. Currently, data can be aggregated from "valid pairs" files stored as tabix-indexed text as well as sorted valid pairs HDF5 files (.frag) obtained from `hiclib <https://bitbucket.org/mirnylab/hiclib>`_.


HDF5 bindings in other languages
--------------------------------


- canonical C-library `libhdf5 <https://www.hdfgroup.org/HDF5/>`_
- C++: `C++ API <https://www.hdfgroup.org/HDF5/doc/cpplus_RM/>`_
- IDL: `bindings <http://www.harrisgeospatial.com/docs/routines-102.html>`_
- Java: `Java HDF5 Interface <https://www.hdfgroup.org/products/java/JNI3/jhi5/index.html>`_
- Julia: `HDF5.jl <https://github.com/JuliaIO/HDF5.jl>`_
- Mathematica: `API <http://reference.wolfram.com/language/ref/format/HDF.html>`_
- MATLAB: `high and low level API <http://www.mathworks.com/help/matlab/hdf5-files.html>`_
- node.js: `hdf5.node <https://github.com/HDF-NI/hdf5.node>`_
- Perl: `PDL::IO::HDF5 <http://search.cpan.org/~chm/PDL-IO-HDF5-0.6501/hdf5.pd>`_
- R: `rhdf5 <http://bioconductor.org/packages/release/bioc/html/rhdf5.html>`_
- Apache `Spark <https://hdfgroup.org/wp/2015/03/from-hdf5-datasets-to-apache-spark-rdds/>`_


Caveats
-------

HDF5 is not a database system and is not journalled. It supports concurrent read access but not simultaneous reads and writes (with upcoming support for the `SWMR <http://docs.h5py.org/en/latest/swmr.html>`_ access pattern). One must be careful using multi-process concurrency based on Unix ``fork()``: if a file is already open before the fork, the child processes will inherit state such that they won't play well with each other on that file. HDF5 will work fine with Python's ``multiprocessing`` as long as you make sure to close file handles before creating a process pool. Otherwise, you'll need to use locks, switch to the MPI parallel programming pattern using `Parallel HDF5 <http://docs.h5py.org/en/latest/mpi.html>`_, or avoid opening the file in worker processes completely (see this `blog post <http://assorted-experience.blogspot.ca/2013/11/h5py-and-multiprocessing.html>`_ for a simple workaround). For more information on using multiprocessing safely, see this `discussion <https://groups.google.com/forum/#!topic/h5py/bJVtWdFtZQM>`_.


Glossary
--------

HDF5 is a general purpose binary container format for large scientific datasets.

h5py is a Python library providing low-level bindings to the libhdf5 C-library and a high-level, numpy-aware API to interact with HDF5 files on disk.

Cooler is a flexible binary schema for Hi-C data based on a two-table sparse data model.

Cooler [Cooler-HDF5?] is also the name of an implementation of the Cooler schema using HDF5 as the container format.

cooler [pycooler?] is a Python package providing an API to create Cooler-HDF5 files and to interact with them both as data frames and sparse matrices.
