.. _data-model:

===============
What is cooler?
===============

Cooler is the :ref:`implementation <schema>` of a data model for genomically-labeled sparse 2D arrays (matrices) with identical axes in HDF5. It is also the name of the `Python package <https://github.com/mirnylab/cooler>`_ that supports the format.

We use the term genomically-labeled array to refer to a data structure that assigns unique quantitative values to tuples of *genomic bins* obtained from an interval partition of a reference genome assembly. The tuples of bins make up the coordinates of the array’s elements. By omitting elements possessing zero or no value, the representation becomes sparse.

Cooler was designed for the storage and manipulation of extremely large Hi-C datasets at any resolution, but is not limited to Hi-C data in any way.

Genomically-labeled arrays
==========================

We can describe two tabular representations of such data.

**BG2**

By extending the `bedGraph <https://genome.ucsc.edu/goldenPath/help/bedgraph.html>`_ format, we can encode a 2D array with the following header.

+--------+--------+------+--------+--------+------+-------+
| chrom1 | start1 | end1 | chrom2 | start2 | end2 | value | 
+--------+--------+------+--------+--------+------+-------+

Other bin-related attributes (e.g. X and Y) and be appended as columns X1, X2, Y1, Y2, and so on. One problem with this representation is that each bin-related attribute can be repeated many times throughout the table, leading to great redundancy.

.. note :: bedGraph is technically different from `BED <https://bedtools.readthedocs.io/en/latest/content/general-usage.html?highlight=bedpe#bed-format>`_: the former describes a quantitative track supported by **non-overlapping** intervals (a step function), while the latter describes genomic intervals with no such restrictions. BG2 is different from `BEDPE <https://bedtools.readthedocs.io/en/latest/content/general-usage.html?highlight=bedpe#bedpe-format>`_ in the same way: intervals on the same axis are non-overlapping and interval pairs are not repeated (describing a heatmap).


**COO**

A simple solution is to decompose or "normalize" the single table into two files. The first is a bin table that describes the genomic bin segmentation on both axes of the matrix
(in the one-dimensional bedGraph style). The second table contains single columns that reference the rows of the bin table, providing a condensed representation of the nonzero elements of the array. Conveniently, this corresponds to the classic coordinate list (COO) sparse matrix representation. This two-table representation is used as a text format by `HiC-Pro <http://nservant.github.io/HiC-Pro/RESULTS.html>`_.

+------------------------+
| bins                   |
+--------+--------+------+
| chrom  | start  | end  |
+--------+--------+------+

+---------------------------+
| elements                  |
+---------+---------+-------+
| bin1_id | bin2_id | value |
+---------+---------+-------+

The table of elements (non-zero pixels) is often too large to hold in memory, but for any small selection of elements we can reconstitute the bin-related attributes by "joining" the bin IDs against the bin table. We refer to this process as element *annotation*.

Data model
==========

We model a genomically-labeled sparse matrix using three tables. It corresponds to the bin and element (pixel) tables above. We include a third chromosome description table for completeness, and indexes to support random access.

Tables
------

chroms
^^^^^^

+ Required columns: ``name[, length]``
+ Order: *enumeration*

An semantic ordering of the chromosomes, scaffolds or contigs of the assembly that the data is mapped to. This information can be extracted from the bin table below, but is included separately for convenience. This enumeration is the intended ordering of the chromosomes as they would appear in a global genomic matrix. Additional columns can provide metadata on the chromosomes, such as their length.

bins
^^^^

+ Required columns: ``chrom, start, end [, weight]``
+ Order: ``chrom`` (*enum*), ``start``

An enumeration of the concatenated genomic bins that make up a single dimension or axis of the global genomic matrix. Genomic bins can be of fixed size or variable sizes (e.g. restriction fragments). A genomic bin is defined by the triple (chrom, start, end), where start is zero-based and end is 1-based. The order is significant: the bins are sorted by chromosome (based on the chromosome enumeration) then by start, and each genomic bin is implicitly endowed with a 0-based bin ID from this ordering (i.e., its row number in the table). A reserved but optional column called ``weight`` can store weights for normalization or matrix balancing. Additional columns can be added to describe other bin-associated properties such as additional normalization vectors or bin-level masks.

pixels
^^^^^^

+ Required columns: ``bin1_id, bin2_id, count``
+ Order: ``bin1_id, bin2_id``

The array is stored as a single table containing only the nonzero upper triangle elements, assuming the ordering of the bins given by the bin table. Each row defines a non-zero element of the genomic matrix. Additional columns can be appended to store pixel-associated properties such as pixel-level masks or filtered and transformed versions of the data. Currently, the pixels are sorted lexicographically by the bin ID of the 1st axis (matrix row) then the bin ID of the 2nd axis (matrix column).


Indexes
-------

The sort order on the pixels and types of indexing strategies that can be used are strongly related.
We stipulate that the records of the pixel table must be sorted lexicographically by the bin
ID along the first axis, then by the bin ID along the second axis. This way, the ``bin1_id`` column can
be substituted with its run length encoding, which serves as a lookup index for the rows of the ma-
trix. With this index, we obtain a compressed sparse row (CSR) sparse matrix representation.

Given an enumeration of chromosomes, the bin table must also be lexicographically sorted by chromosome then by start coordinate. Then similarly, the chrom column of the bin table will reference the rows of the chrom table, and can also be substituted with a run length encoding.


Container
=========

The reference implementation of this data model uses `HDF5 <https://www.hdfgroup.org/HDF5/>`_ as the container format. HDF5 is a hierarchical data format for homongenenously typed multidimensional arrays, which supports chunking, compression, and random access. The HDF5 file specification and open source standard library is maintained by the nonprofit HDF Group.

HDF5 files consist of three fundamental entities: groups, datasets, and attibutes. The hierarchical organization of an HDF5 file is conceptually analogous to a file system: *groups* are akin to directories and *datasets* (arrays) are akin to files. Additionally, key-value metadata can be attached to groups and datasets using *attributes*. The standard library provides the ability to access and manipulate these entities. There are bindings for virtually every platform and programming environment. To learn more in detail about HDF5, I recommend the book `HDF5 and Python <https://www.safaribooksonline.com/library/view/python-and-hdf5/9781491944981/ch01.html>`_ by Andrew Collette, the author of ``h5py``.

To implement the data model in HDF5, data tables are stored in a columnar representation as HDF5 groups of 1D array datasets of equal length. Metadata is stored using top-level attributes. See the :ref:`schema <schema>`.


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
- R: `rhdf5 <http://bioconductor.org/packages/release/bioc/html/rhdf5.html>`_, `h5 <https://cran.r-project.org/web/packages/h5/>`_
- Apache `Spark <https://hdfgroup.org/wp/2015/03/from-hdf5-datasets-to-apache-spark-rdds/>`_


Caveats
-------

HDF5 is not a database system and is not journalled. It supports concurrent read access but not simultaneous reads and writes (with upcoming support for the `SWMR <http://docs.h5py.org/en/latest/swmr.html>`_ access pattern). One must be careful using multi-process concurrency based on Unix ``fork()``: if a file is already open before the fork, the child processes will inherit state such that they won't play well with each other on that file. HDF5 will work fine with Python's ``multiprocessing`` as long as you make sure to close file handles before creating a process pool. Otherwise, you'll need to use locks or avoid opening the file in worker processes completely (see this `blog post <http://assorted-experience.blogspot.ca/2013/11/h5py-and-multiprocessing.html>`_ for a simple workaround). For more information on using multiprocessing safely, see this `discussion <https://groups.google.com/forum/#!topic/h5py/bJVtWdFtZQM>`_.


.. comment:

	Why model it this way?

	To balance the tradeoff between simplicity, terseness and flexibility in an attempt to stay `Zen <https://www.python.org/dev/peps/pep-0020/>`_. 

	+ The schema is flexible enough to describe a whole genome contact matrix, or any subset of a contact matrix, including single contig-contig tiles.
	+ Given the variety of ways we might want to read the data or add new columns, flatter is better than nested.
	+ For one, it makes the data much easier to stream and process in chunks, which ideal for many types of out-of-core algorithms on very large contact matrices.
	+ Separating bins (annotations of the axis labels) from pixels (the matrix data) allows for easy inclusion of bin-level properties without introducing redundancy.


	Note that this flat structure [combination of bin + pixel tables] also defines a companion plain text format, a simple serialization of the binary format. Two forms are possible:

	- Two-file: The bin table and pixel table are stored as separate tab-delimited files (BED file + sparse triple file). See the output format from `Hi-C Pro <http://nservant.github.io/HiC-Pro/RESULTS.html#intra-and-inter-chromosomal-contact-maps>`_.

	- Single-file ("merged"): The ``bin1_id`` and ``bin2_id`` columns of the pixel table are replaced with annotations from the bin table, suffixed with `1` or `2` accordingly (e.g. ``chrom1``, ``start1``, ``end1``, ``weight1``, etc.). The result is a 2D extension of the `bedGraph <https://genome.ucsc.edu/goldenpath/help/bedgraph.html>`_ track format.


.. comment:
	Notes
	~~~~~

	Column-oriented vs record-oriented tables
	^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	Why is the reference schema column-oriented?

	- Cheap column addition/removal.
	- Better compression ratios.
	- Blazingly fast I/O speed can be achieved with new compressors such as `blosc <http://www.blosc.org/>`_.
	- Easy to migrate to other column stores such as `bcolz <https://github.com/Blosc/bcolz>`_, Apache `Parquet <https://parquet.apache.org/>`_, and Apache `Arrow <http://blog.cloudera.com/blog/2016/02/introducing-apache-arrow-a-fast-interoperable-in-memory-columnar-data-structure-standard/>`_.

	There is a tradeoff between flexibility and number of read cycles required to fetch all columns of a table, however, a column-oriented schema is fully interchangeable with a record-oriented representation (e.g., traditional SQL databases, CSV files).


	Supporting a matrix “view”
	^^^^^^^^^^^^^^^^^^^^^^^^^^

	Indexes are stored as 1D datasets in a separate group. The current indexes can be thought of as run-length encodings of the ``bins/chrom`` and ``pixels/bin1_id`` columns, respectively.


	Limitations
	^^^^^^^^^^^

	A complete rectangular matrix “view” of the data must be modeled on top of this representation. 2D range queries must be computed with the help of indexes. The sort order on the pixels and types of indexing strategies that can be used are strongly related. This could be changed in future versions of the schema.


.. comment:

    genome-assembly : string
        Name of genome assembly;  default: "unknown".

    Good h5py examples:
    https://www.uetke.com/blog/python/how-to-use-hdf5-files-in-python/

.. comment:
  Implementation Notes
  ====================

  Having the ``bin1_offset`` index, the ``bin1_id`` column becomes redundant, but we keep it for convenience as it is extremely compressible. It may be dropped in future versions.
