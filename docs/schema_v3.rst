.. _version-3:

+------------------------+-----+
| **Schema Version**     |  3  |
+------------------------+-----+

The following document describes a `compressed sparse row (CSR) <https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_.28CSR.2C_CRS_or_Yale_format.29>`_ storage scheme for a matrix (i.e., a quantitative heatmap) with genomically labeled dimensions/axes.

HDF5 does not natively implement sparse arrays or relational data structures: its datasets are dense multidimensional arrays. We implement tables and sparse array indexes in HDF5 using groups of 1D arrays. The descriptions of tables and indexes in this document specify required groups and arrays, conventional column orders, and default data types.

.. admonition:: Summary of changes

    * Version 3 introduces the ``storage-mode`` metadata attribute to accomodate square matrices that are non-symmetric. Version 2 files which lack the ``storage-mode`` attribute should be interpreted as using the "symmetric-upper" storage mode. See `Storage mode`_.
    * The multi-resolution cooler file layout has been standardized. See `File flavors`_.



Data collection
===============

We refer to the object hierarchy describing a single matrix as a cooler *data collection*. A cooler data collection consists of **tables**, **indexes** and **metadata** describing a genomically-labelled sparse matrix.

A typical data collection has the following structure. At the top level, there are four `HDF5 Groups <http://docs.h5py.org/en/stable/high/group.html>`_, each containing 1D arrays (`HDF5 Datasets <http://docs.h5py.org/en/stable/high/dataset.html>`_). The depiction below shows an example group hierarchy as a tree, with arrays at the leaves, printed with their shapes in parentheses and their data type symbols.

::

  /
   ├── chroms
   │   ├── length (24,) int32
   │   └── name (24,) |S64
   ├── bins
   │   ├── chrom (3088281,) int32
   │   ├── start (3088281,) int32
   │   ├── end (3088281,) int32
   │   └── weight (3088281,) float64
   ├── pixels
   │   ├── bin1_id (271958554,) int64
   │   ├── bin2_id (271958554,) int64
   │   └── count (271958554,) int32
   └── indexes
       ├── bin1_offset (3088282,) int64
       └── chrom_offset (25,) int64

URI syntax
==========

We identify a cooler data collection using a **URI string** to its top-level group, separating the system path to the container file from the **group path** within the container file by a double colon ``::``.

::
  
  path/to/container.cool::/path/to/cooler/group

For any URI, the leading slash after the ``::`` may be omitted. To reference the root group ``/``, the entire ``::/`` suffix may be omitted (i.e., just a file path).

Tables
======

A **table** is a group of equal-length 1D arrays representing **columns**.

Additional groups and tables may be added to a data collection as long as they are not nested under the group of another table.

This storage mode does not enforce specific **column orders**, but conventional orders for *required* columns is provided in the listings below.

This storage mode does not set limits on the **number or length of columns**. Additional arrays may be inserted into a table to form new columns, but they must conform to the common length of the table.

The table descriptions below are given in the `datashape <http://datashape.readthedocs.org/en/latest/>`_ layout language. The column **data types** are given as numpy equivalents. They are only defaults and may be altered as desired.

GZIP is chosen as the default **compression** filter for all columns. This is for portability reasons, since all versions of the HDF5 library ship with it.

chroms
------

::

    chroms: {
      # REQUIRED
      name:     typevar['Nchroms'] * string['ascii'],
      length:   typevar['Nchroms'] * int32
    }

In HDF5, ``name`` is a null-padded, fixed-length ASCII array, which maps to numpy's ``S`` dtype.

bins
----

::

    bins: {
      # REQUIRED
      chrom:    typevar['Nbins'] * categorical[typevar['name'], type=string, ordered=True],
      start:    typevar['Nbins'] * int32,
      end:      typevar['Nbins'] * int32,

      # RESERVED
      weight:   typevar['Nbins'] * float64
    }

In HDF5, we use the integer-backed ENUM type to encode the ``chrom`` column. For data collections with a very large number of scaffolds, the ENUM type information may be too large to fit in the object's metadata header. In that case, the ``chrom`` column is stored using raw integers and the enumeration is inferred from the ``chrom`` table.

Genomic intervals are stored using a `0-start, half-open <http://genome.ucsc.edu/blog/the-ucsc-genome-browser-coordinate-counting-systems>`_ representation. The first interval in a scaffold should have ``start`` = 0 and the last interval should have ``end`` = the chromosome length. Intervals are sorted by ``chrom``, then by ``start``.

Because they measure the same quantity in the same units, the coordinate columns ``chroms/length``, ``bins/start`` and ``bins/end`` should be encoded using the same data type.

The :command:`cooler balance` command stores balancing weights in a column called ``weight`` by default. NaN values indicate genomic bins that were blacklisted during the balancing procedure.

pixels
------

::

    pixels: {
      # REQUIRED
      bin1_id:  typevar['Nnz'] * int64,
      bin2_id:  typevar['Nnz'] * int64,

      # RESERVED
      count:    typevar['Nnz'] * int32
    }

In the matrix coordinate system, ``bin1_id`` refers to the ith axis and ``bin2_id`` refers to the jth. Bin IDs are zero-based, i.e. we start counting at 0. Pixels are sorted by ``bin1_id`` then by ``bin2_id``.

The ``count`` column is integer by default, but floating point types can be substituted. Additional columns are to be interpreted as supplementary value columns.

.. warning:: `float16 <https://github.com/hetio/hetio/pull/15>`_ has limited support from 3rd party libraries and is not recommended. For floating point value columns consider using either single- (float32) or double-precision (float64).

Indexes
=======

Indexes are stored as 1D arrays in a separate group called ``indexes``. They can be thought of as run-length encodings of the ``bins/chrom`` and ``pixels/bin1_id`` columns, respectively. Both arrays are required.

::

    indexes: {
      chrom_offset:  (typevar['Nchroms'] + 1) * int64,
      bin1_offset:   (typevar['Nbins'] + 1) * int64
    }

* ``chrom_offset``: indicates which row in the bin table each chromosome first appears. The last element stores the length of the bin table.
* ``bin1_offset``: indicates which row in the pixel table each bin1 ID first appears. The last element stores the length of the pixel table. This index is usually called *indptr* in CSR data structures. 

Storage mode
============

Storing a symmetric matrix requires only the *upper triangular part, including the diagonal*, since the remaining elements can be reconstructed from the former ones. To indicate the use of this **mode of matrix storage** to client software, the value of the metadata attribute ``storage-mode`` must be set to ``"symmetric-upper"`` (see `Metadata`_). 

.. versionadded:: 3

    To indicate the absence of a special storage mode, e.g. for **non-symmetric** matrices, ``storage-mode`` must be set to ``"square"``.  This storage mode indicates to client software that 2D range queries should not be symmetrized.

.. warning:: In schema v2 and earlier, the symmetric-upper storage mode is always assumed.


Metadata
========

Essential key-value properties are stored as root-level `HDF5 attributes <http://docs.h5py.org/en/stable/high/attr.html>`_ in the data collection.

.. rubric:: Required attributes

.. describe:: format : string (constant)

    "HDF5::Cooler"

.. describe:: format-version : int

    The schema version used.

.. describe:: bin-type : { "fixed", "variable" }

    Indicates whether the resolution is constant along both axes.

.. describe:: bin-size : int or "null"

    Size of genomic bins in base pairs if bin-type is "fixed". Otherwise, "null".

.. describe:: storage-mode : { "symmetric-upper", "square" }

    Indicates whether ordinary sparse matrix encoding is used ("square") or whether a symmetric matrix is encoded by storing only the upper triangular elements ("symmetric-upper").

.. rubric:: Reserved, but optional

.. describe:: assembly : string

    Name of the genome assembly, e.g. "hg19".

.. describe:: generated-by : string

    Agent that created the file, e.g. "cooler-x.y.z".

.. describe:: creation-date : datetime string

    The moment the collection was created.

.. describe:: metadata : JSON

    Arbitrary JSON-compatible **user metadata** about the experiment.


All scalar string attributes, including serialized JSON, must be stored as **variable-length UTF-8 encoded strings**. 

.. warning:: When assigning scalar string attributes in Python 2, always store values having ``unicode`` type. In h5py, assigning a Python text string (Python 3 ``str`` or Python 2 ``unicode``) to an HDF5 attribute results in variable-length UTF-8 storage.

Additional metadata may be stored in other top-level attributes and the attributes of table groups and columns.


File flavors
============

Many cooler data collections can be stored in a single file. We recognize two common **layouts**:

* A single-resolution cooler file that contains a single data collection under the ``/`` group. Conventional file extension: ``.cool``.

::
  
  XYZ.1000.cool
  /
   ├── bins
   ├── chroms
   ├── pixels
   └── indexes


* A multi-resolution cooler file that contains multiple "coarsened" resolutions or "zoom-levels" derived from the same dataset. Multires cooler files should store each data collection underneath a group called ``/resolutions`` within a sub-group whose name is the bin size (e.g, ``XYZ.1000.mcool::resolutions/10000``). If the base cooler has variable-length bins, then use ``1`` to designate the base resolution, and the use coarsening multiplier (e.g. ``2``, ``4``, ``8``, etc.) to name the lower resolutions. Conventional file extension: ``.mcool``.

:: 

  XYZ.1000.mcool
  /
   └── resolutions
       ├── 1000
       │   ├── bins
       │   ├── chroms
       │   ├── pixels
       │   └── indexes
       ├── 2000
       │   ├── bins
       │   ├── chroms
       │   ├── pixels
       │   └── indexes
       ├── 5000
       │   ├── bins
       │   ├── chroms
       │   ├── pixels
       │   └── indexes
       ├── 10000
       │   ├── bins
       │   ├── chroms
       │   ├── pixels
       │   └── indexes
       .
       .
       .

.. note:: 

  The old multi-resolution layout used resolutions strictly in increments of *powers of 2*. In this layout, the data collections are named by zoom level, starting with ``XYZ.1000.mcool::0`` being the coarsest resolution up until the finest or "base" resolution (e.g., ``XYZ.1000.mcool::14`` for 14 levels of coarsening). 

  .. versionchanged:: 0.8
    Both the legacy layout and the new mcool layout are supported by `HiGlass <http://higlass.io/app/>`_. Prior to cooler 0.8, the new layout was produced only when requesting a specific list of resolutions. As of cooler 0.8, the new layout is always produced by the :command:`cooler zoomify` command unless the ``--legacy`` option is given. Files produced by :py:func:`cooler.zoomify_cooler`, `hic2cool <https://github.com/4dn-dcic/hic2cool/>`_, and the mcools from the `4DN data portal <https://data.4dnucleome.org/>`_ also follow the new layout.
