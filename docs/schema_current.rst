.. _current-version:

Schema Version: 3
-----------------

This schema describes a `compressed sparse row <https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_.28CSR.2C_CRS_or_Yale_format.29>`_ storage scheme (CSR) for a matrix with genomic dimension/axis annotations.

Data collection
~~~~~~~~~~~~~~~

We refer to the object hierarchy describing a single contact matrix as a cooler *data collection*.


URI syntax
~~~~~~~~~~

We identify a cooler data collection using a **URI string**, separating the system path to the container file from the data path within the container file by a double colon ``::``.

::
  
  /path/to/container.cool::/path/to/cooler/group

Tables
~~~~~~

HDF5 does not natively support sparse arrays or relational data structures: its datasets are dense multidimensional arrays. As groups of 1D arrays, tables and indexes can be represented using the `datashape <http://datashape.readthedocs.org/en/latest/>`_ layout language. The descriptions below describe required groups and arrays, conventional column orders, and default data types.

A **table** is a group of equal-length 1D arrays (HDF5 datasets) representing **columns**.

Additional groups and tables may be added to a data collection as long as they are not nested under the group of another table.

This storage mode does not enforce specific **column orders**, but conventional orders for *required* columns is provided in the listings below.

This storage mode does not set limits on the **number or length of columns**. Additional column arrays may be inserted into a table, but they must conform to the common length of the table.

The column **data types** are listed below as numpy equivalents. They are only defaults and may be altered as desired.

GZIP is chosen as the default **compression** filter for all columns. This is for portability reasons, since all versions of the HDF5 library ship with it.


chroms
""""""

::

    chroms: {
      # REQUIRED
      name:     typevar['Nchroms'] * string['ascii'],
      length:   typevar['Nchroms'] * int32
    }

In HDF5, ``name`` is a null-padded, fixed-length ASCII array, which maps to numpy's ``S`` dtype.

bins
""""

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

The ``cooler balance`` command by default stores balancing weights in a column called ``weight``. NaN values indicate genomic bins that were blacklisted during the balancing procedure.

pixels
""""""

::

    pixels: {
      # REQUIRED
      bin1_id:  typevar['Nnz'] * int64,
      bin2_id:  typevar['Nnz'] * int64,

      # RESERVED
      count:    typevar['Nnz'] * int32
    }

The ``count`` column is integer by default, but floating point types can be substituted. Additional columns are to be interpreted as supplementary value columns.

Indexes
~~~~~~~

Indexes are stored as 1D arrays in a separate group called ``indexes``. They can be thought of as run-length encodings of the ``bins/chrom`` and ``pixels/bin1_id`` columns, respectively. ``chrom_offset`` : indicates what row in the bin table each chromosome first appears. ``bin1_offset`` : indicates what row in the pixel table each bin1 ID appears. This is often called *indptr* in CSR data structures. Both arrays are required.

::

    indexes: {
      chrom_offset:  (typevar['Nchroms'] + 1) * int64,
      bin1_offset:   (typevar['Nbins'] + 1) * int64
    }

Sparse array interface and symmetry
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TODO


Metadata
~~~~~~~~

Essential key-value properties are stored as root-level HDF5 attributes. A specific bucket called ``metadata`` is reserved for arbitrary JSON-compatible user metadata.

All scalar string attributes, including serialized JSON, must be stored as variable-length UTF-8 encoded strings [*]_. 

::

    format : string
        'HDF5::Cooler'

    format-version : int
        The version of the format used

    bin-type : { "fixed" | "variable" }
        Indicates whether the resolution is constant along both axes.

    bin-size : int or "none"
        Size of bins in base pairs if bin-type is "fixed".

    symmetric-storage-mode : { "upper" | "none" }
        Indicates whether a symmetric matrix is stored using only upper triangular elements, including the diagonal.

    generated-by : string
        Agent that created the file (e.g. 'cooler-x.y.z').

    creation-date : datetime string
        Moment the file was built.

    metadata : JSON
        Custom user metadata about the experiment.

Additional metadata may be stored in the attributes of table columns or groups.

.. [*] In h5py, assigning a Python text string (Python 3 ``str`` or Python 2 ``unicode``) to an HDF5 attribute results in variable-length UTF-8 storage. When assigning attributes from h5py in Python 2, always use the ``unicode`` type.

Additional Notes
~~~~~~~~~~~~~~~~

Having the ``bin1_offset`` index, the ``bin1_id`` column becomes redundant, but we keep it for convenience as it is extremely compressible. It may be dropped in future versions.


Flavors
~~~~~~~

MCOOL


.. comment:

    genome-assembly : string
        Name of genome assembly;  default: "unknown".

    Good h5py examples:
    https://www.uetke.com/blog/python/how-to-use-hdf5-files-in-python/
