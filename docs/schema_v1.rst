:orphan:

.. _version-1:

Schema
======

**Version: 1**

This schema describes a compressed sparse row storage scheme (CSR) for a *symmetric* matrix with genomic dimension/axis annotations.

Notes:

- Any number of additional optional columns can be added to each table. (e.g. normalization vectors, quality masks).
- Genomic coordinates are assumed to be 0-based and intervals half-open (1-based ends).


Contact matrix
~~~~~~~~~~~~~~

The tables and indexes can be represented in the `Datashape <http://datashape.readthedocs.org/en/latest/>`_ layout language:

::

    {
      chroms: {
        name:     typevar['Nchroms'] * string[32, 'ascii'],
        length:   typevar['Nchroms'] * int64,
      },
      bins: {
        chrom_id: typevar['Nbins'] * int32,
        start:    typevar['Nbins'] * int64,
        end:      typevar['Nbins'] * int64,
        weight:   typevar['Nbins'] * float64
      },
      pixels: {
        bin1_id:  typevar['Nnz'] * int32,
        bin2_id:  typevar['Nnz'] * int32,
        count:    typevar['Nnz'] * int32
      },
      indexes: {
        chrom_offset:  (typevar['Nchroms'] + 1) * int32,
        bin1_offset:   (typevar['Nbins'] + 1) * int32
      }
    }

Notes:

- Having the ``bin1_offset`` index, the ``bin1_id`` column becomes redundant, but we keep it for convenience as it is extremely compressible. It may be dropped in future versions.

Metadata
~~~~~~~~~

Essential key-value properties are stored as root-level HDF5 attributes. A specific bucket called ``metadata`` is reserved for arbitrary JSON-compatible user metadata.

::

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
    metadata        : <json> custom user metadata about the experiment


Indexes
~~~~~~~

Indexes are stored as 1D datasets in a separate group. The current indexes can be thought of as run-length encodings of the ``bins/chrom`` and ``pixels/bin1_id`` columns, respectively.

- ``chrom_offset`` : indicates what row in the bin table each chromosome first appears.
- ``bin1_offset`` : indicates what row in the pixel table each bin1 ID appears. This is often called *indptr* in CSR data structures.

