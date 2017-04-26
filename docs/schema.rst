Schema
======

.. include:: schema_current.rst


Previous versions
~~~~~~~~~~~~~~~~~

* :ref:`v1 <version-1>`

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
