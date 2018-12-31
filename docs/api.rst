.. _api-reference:

API Reference
=============

.. toctree::
   :maxdepth: 1


Quick reference
---------------

Cooler class
~~~~~~~~~~~~
.. autosummary:: 
    cooler.Cooler
    cooler.Cooler.binsize
    cooler.Cooler.chromnames
    cooler.Cooler.chromsizes
    cooler.Cooler.bins
    cooler.Cooler.pixels
    cooler.Cooler.matrix
    cooler.Cooler.open
    cooler.Cooler.info
    cooler.Cooler.offset
    cooler.Cooler.extent

Creation/reduction
~~~~~~~~~~~~~~~~~~
.. autosummary:: 
    cooler.create_cooler
    cooler.merge_coolers
    cooler.coarsen_cooler
    cooler.zoomify_cooler

Manipulation
~~~~~~~~~~~~
.. autosummary:: 
    cooler.annotate
    cooler.balance_cooler
    cooler.rename_chroms

File operations
~~~~~~~~~~~~~~~
.. autosummary::
    cooler.fileops.is_cooler
    cooler.fileops.is_multires_file
    cooler.fileops.list_coolers
    cooler.fileops.cp
    cooler.fileops.mv
    cooler.fileops.ln

Sandbox
~~~~~~~

.. autosummary::
    cooler.sandbox.dask.read_table

----

cooler
------

.. autoclass:: cooler.Cooler
    :members:
.. autofunction:: cooler.annotate
.. autofunction:: cooler.create_cooler
.. autofunction:: cooler.merge_coolers
.. autofunction:: cooler.coarsen_cooler
.. autofunction:: cooler.zoomify_cooler

----

cooler.create
---------

.. autofunction:: cooler.create.sanitize_pixels
.. autofunction:: cooler.create.sanitize_records

cooler.fileops
---------

.. autofunction:: cooler.fileops.is_cooler
.. autofunction:: cooler.fileops.is_multires_file
.. autofunction:: cooler.fileops.list_coolers
.. autofunction:: cooler.fileops.cp
.. autofunction:: cooler.fileops.mv
.. autofunction:: cooler.fileops.ln

cooler.util
-----------

.. autofunction:: cooler.util.partition
.. autofunction:: cooler.util.fetch_chromsizes
.. autofunction:: cooler.util.read_chromsizes
.. autofunction:: cooler.util.binnify
.. autofunction:: cooler.util.digest

cooler.sandbox
--------------

.. autofunction:: cooler.sandbox.dask.read_table