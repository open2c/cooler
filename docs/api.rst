.. _api-reference:

API Reference
=============

.. toctree::
   :maxdepth: 1


Quick reference
---------------

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
    cooler.annotate

.. autosummary:: 
    cooler.ice.iterative_correction
    cooler.reduce.merge
    cooler.reduce.coarsen
    cooler.reduce.zoomify

.. autosummary::
    cooler.io.create
    cooler.io.create_from_unordered
    cooler.io.rename_chroms
    cooler.io.sanitize_records
    cooler.io.sanitize_pixels

.. autosummary::
    cooler.io.is_cooler
    cooler.io.ls
    cooler.io.cp
    cooler.io.mv
    cooler.io.ln


Sandbox
~~~~~~~

.. autosummary::
    cooler.contrib.dask.daskify


cooler
------

.. autoclass:: cooler.Cooler
    :members:
.. autofunction:: cooler.annotate


cooler.io
---------

.. automodule:: cooler.io
    :members:
    :undoc-members:
    :show-inheritance:
.. autofunction:: cooler.io.rename_chroms
.. autofunction:: cooler.io.sanitize_pixels
.. autofunction:: cooler.io.sanitize_records


cooler.reduce
-------------

.. autofunction:: cooler.reduce.merge
.. autofunction:: cooler.reduce.coarsen
.. autofunction:: cooler.reduce.zoomify


cooler.ice
----------

.. autofunction:: cooler.ice.iterative_correction


cooler.util
-----------

.. autofunction:: cooler.util.fetch_chromsizes
.. autofunction:: cooler.util.read_chromsizes
.. autofunction:: cooler.util.binnify
.. autofunction:: cooler.util.digest
.. autofunction:: cooler.util.open_hdf5

cooler.contrib
--------------

.. automodule:: cooler.contrib.dask
    :members:
