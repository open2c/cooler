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
    cooler.create_scool

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

Logging
~~~~~~~
Cooler logs progress and diagnostics under the ``cooler`` logger, and doesn't
configure where its records go: they propagate to the root logger, so they
follow your own logging setup (e.g. ``logging.basicConfig``). Problems you can
act on are issued as warnings with :mod:`warnings`. To print cooler's log
messages without configuring logging, call :func:`cooler.log_to_stderr`.

.. autosummary::
    cooler.log_to_stderr

.. Sandbox
.. ~~~~~~~

.. .. autosummary::
..     cooler.sandbox.dask.read_table

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
.. autofunction:: cooler.balance_cooler
.. autofunction:: cooler.rename_chroms
.. autofunction:: cooler.create_scool
.. autofunction:: cooler.log_to_stderr


----

cooler.create
-------------

.. autofunction:: cooler.create.sanitize_pixels
.. autofunction:: cooler.create.sanitize_records

cooler.fileops
--------------

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

.. cooler.sandbox
.. --------------

.. .. autofunction:: cooler.sandbox.dask.read_table
