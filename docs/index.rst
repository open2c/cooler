.. cooler documentation master file, created by
   sphinx-quickstart on Sun Jan 17 11:53:23 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Cooler
======

Cooler is a support library for a sparse, compressed, binary persistent storage format for Hi-C contact matrices, called cool or COOL.

Cooler aims to provide the following functionality:

* Generate contact matrices from contact lists at arbitrary resolutions.
* Store contact matrices efficiently in ``cool`` format.
* Perform out-of-core genome wide contact matrix normalization (a.k.a. balancing).
* Perform fast range queries on a contact matrix.
* Convert contact matrices between formats.
* Provide a clean and well-documented Python API to work with Hi-C data.


Contents:

.. toctree::
   :maxdepth: 2

   quickstart
   datamodel
   schema
   container
   api
   cli


* :ref:`genindex`
* :ref:`Glossary`