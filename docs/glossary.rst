.. _Glossary:

Glossary
--------

* HDF5 is a general purpose binary container format for large scientific datasets.

* h5py is a Python library providing low-level bindings to the libhdf5 C-library and a high-level, numpy-aware API to interact with HDF5 files on disk.

* The cooler **data model** is a flexible sparse data model for Hi-C and other genomically-labeled arrays.

* The cooler **schema** describes an implementation of the cooler data model using HDF5 as the underlying storage layer.

* Cooler files store one or more cooler **data collections**, each representing a genomically-labeled sparse array.

* Single-resolution cooler files are conventionally given the extension ``.cool``. Multi-resolution files are usually suffixed ``.mcool``.

* The *cooler* Python package provides an API to create cooler files and to interact with them both as data frames and sparse matrices.

* A genomic **pairs** list provides pointwise 2-tuples of single-bp genomic locations. In Hi-C this is also called a contact list.

* A genomic **matrix**, 2D array or heatmap assigns unique quantitative values to pairs of genomic intervals taken from a bin segmentation of a genome assembly. In Hi-C, a contact matrix is obtained by aggregating pairs.
