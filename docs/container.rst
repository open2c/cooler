Container
=========

The reference implementation of the Cooler schema uses `HDF5 <https://www.hdfgroup.org/HDF5/>`_ as the container format and uses the extension ``.cool``. HDF5 is a hierarchical data format for homongenenously typed multidimensional arrays, which supports chunking, compression, and random access. The HDF5 file specification and open source standard library is maintained by the nonprofit HDF Group.

HDF5 files consist of three fundamental entities: groups, datasets, and attibutes. The hierarchical organization of an HDF5 file is conceptually analogous to a file system: *groups* are akin to directories and *datasets* (arrays) are akin to files. Additionally, key-value metadata can be attached to groups and datasets using *attributes*. The standard library provides the ability to access and manipulate these entities. There are bindings for virtually every platform and programming environment. To learn more in detail about HDF5, I recommend the book `HDF5 and Python <https://www.safaribooksonline.com/library/view/python-and-hdf5/9781491944981/ch01.html>`_ by Andrew Collette, the author of ``h5py``.

To implement the Cooler data model in HDF5, data tables are stored in a columnar representation as HDF5 groups of 1D array datasets of equal length. Metadata is stored using top-level attributes.



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
- R: `rhdf5 <http://bioconductor.org/packages/release/bioc/html/rhdf5.html>`_
- Apache `Spark <https://hdfgroup.org/wp/2015/03/from-hdf5-datasets-to-apache-spark-rdds/>`_


Caveats
-------

HDF5 is not a database system and is not journalled. It supports concurrent read access but not simultaneous reads and writes (with upcoming support for the `SWMR <http://docs.h5py.org/en/latest/swmr.html>`_ access pattern). One must be careful using multi-process concurrency based on Unix ``fork()``: if a file is already open before the fork, the child processes will inherit state such that they won't play well with each other on that file. HDF5 will work fine with Python's ``multiprocessing`` as long as you make sure to close file handles before creating a process pool. Otherwise, you'll need to use locks, switch to the MPI parallel programming pattern using `Parallel HDF5 <http://docs.h5py.org/en/latest/mpi.html>`_, or avoid opening the file in worker processes completely (see this `blog post <http://assorted-experience.blogspot.ca/2013/11/h5py-and-multiprocessing.html>`_ for a simple workaround). For more information on using multiprocessing safely, see this `discussion <https://groups.google.com/forum/#!topic/h5py/bJVtWdFtZQM>`_.