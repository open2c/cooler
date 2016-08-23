Library
=======


The excellent `h5py <http://docs.h5py.org/en/latest/>`_ Python interface to HDF5 provides direct access to the group and dataset structure of a cooler file. h5py translates HDF5 dataset queries directly into NumPy arrays.

The cooler library provides an additional thin wrapper over h5py to support creation and conversion of cooler files as well as both tabular and sparse matrix views on the data. Range queries can be made using either integer bin indexing or genomic interval strings. Table range queries are retrieved as Pandas DataFrames and Series. Matrix range queries are retrieved as SciPy sparse matrices. Metadata is retrieved as a json-compatible Python dictionary. The cooler library also includes utilities for performing contact matrix balancing on a cooler file of any resolution.

Try it out in a Jupyter notebook using `Binder <https://github.com/mirnylab/cooler-binder>`_.


- Create a cooler file
- Access a cooler file
- Balance a cooler matrix
- Dump output


Scripts
-------

See the `scripts <https://github.com/mirnylab/cooler/tree/master/scripts>`_ folder in the git repository for examples of how to aggregate, load, dump and balance contact matrices. Currently, data can be aggregated from "valid pairs" files stored as tabix-indexed text as well as sorted valid pairs HDF5 files (.frag) obtained from `hiclib <https://bitbucket.org/mirnylab/hiclib>`_.


- Start with contact list
- Choose chromosome ordering
- Generate a table of bins
- Order, sort the contact list and index it
- Aggregate