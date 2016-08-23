Intro
=====

Cooler is a sparse data model, schema and HDF5-based file format for high resolution Hi-C contact maps. 

The cooler format implements a simple data model that stores a high resolution contact matrix along with important auxiliary data such as scaffold information, genomic bin annotations, and basic metadata.

Why? As published Hi-C datasets increase in sequencing depth and resolution, a simple sparse representation lends itself better not only to the increasing storage demands but also to performing streaming and `out-of-core <https://en.wikipedia.org/wiki/Out-of-core_algorithm>`_ algorithms for analysis.