# Cooler

## A sparse binary format for Hi-C contact maps

The `cooler` format is based on HDF5. It will use the file extension `.coo`.

### Top-level attributes (metadata)
...Copied most of this from biom format spec.

```
id              : <string or null> a field that can be used to id a table (or null)
type            : <string> Table type (a controlled vocabulary)
bintype         : <string> "fixed" or "variable"
binsize         : <int or null> Size of bins in bp if bintype is fixed.
format-url      : <url> A string with a static URL providing format details
format-version  : <tuple> The version of the current format, major and minor
generated-by    : <string> Package and revision that built the table
creation-date   : <datetime> Date the table was built (ISO 8601 format)
```

### Required groups

```
contigs/              : <HDF5 group> chromosome table
bins/                 : <HDF5 group> genomic bin table
tracks/               : <HDF5 group> additional columns along genomic bins
matrix/               : <HDF5 group> contact matrix in COO format
indexes/              : <HDF5 group> stores indexes for fast lookup
    contig_to_bin/    : <HDF5 group> maps chromosome IDs to ranges of bin IDs
    bin_to_matrix/    : <HDF5 group> maps bin IDs to ranges of matrix record IDs
```

### Required datasets

All datasets are 1D arrays that represent table columns. Datasets in the same group must have the same length. The implicit primary key (ID) for each table is the 0-based array index. Genomic coordinates are assumed to be 0-based and intervals half-open (1-based ends).

```
contigs/name                    : <HDF5 dataset> <S32> chromosome name
contigs/length                  : <HDF5 dataset> <int32> chromosome length in bp

bins/chrom_id                   : <HDF5 dataset> <int32> bin chromosome id
bins/start                      : <HDF5 dataset> <int64> bin start coordinate (bp)
bins/end                        : <HDF5 dataset> <int64> bin end coorindate (bp)

matrix/bin1_id                  : <HDF5 dataset> <int64> matrix pixel index along 1st axis
matrix/bin2_id                  : <HDF5 dataset> <int64> matrix pixel index along 2nd axis
matrix/count                    : <HDF5 dataset> <float64> matrix pixel value

indexes/contig_to_bin/bin_lo    : <HDF5 dataset> <int64> start of range in bin table
indexes/contig_to_bin/bin_hi    : <HDF5 dataset> <int64> end of range in bin table

indexes/bin_to_matrix/mat_lo    : <HDF5 dataset> <int64> start of range in matrix table
indexes/bin_to_matrix/mat_hi    : <HDF5 dataset> <int64> end of range in matrix table
```

