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
format-version  : <tuple> The version of the current biom format, major and minor
generated-by    : <string> Package and revision that built the table
creation-date   : <datetime> Date the table was built (ISO 8601 format)
```

### Required groups

```
chromosomes/        : <HDF5 group> chromosome table
bins/               : <HDF5 group> genomic bin table
tracks/             : <HDF5 group> additional columns along genomic bins
heatmap/            : <HDF5 group> contact matrix in COO format
indexes/            : <HDF5 group> stores indexes for fast lookup
    chrom2bin/      : <HDF5 group> maps chromosome IDs to ranges of bin IDs
    bin2heatmap/    : <HDF5 group> maps bin IDs to ranges of heatmap record IDs
```

### Required datasets

All datasets are 1D arrays that represent table columns. Datasets in the same group must have the same length. The implicit primary key (ID) for each table is the 0-based array index. Genomic coordinates are assumed to be 0-based and intervals half-open (1-based ends).

```
chromosomes/name                : <HDF5 dataset> chromosome name
chromosomes/length              : <HDF5 dataset> chromosome length in bp

bins/chrom_id                   : <HDF5 dataset> bin chromosome id
bins/start                      : <HDF5 dataset> bin start coordinate (bp)
bins/end                        : <HDF5 dataset> bin end coorindate (bp)

heatmap/bin1_id                 : <HDF5 dataset> heatmap pixel index along 1st axis
heatmap/bin2_id                 : <HDF5 dataset> heatmap pixel index along 2nd axis
heatmap/count                   : <HDF5 dataset> heatmap pixel value

indexes/chrom2bin/bin_lo        : <HDF5 dataset> start of range in bin table
indexes/chrom2bin/bin_hi        : <HDF5 dataset> end of range in bin table

indexes/bin2heatmap/heatmap_lo  : <HDF5 dataset> start of range in heatmap table
indexes/bin2heatmap/heatmap_hi  : <HDF5 dataset> end of range in heatmap table
```

