from __future__ import absolute_import, print_function, division
import warnings
import numpy as np

MAGIC = u"HDF5::Cooler"
URL = u"https://github.com/mirnylab/cooler"
CHROM_DTYPE = np.dtype('S')
CHROMID_DTYPE = np.int32
CHROMSIZE_DTYPE = np.int32
COORD_DTYPE = np.int32
BIN_DTYPE = np.int64
COUNT_DTYPE = np.int32
CHROMOFFSET_DTYPE = np.int64
BIN1OFFSET_DTYPE = np.int64
PIXEL_FIELDS = ('bin1_id', 'bin2_id', 'count')
PIXEL_DTYPES = (('bin1_id', BIN_DTYPE),
                ('bin2_id', BIN_DTYPE),
                ('count', COUNT_DTYPE))


from ._ingest import (
    sanitize_pixels, validate_pixels, sanitize_records, aggregate_records,
    BadInputError, HDF5Aggregator, TabixAggregator, PairixAggregator,
    ArrayLoader, ContactBinner
)

from ._create import (
    create_cooler,
    create,
    create_from_unordered,
    append,
    rename_chroms
)
