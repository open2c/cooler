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
    append,
    rename_chroms
)


def ls(*args, **kwargs):
    warnings.warn(
        "`cooler.io` is deprecated in 0.8, will be removed in 0.9. "
        "Use `cooler.fileops.list_coolers()` instead.",
        category=FutureWarning, stacklevel=2)
    from ..fileops import list_coolers
    return list_coolers(*args, **kwargs)


def is_cooler(*args, **kwargs):
    warnings.warn(
        "`cooler.io` is deprecated in 0.8, will be removed in 0.9. "
        "Use `cooler.fileops.is_cooler()` instead.",
        category=FutureWarning, stacklevel=2)
    from ..fileops import is_cooler
    return is_cooler(*args, **kwargs)


def create(*args, **kwargs):
    warnings.warn(
        "`cooler.io.create()` is deprecated in 0.8, will be removed in 0.9. "
        "Use `cooler.create_cooler()` with ordered=True instead.",
        category=FutureWarning, stacklevel=2)
    from ._create import create
    return create(*args, **kwargs)


def create_from_unordered(*args, **kwargs):
    warnings.warn(
        "`cooler.io.create_from_unordered()` is deprecated in 0.8, "
        "will be removed in 0.9. Use `cooler.create_cooler()` instead.",
        category=FutureWarning, stacklevel=2)
    from ._create import create_from_unordered
    return create_from_unordered(*args, **kwargs)
