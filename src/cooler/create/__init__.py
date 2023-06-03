import numpy as np

MAGIC = "HDF5::Cooler"
MAGIC_SCOOL = "HDF5::SCOOL"
MAGIC_MCOOL = "HDF5::MCOOL"

URL = "https://github.com/open2c/cooler"
CHROM_DTYPE = np.dtype("S")
CHROMID_DTYPE = np.int32
CHROMSIZE_DTYPE = np.int32
COORD_DTYPE = np.int32
BIN_DTYPE = np.int64
COUNT_DTYPE = np.int32
CHROMOFFSET_DTYPE = np.int64
BIN1OFFSET_DTYPE = np.int64
PIXEL_FIELDS = ("bin1_id", "bin2_id", "count")
PIXEL_DTYPES = (("bin1_id", BIN_DTYPE), ("bin2_id", BIN_DTYPE), ("count", COUNT_DTYPE))


from ._create import (
    append,
    create,
    create_cooler,
    create_from_unordered,
    create_scool,
    rename_chroms,
)
from ._ingest import (
    ArrayLoader,
    BadInputError,
    ContactBinner,
    HDF5Aggregator,
    PairixAggregator,
    TabixAggregator,
    aggregate_records,
    sanitize_pixels,
    sanitize_records,
    validate_pixels,
)
