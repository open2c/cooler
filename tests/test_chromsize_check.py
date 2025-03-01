import io
from typing import TextIO  # Import TextIO for type hinting

import pandas as pd
import pytest
from numpy import argsort as argnatsort


# Updated read_chromsizes function with TextIO for Python 3.11+
def read_chromsizes(
    filepath_or: str | TextIO,  # Use TextIO for file-like objects
    name_patterns: tuple[str, ...] = (r"^chr[0-9]+$", r"^chr[XY]$", r"^chrM$"),
    all_names: bool = False,
    **kwargs,
) -> pd.Series:
    """
    Parse a `<db>.chrom.sizes or <db>.chromInfo.txt file from the UCSC
    database, where `db` is a genome assembly name.
    """
    # Handle .gz files
    if isinstance(filepath_or, str) and filepath_or.endswith(".gz"):
        kwargs.setdefault("compression", "gzip")

    # Read the chromosome size file into a DataFrame
    chromtable = pd.read_csv(
        filepath_or,
        sep="\t",
        usecols=[0, 1],
        names=["name", "length"],
        dtype={"name": str},
        **kwargs,
    )

    # Convert the "length" column to numeric values, coercing errors to NaN
    chromtable["length"] = pd.to_numeric(chromtable["length"], errors="coerce")

    # Raise an error if there are any invalid (NaN) lengths
    if chromtable["length"].isnull().any():
        raise ValueError(
            f"Chromsizes file '{filepath_or}' contains missing or invalid "
            "length values. Please ensure that the file is properly formatted "
            "as tab-delimited with two columns: sequence name and integer "
            "length. Check for extraneous spaces or hidden characters."
        )

    # Filter chromosomes by pattern and sort them
    if not all_names:
        parts = []
        for pattern in name_patterns:
            part = chromtable[chromtable["name"].str.contains(pattern)]
            part = part.iloc[argnatsort(part["name"])]
            parts.append(part)
        chromtable = pd.concat(parts, axis=0)

    # Set the chromosome names as the index
    chromtable.index = chromtable["name"].values
    return chromtable["length"]


# Test for the read_chromsizes function with invalid length value
def test_read_chromsizes_bad_input():
    broken_data = "chr1\t1000\nchr2\tbad_value\nchr3\t2000\n"
    broken_file = io.StringIO(broken_data)

    # Expect a ValueError due to the non-numeric value in the "length" column.
    with pytest.raises(ValueError, match="Chromsizes file '.*' invalid length values"):
        read_chromsizes(broken_file)


# Main function to run the test
if __name__ == "__main__":
    pytest.main([__file__])
