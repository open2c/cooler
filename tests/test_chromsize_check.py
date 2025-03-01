import io

import pytest

from cooler.util import (
    read_chromsizes,  # Correctly import from the 'cooler.util' module
)


# Test for the read_chromsizes function with invalid length value (non-numeric value)
def test_read_chromsizes_bad_input():
    # Simulating a bad .chrom.sizes file (with non-numeric value in the 'length' column)
    broken_data = "chr1\t1000\nchr2\tbad_value\nchr3\t2000\n"
    broken_file = io.StringIO(broken_data)

    # Expect a ValueError due to the non-numeric value in the "length" column.
    with pytest.raises(ValueError, match=r"Chromsizes file '.*' contains missing or invalid length values"):
        read_chromsizes(broken_file)


# Test for the read_chromsizes function with space delimiter instead of tab delimiter
def test_read_chromsizes_bad_delimiter():
    # Simulating a .chrom.sizes file with space delimiters (instead of tabs)
    broken_data = "chr1 1000\nchr2 bad_value\nchr3 2000\n"
    broken_file = io.StringIO(broken_data)

    # Expect a ValueError due to space delimiter (not tab) in the file
    with pytest.raises(ValueError, match=r"Chromsizes file '.*' uses spaces instead of tabs as delimiters. Please use tabs."):
        read_chromsizes(broken_file)


# Main function to run the tests
if __name__ == "__main__":
    pytest.main([__file__])
