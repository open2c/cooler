# -*- coding: utf-8 -*-
from __future__ import division, print_function

from nose.tools import assert_raises
from cooler.util import parse_region_string


def test_region_string_parser():
    # UCSC-style names
    assert parse_region_string('chr21') == ('chr21', None, None)
    assert parse_region_string('chr21:1000-2000') == ('chr21', 1000, 2000)
    assert parse_region_string('chr21:1,000-2,000') == ('chr21', 1000, 2000)

    # Ensembl style names
    assert parse_region_string('6') == ('6', None, None)
    assert parse_region_string('6:1000-2000') == ('6', 1000, 2000)
    assert parse_region_string('6:1,000-2,000') == ('6', 1000, 2000)

    # FASTA style names
    assert parse_region_string('gb|accession|locus') == ('gb|accession|locus', None, None)
    assert parse_region_string('gb|accession|locus:1000-2000') == ('gb|accession|locus', 1000, 2000)
    assert parse_region_string('gb|accession|locus:1,000-2,000') == ('gb|accession|locus', 1000, 2000)

    assert_raises(ValueError, parse_region_string, 'chr1:2,000-1,000')  # reverse selection
    assert_raises(ValueError, parse_region_string, 'name-with-hyphen')