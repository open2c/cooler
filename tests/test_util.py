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

    # Punctuation in names (aside from :)
    assert parse_region_string('name-with-hyphens-') == ('name-with-hyphens-', None, None)
    assert parse_region_string('GL000207.1') == ('GL000207.1', None, None)
    assert parse_region_string('GL000207.1:1000-2000') == ('GL000207.1', 1000, 2000)

    # Trailing dash
    assert parse_region_string('chr21:1000-') == ('chr21', 1000, None)

    # Humanized units
    assert parse_region_string('6:1kb-2kb') == ('6', 1000, 2000)
    assert parse_region_string('6:1k-2000') == ('6', 1000, 2000)
    assert parse_region_string('6:1kb-2M') == ('6', 1000, 2000000)
    assert parse_region_string('6:1Gb-') == ('6', 1000000000, None)

    assert_raises(ValueError, parse_region_string, 'chr1:2,000-1,000')  # reverse selection
    assert_raises(ValueError, parse_region_string, 'chr1::1000-2000')  # more than one colon
