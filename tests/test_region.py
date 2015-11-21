# -*- coding: utf-8 -*-
from cooler.region import *
import nose


def test_bool_ops():
    a, b = parse_region(('chr1', 5, 10)), parse_region(('chr1', 15, 20))
    assert comes_before(a, b) == True
    assert comes_after(a, b) == False
    assert contains(a, b) == False
    assert overlaps(a, b) == False

    a, b = parse_region(('chr1', 5, 10)), parse_region(('chr1', 10, 20))
    assert comes_before(a, b) == True
    assert comes_after(a, b) == False
    assert contains(a, b) == False
    assert overlaps(a, b) == False

    a, b = parse_region(('chr1', 5, 10)), parse_region(('chr1', 6, 10))
    assert comes_before(a, b) == True
    assert comes_after(a, b) == False
    assert contains(a, b) == False
    assert overlaps(a, b) == True

    a, b = parse_region(('chr1', 5, 10)), parse_region(('chr1', 5, 10))
    assert comes_before(a, b) == False
    assert comes_after(a, b) == False
    assert contains(a, b) == False
    assert overlaps(a, b) == True

    a, b = parse_region(('chr1', 5, 10)), parse_region(('chr1', 0, 6))
    assert comes_before(a, b) == False
    assert comes_after(a, b) == True
    assert contains(a, b) == False
    assert overlaps(a, b) == True

    a, b = parse_region(('chr1', 5, 10)), parse_region(('chr1', 0, 5))
    assert comes_before(a, b) == False
    assert comes_after(a, b) == True
    assert contains(a, b) == False
    assert overlaps(a, b) == False

    a, b = parse_region(('chr1', 5, 10)), parse_region(('chr1', 0, 15))
    assert comes_before(a, b) == False
    assert comes_after(a, b) == False
    assert contains(a, b) == False
    assert overlaps(a, b) == True


def test_set_ops():
    a, b = parse_region(('chr1', 5, 15)), parse_region(('chr1', 10, 20))
    assert intersection(a, b) == Region('chr1', 10, 15)

    a, b = parse_region(('chr1', 5, 15)), parse_region(('chr1', 10, 20))
    assert union(a, b) ==  Region('chr1', 5, 20)

    a, b = parse_region(('chr1', 5, 10)), parse_region(('chr1', 15, 20))
    assert hull(a, b) == Region('chr1', 5, 20)

    a, b = parse_region(('chr1', 5, 15)), parse_region(('chr1', 10, 20))
    assert diff(a, b) == Region('chr1', 5, 10)

    a, b = parse_region(('chr1', 5, 15)), parse_region(('chr1', 10, 20))
    x, y, z = partition(a, b)
    assert x == Region('chr1', 5, 10)
    assert y == Region('chr1', 10, 15)
    assert z == Region('chr1', 15, 20)

