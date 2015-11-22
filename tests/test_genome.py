# -*- coding: utf-8 -*-
from cooler.genome import Region
import nose


def test_bool_ops():
    a, b = Region(('chr1', 5, 10)), Region(('chr1', 15, 20))
    assert a.comes_before(b) == True
    assert a.comes_after(b) == False
    assert a.contains(b) == False
    assert a.overlaps(b) == False

    a, b = Region(('chr1', 5, 10)), Region(('chr1', 10, 20))
    assert a.comes_before(b) == True
    assert a.comes_after(b) == False
    assert a.contains(b) == False
    assert a.overlaps(b) == False

    a, b = Region(('chr1', 5, 10)), Region(('chr1', 6, 10))
    assert a.comes_before(b) == True
    assert a.comes_after(b) == False
    assert a.contains(b) == True
    assert a.strictly_contains(b) == False
    assert a.overlaps(b) == True

    a, b = Region(('chr1', 5, 10)), Region(('chr1', 5, 10))
    assert a.comes_before(b) == False
    assert a.comes_after(b) == False
    assert a.contains(b) == True
    assert a.strictly_contains(b) == False
    assert a.overlaps(b) == True

    a, b = Region(('chr1', 5, 10)), Region(('chr1', 0, 6))
    assert a.comes_before(b) == False
    assert a.comes_after(b) == True
    assert a.contains(b) == False
    assert a.overlaps(b) == True

    a, b = Region(('chr1', 5, 10)), Region(('chr1', 0, 5))
    assert a.comes_before(b) == False
    assert a.strictly_comes_before(b) == False
    assert a.comes_after(b) == True
    assert a.strictly_comes_after(b) == True
    assert a.contains(b) == False
    assert a.overlaps(b) == False

    a, b = Region(('chr1', 5, 10)), Region(('chr1', 0, 15))
    assert a.comes_before(b) == False
    assert a.comes_after(b) == False
    assert a.contains(b) == False
    assert a.overlaps(b) == True


def test_set_ops():
    a, b = Region(('chr1', 5, 15)), Region(('chr1', 10, 20))
    assert a.intersection(b) == Region(('chr1', 10, 15))

    a, b = Region(('chr1', 5, 15)), Region(('chr1', 10, 20))
    assert a.combined(b) ==  Region(('chr1', 5, 20))

    a, b = Region(('chr1', 5, 10)), Region(('chr1', 15, 20))
    assert a.hull(b) == Region(('chr1', 5, 20))

    a, b = Region(('chr1', 5, 15)), Region(('chr1', 10, 20))
    assert a.diff(b) == Region(('chr1', 5, 10))


