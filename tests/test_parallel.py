import os.path as op
from operator import add

import cooler
from cooler import parallel

testdir = op.realpath(op.dirname(__file__))
datadir = op.join(testdir, "data")


def test_datapipe():
    inputs = {'a': 1, 'b': 2, 'c': 3, 'd': 4}
    keys = ['a', 'b', 'c', 'd']

    dp = parallel.MultiplexDataPipe(inputs.get, keys, map)
    dp = dp.pipe(lambda x: x)
    dp = dp.pipe((lambda x, const: x * const), 10)
    dp = dp.pipe([
        lambda x: x + 1,
        lambda x: x - 1,
    ])
    out = dp.gather()
    assert out == [10, 20, 30, 40]
    out = dp.reduce(add, 0)
    assert out == 100
    assert sum(i for i in dp) == 100

    dp = parallel.MultiplexDataPipe(inputs.get, keys, map)
    dp = dp.prepare(lambda x: x)
    # prepare initializer modifies the function signature
    dp = (
        dp
        .pipe(lambda x0, x: x + 100)
        .pipe(lambda x0, x: x0)
    )
    out = dp.gather()
    assert out == [1, 2, 3, 4]


def test_chunkgetter():
    path = op.join(datadir, "toy.symm.upper.2.cool")
    clr = cooler.Cooler(path)
    lo, hi = 1, 3

    getter = parallel.chunkgetter(clr)
    chunk = getter((lo, hi))
    assert isinstance(chunk, dict)
    assert 'chroms' not in chunk
    assert 'bins' in chunk
    assert 'pixels' in chunk
    assert len(chunk['pixels']['bin1_id']) == 2

    getter = parallel.chunkgetter(clr, include_chroms=True)
    chunk = getter((lo, hi))
    assert isinstance(chunk, dict)
    assert 'chroms' in chunk
    assert 'bins' in chunk
    assert 'pixels' in chunk

    getter = parallel.chunkgetter(clr, use_lock=True)
    chunk = getter((lo, hi))
    assert isinstance(chunk, dict)
    assert len(chunk['pixels']['bin1_id']) == 2


def test_split():
    path = op.join(datadir, "toy.symm.upper.2.cool")
    clr = cooler.Cooler(path)
    parallel.split(clr, map, chunksize=2)
    parallel.split(clr, map, spans=[(0, 2), (2, 4)])
