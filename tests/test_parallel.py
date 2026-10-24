import os.path as op
import subprocess
import sys
from operator import add

import cooler
from cooler import parallel

testdir = op.realpath(op.dirname(__file__))
datadir = op.join(testdir, "data")


def test_datapipe():
    inputs = {"a": 1, "b": 2, "c": 3, "d": 4}
    keys = ["a", "b", "c", "d"]

    dp = parallel.MultiplexDataPipe(inputs.get, keys, map)
    dp = dp.pipe(lambda x: x)
    dp = dp.pipe((lambda x, const: x * const), 10)
    dp = dp.pipe(
        [
            lambda x: x + 1,
            lambda x: x - 1,
        ]
    )
    out = dp.gather()
    assert out == [10, 20, 30, 40]
    out = dp.reduce(add, 0)
    assert out == 100
    assert sum(i for i in dp) == 100

    dp = parallel.MultiplexDataPipe(inputs.get, keys, map)
    dp = dp.prepare(lambda x: x)
    # prepare initializer modifies the function signature
    dp = dp.pipe(lambda x0, x: x + 100).pipe(lambda x0, x: x0)
    out = dp.gather()
    assert out == [1, 2, 3, 4]


def test_chunkgetter():
    path = op.join(datadir, "toy.symm.upper.2.cool")
    clr = cooler.Cooler(path)
    lo, hi = 1, 3

    getter = parallel.chunkgetter(clr)
    chunk = getter((lo, hi))
    assert isinstance(chunk, dict)
    assert "chroms" not in chunk
    assert "bins" in chunk
    assert "pixels" in chunk
    assert len(chunk["pixels"]["bin1_id"]) == 2

    getter = parallel.chunkgetter(clr, include_chroms=True)
    chunk = getter((lo, hi))
    assert isinstance(chunk, dict)
    assert "chroms" in chunk
    assert "bins" in chunk
    assert "pixels" in chunk

    getter = parallel.chunkgetter(clr, use_lock=True)
    chunk = getter((lo, hi))
    assert isinstance(chunk, dict)
    assert len(chunk["pixels"]["bin1_id"]) == 2


def test_split():
    path = op.join(datadir, "toy.symm.upper.2.cool")
    clr = cooler.Cooler(path)
    parallel.split(clr, map, chunksize=2)
    parallel.split(clr, map, spans=[(0, 2), (2, 4)])


def test_lock_is_lazy():
    assert parallel.get_lock() is parallel.get_lock()
    assert parallel.lock is parallel.get_lock()


def test_import_without_semaphores():
    # Simulate environments without process-shared semaphores, like AWS Lambda
    script = """
import errno
import multiprocess.synchronize as sync

def _no_semaphores(self, *args, **kwargs):
    raise OSError(errno.ENOSYS, "Function not implemented")

sync.SemLock.__init__ = _no_semaphores

import cooler
import cooler.parallel

assert cooler.parallel._lock is None
"""
    result = subprocess.run(
        [sys.executable, "-c", script], check=False, capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr
