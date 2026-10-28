import os
import os.path as op
import subprocess
import sys
import threading
from functools import partial
from operator import add

import pytest

import cooler
import cooler.cli.balance
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
import multiprocessing.synchronize as sync

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


# Functions run in pool workers must be importable, so they're defined here


def _square(x):
    return x * x


def _fail(x):
    raise ValueError(f"bad item {x}")


def _exit_on(target, func, item):
    # Simulate a worker killed mid-task, e.g. by the OOM killer
    if item == target:
        os._exit(1)
    return func(item)


def _try_lock(_):
    lock = parallel.get_lock()
    got = lock.acquire(timeout=0.5)
    if got:
        lock.release()
    return got


def _run_with_timeout(func, timeout=60):
    # A pool that waits forever for a dead worker would hang the test suite
    result = {}

    def target():
        try:
            result["value"] = func()
        except Exception as e:  # noqa: BLE001 (re-raised below)
            result["error"] = e

    thread = threading.Thread(target=target, daemon=True)
    thread.start()
    thread.join(timeout)
    assert not thread.is_alive(), "timed out: the pool is waiting for a dead worker"
    if "error" in result:
        raise result["error"]
    return result["value"]


def test_process_pool():
    with parallel._ProcessPool(2) as pool:
        assert list(pool.map(_square, range(10))) == [x * x for x in range(10)]
        assert sorted(pool.imap_unordered(_square, range(10))) == [
            x * x for x in range(10)
        ]
        with pytest.raises(ValueError, match="bad item 0"):
            list(pool.map(_fail, range(3)))
        with pytest.raises(ValueError, match="bad item"):
            list(pool.imap_unordered(_fail, range(3)))


def test_process_pool_is_lazy():
    # At most max_pending inputs are taken before the first result is used
    taken = []

    def items():
        for i in range(100):
            taken.append(i)
            yield i

    with parallel._ProcessPool(2, max_pending=3) as pool:
        for method in [pool.map, pool.imap_unordered]:
            taken.clear()
            results = method(_square, items())
            next(results)
            assert len(taken) <= 3
            results.close()


def test_process_pool_shares_lock():
    # Workers use the parent's lock, whatever the start method
    lock = parallel.get_lock()
    with parallel._ProcessPool(2) as pool:
        lock.acquire()
        try:
            assert not any(pool.map(_try_lock, range(2)))
        finally:
            lock.release()
        assert all(pool.map(_try_lock, range(2)))


@pytest.mark.parametrize("method", ["map", "imap_unordered"])
def test_process_pool_worker_dies(method):
    def run():
        with parallel._ProcessPool(2) as pool:
            return list(getattr(pool, method)(partial(_exit_on, 3, _square), range(8)))

    with pytest.raises(parallel.WorkerDiedError, match="ran out of memory"):
        _run_with_timeout(run)


def test_process_pool_start_fails(monkeypatch):
    # If the workers can't start, the executor is shut down
    shutdowns = []
    orig = parallel.ProcessPoolExecutor.shutdown

    def shutdown(self, *args, **kwargs):
        shutdowns.append(kwargs)
        return orig(self, *args, **kwargs)

    def fail(future):
        raise parallel.WorkerDiedError("no workers")

    monkeypatch.setattr(parallel.ProcessPoolExecutor, "shutdown", shutdown)
    monkeypatch.setattr(parallel._ProcessPool, "_result", staticmethod(fail))
    with pytest.raises(parallel.WorkerDiedError, match="no workers"):
        parallel._ProcessPool(2)
    assert shutdowns and shutdowns[0]["cancel_futures"]


class _DyingPool(parallel._ProcessPool):
    # Kills the worker that gets the first task
    def map(self, func, iterable):
        items = list(iterable)
        return super().map(partial(_exit_on, items[0], func), items)

    def imap_unordered(self, func, iterable):
        items = list(iterable)
        return super().imap_unordered(partial(_exit_on, items[0], func), items)


def test_coarsen_worker_dies(monkeypatch, tmp_path):
    monkeypatch.setattr(cooler._reduce, "_ProcessPool", _DyingPool)
    path = op.join(datadir, "hg19.GM12878-MboI.matrix.2000kb.cool")
    with pytest.raises(parallel.WorkerDiedError):
        _run_with_timeout(
            lambda: cooler.coarsen_cooler(
                path, str(tmp_path / "out.cool"), 2, chunksize=10_000, nproc=2
            )
        )


def test_balance_worker_dies(run_cli, monkeypatch):
    monkeypatch.setattr(cooler.cli.balance, "_ProcessPool", _DyingPool)
    path = op.join(datadir, "hg19.GM12878-MboI.matrix.2000kb.cool")
    result = _run_with_timeout(
        lambda: run_cli("balance", path, "--stdout", "-p", 2, exit_code=1)
    )
    assert "Error: A worker process exited unexpectedly" in result.stderr
    assert "Try fewer processes (-p) or smaller chunks." in result.stderr
