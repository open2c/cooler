"""
Experimental API for developing split-apply-combine style algorithms on
coolers.

"""

from __future__ import annotations

import builtins
from collections import deque
from collections.abc import Iterable, Iterator, Sequence
from concurrent.futures import FIRST_COMPLETED, Future, ProcessPoolExecutor, wait
from concurrent.futures.process import BrokenProcessPool
from functools import partial, reduce
from multiprocessing import Lock, get_start_method
from typing import Any, Callable

from ._typing import MapFunctor
from .api import Cooler
from .core import get
from .util import partition

__all__ = ["get_lock", "partition", "split"]

"""
Two possible reasons for using a lock

(1) Prevent a concurrent process from opening an HDF5 file while the same
file is open for writing. In order for reading processes to obtain the correct
state, make sure the writing process finishes writing (flushes its buffers and
actually closes the file) before reading processes attempt to open it.
This explicit synchronization shouldn't be necessary if using the file in
SWMR mode.

See also:
* <https://support.hdfgroup.org/HDF5/hdf5-quest.html#grdwt>
* <https://support.hdfgroup.org/projects/SWMR>

(2) Synchronize file access when opened before a fork(). Fork-based (Unix)
multiprocessing and concurrent reading are compatible as long as the fork
happens before the child processes open the file. If an HDF5 file is already
open before forking, the child processes inherit the same global HDF5 state,
which leads to a race condition that causes simultaneous access to fail. One
can either use a lock to prevent the race condition, or close and re-open the
file in the workers after the fork.

See also:
* <https://groups.google.com/forum/#!topic/h5py/bJVtWdFtZQM>
* <https://github.com/h5py/h5py/issues/591#issuecomment-116785660>.

"""
_lock = None


def get_lock() -> Any:
    """
    Return the global lock used to synchronize HDF5 file access between
    processes.

    The lock is created on first use rather than on import, so that cooler can
    be imported in environments that don't support process-shared semaphores
    (e.g. AWS Lambda). To share it with fork-based worker processes, call this
    before creating the process pool.

    """
    global _lock
    if _lock is None:
        _lock = Lock()
    return _lock


def __getattr__(name: str) -> Any:
    # Keep `cooler.parallel.lock` working, but create it lazily
    if name == "lock":
        return get_lock()
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


class WorkerDiedError(RuntimeError):
    """A worker process exited without finishing its task."""


class WorkerStartError(WorkerDiedError):
    """The worker processes exited before running any task."""


def _start_error_message() -> str:
    msg = "Worker processes failed to start; see their error output above."
    method = get_start_method()
    if method in ("spawn", "forkserver"):
        # Each worker imports the main script, and an unguarded script starts
        # another pool from inside the worker
        msg += (
            f" With the '{method}' start method, a script that runs cooler with "
            'several processes must start them inside `if __name__ == "__main__":`.'
        )
    return msg


def _init_worker(lock: Any) -> None:  # pragma: no cover (runs in workers)
    # Workers started with spawn or forkserver don't inherit the parent's lock
    global _lock
    _lock = lock


def _noop() -> None:  # pragma: no cover (runs in workers)
    pass


class _ProcessPool:
    """
    A process pool for cooler's parallel tasks.

    Unlike ``multiprocessing.Pool``, which waits forever for the task of a
    worker that dies (e.g. killed for using too much memory), this raises
    :class:`WorkerDiedError`. Its ``map`` and ``imap_unordered`` are lazy and
    keep at most ``max_pending`` tasks in flight, so neither the inputs nor
    the results pile up in memory.

    Workers share the parent's lock (see :func:`get_lock`) whatever the start
    method, and are started when the pool is created: with fork, before the
    parent opens any more HDF5 files.
    """

    def __init__(self, nproc: int, max_pending: int | None = None):
        self.nproc = nproc
        self.max_pending = max_pending or 2 * nproc
        self._executor = ProcessPoolExecutor(
            nproc, initializer=_init_worker, initargs=(get_lock(),)
        )
        try:
            self._executor.submit(_noop).result()
        except BrokenProcessPool as e:
            # Not out of memory: a worker that fails before running anything
            # couldn't start, e.g. it failed to import the main script
            self._executor.shutdown(wait=True)
            raise WorkerStartError(_start_error_message()) from e
        except BaseException:
            self._executor.shutdown(wait=False, cancel_futures=True)
            raise

    def _result(self, future: Future) -> Any:
        try:
            return future.result()
        except BrokenProcessPool as e:
            # The executor fails the pending tasks before it stops the other
            # workers. Wait for them to exit, so none still has a file open
            # when the caller cleans up, e.g. removes a partly written cooler.
            self._executor.shutdown(wait=True)
            raise WorkerDiedError(
                "A worker process exited unexpectedly, most likely because it "
                "ran out of memory and was killed."
            ) from e

    def map(self, func: Callable, iterable: Iterable) -> Iterator[Any]:
        """Apply ``func`` to each item, yielding results in order."""
        items = iter(iterable)
        pending = deque()
        try:
            for item in items:
                pending.append(self._executor.submit(func, item))
                if len(pending) >= self.max_pending:
                    yield self._result(pending.popleft())
            while pending:
                yield self._result(pending.popleft())
        finally:
            self._stop(pending)

    def imap_unordered(self, func: Callable, iterable: Iterable) -> Iterator[Any]:
        """Apply ``func`` to each item, yielding results as they finish."""
        items = iter(iterable)
        pending = set()
        try:
            for item in items:
                pending.add(self._executor.submit(func, item))
                if len(pending) >= self.max_pending:
                    done, pending = wait(pending, return_when=FIRST_COMPLETED)
                    for future in done:
                        yield self._result(future)
            while pending:
                done, pending = wait(pending, return_when=FIRST_COMPLETED)
                for future in done:
                    yield self._result(future)
        finally:
            self._stop(pending)

    def _stop(self, pending: Iterable[Future]) -> None:
        # On an error or an early exit, cancel the queued tasks and wait for
        # the running ones, so no worker is still reading a file when the
        # caller handles the error
        running = [future for future in pending if not future.cancel()]
        wait(running)
        # If a worker died meanwhile, the executor failed these tasks before
        # stopping the other workers: wait for them to exit too
        if any(
            not future.cancelled() and isinstance(future.exception(), BrokenProcessPool)
            for future in running
        ):
            self._executor.shutdown(wait=True)

    def close(self) -> None:
        """Wait for running tasks, cancel queued ones, and stop the workers."""
        self._executor.shutdown(wait=True, cancel_futures=True)

    def __enter__(self) -> _ProcessPool:  # noqa: PYI034
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()


KeyType = Any


def apply_pipeline(
    funcs: list[Callable],
    prepare: Callable[[Any], Any] | None,
    get: Callable[[KeyType], Any],
    key: KeyType,
) -> Any:
    chunk = get(key)
    if prepare is not None:
        data = prepare(chunk)
        for func in funcs:
            data = func(chunk, data)
    else:
        data = chunk
        for func in funcs:
            data = func(data)
    return data


class MultiplexDataPipe:
    """
    Create an extendable pipeline of callables to be applied independently to
    each of a collection of inputs and produce a collection of outputs.

    New tasks are appended with the ``pipe`` method. Pipeline execution can be
    multiplexed using any ``map`` implementation, e.g. multiprocessing Pool.map
    or ipyparallel view.map for distributed execution.

    Depending on the ``map`` implementation results may be

    * yielded sequentially or online:
        Python 3 ``map`` or ``itertools.imap``, ``Pool.imap``,
        ``Pool.imap_unordered``

    * gathered and returned once all outputs are finished:
        Python 2 ``map``, ``Pool.map``, ``Pool.map_async``

    The pipeline can be run using one of:

    * ``gather``:
        Results are gathered and combined after all pipelines complete.

    * ``reduce``:
        Results are sequentially folded using a binary operator. This can
        save on memory when using a sequential or online ``map``
        implementation.

    Both methods are blocking, regardless of the blocking behavior of the
    ``map`` implementation (e.g., ``Pool.map_async``).

    Notes
    -----
    With a process pool from Python's multiprocessing module, the pipeline's
    functions are sent to the workers with pickle, so they must be picklable:
    e.g. module-level functions or ``functools.partial`` objects of them, but
    not lambdas or nested functions.

    See also
    --------
    http://stackoverflow.com/a/26521507 for a discussion of the differences
    between multiprocessing Pool implementations.

    Examples
    --------
    >>> X = np.arange(30)
    >>> spans = [(0, 10), (10, 20), (20, 30)]
    >>> dp = MultiplexDataPipe(lambda span: X[span[0]:span[1]], spans, map)
    >>> dp = dp.pipe(lambda x: x - 1).pipe(sum)
    >>> dp.gather()
    [35, 135, 235]
    >>> dp.reduce(lambda x, y: x + y, 0)
    405

    """

    def __init__(
        self,
        get: Callable[[KeyType], Any],
        keys: Iterable[KeyType],
        map: MapFunctor,
    ):
        """

        Parameters
        ----------
        get : callable
            Callable used to be used by workers that fetches the data
            corresponding to any of the provided keys

        keys : iterable
            Keys corresponding to input data

        map : callable
            Implementation of a map functor

        """
        self.get = get
        self.keys = list(keys)
        self.map = map
        self.funcs = []
        self._prepare = None

    def __copy__(self) -> MultiplexDataPipe:
        other = self.__class__(self.get, self.keys, self.map)
        other.funcs = list(self.funcs)
        other._prepare = self._prepare
        return other

    def __getstate__(self) -> dict:
        # The map function (e.g. a pool's) isn't picklable
        d = self.__dict__.copy()
        d.pop("map", None)
        return d

    def __setstate__(self, state: dict) -> None:
        self.__dict__.update(state)
        self.map = builtins.map

    def __iter__(self) -> Iterator[Any]:
        return iter(self.run())

    def prepare(self, func: Callable[[Any], Any]) -> MultiplexDataPipe:
        """
        Prepend a task that initializes the data for transformation.

        This optional step allows one to keep the original data chunk pristine.
        The callable ``func`` should return an initial object to pass along the
        pipeline for transformation. Subsequent callables in the pipeline will
        take two arguments instead of one:

        * chunk: original data chunk
        * data: transformed data passed along the pipeline

        Parameters
        ----------
        func : function/callable

        Returns
        -------
        A new datapipe with the initializer set.

        """
        self._prepare = func
        return self

    def pipe(
        self, func: Callable[[Any], Any] | Callable[[Any, Any], Any], *args, **kwargs
    ) -> MultiplexDataPipe:
        """
        Append new task(s) to the pipeline

        Parameters
        ----------
        func : function/callable or sequence of callables
            If a single function is provided, additional positional and keyword
            arguments can be provided and will be curried into the function.

        Returns
        -------
        A new datapipe with the additional task(s) appended.

        """
        other = self.__copy__()
        if args or kwargs:
            addon = [partial(func, *args, **kwargs)]
        else:
            try:
                addon = list(func)
            except TypeError:
                addon = [func]
        other.funcs += addon
        return other

    def run(self) -> Iterable[Any]:
        """
        Run the pipeline

        Output depends on map implementation.

        """
        pipeline = partial(apply_pipeline, self.funcs, self._prepare, self.get)
        return self.map(pipeline, self.keys)

    def gather(
        self, combine: Callable[[Iterable], Sequence] = list, *args, **kwargs
    ) -> Sequence[Any]:
        """
        Run the pipeline and gather outputs

        Parameters
        ----------
        combine : callable, optional
            Callable to consume the output. Default is builtin list.

        Returns
        -------
        Output of ``combine``

        """
        return combine(iter(self.run()), *args, **kwargs)

    def reduce(
        self,
        binop: Callable[[Any, Any], Any],
        init: Any,
    ) -> Any:
        """
        Run the pipeline and fold outputs cumulatively as they are returned

        Parameters
        ----------
        binop : binary operation
            A function of two arguments that returns a single value.
        init : object
            The initial value of the accumulation.

        Returns
        -------
        Reduced output

        """
        return reduce(binop, iter(self.run()), init)


class chunkgetter:
    def __init__(
        self,
        clr: Cooler,
        include_chroms: bool = False,
        include_bins: bool = True,
        use_lock: bool = False,
    ):
        self.cooler = clr
        self.include_chroms = include_chroms
        self.include_bins = include_bins
        self.use_lock = use_lock

    def __call__(self, span: tuple[int, int]) -> dict[str, Any]:
        lo, hi = span
        chunk = {}
        lock = get_lock() if self.use_lock else None
        try:
            if lock is not None:
                lock.acquire()
            with self.cooler.open("r") as grp:
                if self.include_chroms:
                    chunk["chroms"] = get(grp["chroms"], as_dict=True)
                if self.include_bins:
                    chunk["bins"] = get(grp["bins"], as_dict=True)
                chunk["pixels"] = get(grp["pixels"], lo, hi, as_dict=True)
        finally:
            if lock is not None:
                lock.release()
        return chunk


def split(
    clr: Cooler,
    map: MapFunctor = map,
    chunksize: int = 10_000_000,
    spans: Iterable[tuple[int, int]] | None = None,
    **kwargs,
) -> MultiplexDataPipe:
    if spans is None:
        spans = partition(0, int(clr.info["nnz"]), chunksize)
    return MultiplexDataPipe(get=chunkgetter(clr, **kwargs), keys=spans, map=map)
