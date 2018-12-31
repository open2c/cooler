# -*- coding: utf-8 -*-
"""
Experimental API for developing split-apply-combine style algorithms on coolers.

"""
from __future__ import division, print_function
from functools import partial, reduce
from multiprocess import Pool, Lock

import numpy as np
import pandas
import h5py

from .util import partition
from .core import get

__all__ = ['partition', 'split', 'lock']

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
which leads to a race condition that causes simultaneous access to fail. One can
either use a lock to prevent the race condition, or close and re-open the file
in the workers after the fork.

See also:
* <https://groups.google.com/forum/#!topic/h5py/bJVtWdFtZQM>
* <https://github.com/h5py/h5py/issues/591#issuecomment-116785660>.

"""
lock = Lock()


def apply_pipeline(funcs, prepare, get, key):
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


class MultiplexDataPipe(object):
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
    Python's multiprocessing module uses Pickle for serialization, which has
    several limitations. Consider using a parallel map implementation that uses
    a more versatile serializer, such as dill or cloudpickle.

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
    def __init__(self, get, keys, map):
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

    def __copy__(self):
        other = self.__class__(self.get, self.keys, self.map)
        other.funcs = list(self.funcs)
        other._prepare = self._prepare
        return other

    def __reduce__(self):
        d = self.__dict__.copy()
        d.pop('map', None)
        return d

    def __iter__(self):
        return iter(self.run())

    def prepare(self, func):
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

    def pipe(self, func, *args, **kwargs):
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

    def run(self):
        """
        Run the pipeline

        Output depends on map implementation.

        """
        pipeline = partial(apply_pipeline, self.funcs, self._prepare, self.get)
        return self.map(pipeline, self.keys)

    def gather(self, combine=list, *args, **kwargs):
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

    def reduce(self, binop, init):
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


class chunkgetter(object):
    def __init__(self, clr, include_chroms=False, include_bins=True, use_lock=False):
        self.cooler = clr
        self.include_chroms = include_chroms
        self.include_bins = include_bins
        self.use_lock = use_lock

    def __call__(self, span):
        lo, hi = span
        chunk = {}
        try:
            if self.use_lock:
                lock.acquire()
            with self.cooler.open('r') as grp:
                if self.include_chroms:
                    chunk['chroms'] = get(grp['chroms'], as_dict=True)
                if self.include_bins:
                    chunk['bins'] = get(grp['bins'], as_dict=True)
                chunk['pixels'] = get(grp['pixels'], lo, hi, as_dict=True)
        finally:
            if self.use_lock:
                lock.release()
        return chunk

def split(clr, map=map, chunksize=int(10e6), spans=None, **kwargs):
    if spans is None:
        spans = partition(0, clr.info['nnz'], chunksize)
    return MultiplexDataPipe(chunkgetter(clr, **kwargs), spans, map)
