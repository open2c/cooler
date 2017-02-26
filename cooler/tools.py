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


# Lock prevents race condition if HDF5 file is already open before forking.
# See discussion: <https://groups.google.com/forum/#!topic/h5py/bJVtWdFtZQM>
# And <https://github.com/h5py/h5py/issues/591#issuecomment-116785660>.
lock = Lock()


def get_dict(h5, lo=0, hi=None, fields=None, convert_enum=True):
    grp = h5
    series = False
    if fields is None:
        fields = list(grp.keys())
    elif isinstance(fields, six.string_types):
        fields = [fields]
        series = True

    data = {}
    for field in fields:
        dset = grp[field]

        if convert_enum:
            dt = h5py.check_dtype(enum=dset.dtype)
        else:
            dt = None

        if dt is not None:
            data[field] = pandas.Categorical.from_codes(
                dset[lo:hi],
                sorted(dt, key=dt.__getitem__),
                ordered=True)
        elif dset.dtype.type == np.string_:
            data[field] = dset[lo:hi].astype('U')
        else:
            data[field] = dset[lo:hi]

    return data


class _applymany(object):
    def __init__(self, funcs, get):
        self.funcs = funcs
        self.get = get

    def __call__(self, span):
        chunk = self.get(span)
        for func in self.funcs:
            chunk = func(chunk)
        return chunk


class chunkgetter(object):
    def __init__(self, cooler_path, cooler_root='/', include_chroms=False, include_bins=True, use_lock=False):
        self.cooler_path = cooler_path
        self.cooler_root = cooler_root
        self.include_chroms = include_chroms
        self.include_bins = include_bins
        self.use_lock = use_lock

    def __call__(self, span):
        lo, hi = span
        chunk = {}

        try:
            if self.use_lock:
                lock.acquire()
            with h5py.File(self.cooler_path, 'r') as h5:
                grp = h5[self.cooler_root]
                if self.include_chroms:
                    chunk['chroms'] = get_dict(grp['chroms'])
                if self.include_bins:
                    chunk['bins'] = get_dict(grp['bins'])
                chunk['pixels'] = get_dict(grp['pixels'], lo, hi)
        finally:
            if self.use_lock:
                lock.release()
        
        return chunk


class ChunkedDataPipe(object):
    def __init__(self, spans, getter, mapper):
        self.funcs = []
        self.spans = list(spans)
        self.get = getter
        self.map_impl = mapper

    def __reduce__(self):
        d = self.__dict__.copy()
        d.pop('map_impl', None)
        return d
            
    def __iter__(self): 
        return iter(self.map_impl(
            _applymany(self.funcs, self.get), 
            self.spans))
    
    def __next__(self):
        return next(iter(self))

    def pipe(self, func, *args, **kwargs):
        if args or kwargs:
            self.funcs.append(partial(func, *args, **kwargs))
        else:
            try:
                self.funcs.extend(list(func))
            except TypeError:
                self.funcs.append(func)
        return self

    def pipev(self, *funcs):
        self.funcs.extend(funcs)
        return self

    def combine(self, func=list, *args, **kwargs):
        return func(iter(self), *args, **kwargs)
    
    def reduce(self, binop, init):
        return reduce(binop, iter(self), init)


def partition(start, stop, step):
    return ((i, min(i+step, stop))
                for i in range(start, stop, step))


def split(clr, map=map, chunksize=int(10e6), spans=None, **kwargs):
    if spans is None:
        spans = partition(0, clr.info['nnz'], chunksize)
    return ChunkedDataPipe(spans, chunkgetter(clr.filename, **kwargs), map)
