class _IndexingMixin:
    def _unpack_index(self, key):
        if isinstance(key, tuple):
            if len(key) == 2:
                row, col = key
            elif len(key) == 1:
                row, col = key[0], slice(None)
            else:
                raise IndexError("invalid number of indices")
        else:
            row, col = key, slice(None)
        return row, col

    def _isintlike(self, num):
        try:
            int(num)
        except (TypeError, ValueError):
            return False
        return True

    def _process_slice(self, s, nmax):
        if isinstance(s, slice):
            if s.step not in (1, None):
                raise ValueError("slicing with step != 1 not supported")
            i0, i1 = s.start, s.stop
            if i0 is None:
                i0 = 0
            elif i0 < 0:
                i0 = nmax + i0
            if i1 is None:
                i1 = nmax
            elif i1 < 0:
                i1 = nmax + i1
            return i0, i1
        elif self._isintlike(s):
            if s < 0:
                s += nmax
            if s >= nmax:
                raise IndexError("index is out of bounds")
            return int(s), int(s + 1)
        else:
            raise TypeError("expected slice or scalar")


class RangeSelector1D(_IndexingMixin):
    """
    Selector for out-of-core tabular data. Provides DataFrame-like selection of
    columns and list-like access to rows.

    Examples
    --------

    Passing a column name or list of column names as subscript returns a new
    selector.

    >>> sel[ ['A', 'B'] ]  # doctest: +SKIP
    >>> sel['C']

    Passing a scalar or slice as subscript invokes the slicer.

    >>> sel[0]  # doctest: +SKIP
    >>> sel['A'][50:100]

    Calling the fetch method invokes the fetcher to parse the input into an
    integer range and then invokes the slicer.

    >>> sel.fetch('chr3:10,000,000-12,000,000') # doctest: +SKIP
    >>> sel.fetch(('chr3', 10000000, 12000000))

    """

    def __init__(self, fields, slicer, fetcher, nmax):
        self.fields = fields
        self._slice = slicer
        self._fetch = fetcher
        self._shape = (nmax,)

    @property
    def shape(self):
        return self._shape

    @property
    def columns(self):
        return self._slice(self.fields, 0, 0).columns

    @property
    def dtypes(self):
        return self._slice(self.fields, 0, 0).dtypes

    def keys(self):
        return list(self.columns)

    def __len__(self):
        return self._shape[0]

    def __contains__(self, key):
        return key in self.columns

    def __getitem__(self, key):
        # requesting a subset of columns
        if isinstance(key, (list, str)):
            return self.__class__(key, self._slice, self._fetch, self._shape[0])

        # requesting an interval of rows
        if isinstance(key, tuple):
            if len(key) == 1:
                key = key[0]
            else:
                raise IndexError("too many indices for table")
        lo, hi = self._process_slice(key, self._shape[0])
        return self._slice(self.fields, lo, hi)

    def fetch(self, *args, **kwargs):
        if self._fetch is not None:
            lo, hi = self._fetch(*args, **kwargs)
            return self._slice(self.fields, lo, hi)
        else:
            raise NotImplementedError


class RangeSelector2D(_IndexingMixin):
    """
    Selector for out-of-core sparse matrix data. Supports 2D scalar and slice
    subscript indexing.

    """

    def __init__(self, field, slicer, fetcher, shape):
        self.field = field
        self._slice = slicer
        self._fetch = fetcher
        self._shape = shape

    @property
    def shape(self):
        return self._shape

    def __len__(self):
        return self._shape[0]

    def __getitem__(self, key):
        s1, s2 = self._unpack_index(key)
        i0, i1 = self._process_slice(s1, self._shape[0])
        j0, j1 = self._process_slice(s2, self._shape[1])
        return self._slice(self.field, i0, i1, j0, j1)

    def fetch(self, *args, **kwargs):
        if self._fetch is not None:
            i0, i1, j0, j1 = self._fetch(*args, **kwargs)
            return self._slice(self.field, i0, i1, j0, j1)
        else:
            raise NotImplementedError
