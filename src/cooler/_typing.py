from __future__ import annotations

from collections.abc import Iterable
from typing import Callable, Optional, TypeVar, Union

import numpy as np
import pandas as pd

T = TypeVar('T')
U = TypeVar('U')
MapFunctor = Callable[[Callable[[T], U], Iterable[T]], Iterable[U]]
GenomicRangeSpecifier = Union[str , tuple[str, Optional[int], Optional[int]]]
GenomicRangeTuple = tuple[str, int, int]
Tabular = Union[pd.DataFrame, dict[str, np.ndarray]]

__all__ = ["MapFunctor", "GenomicRangeSpecifier", "GenomicRangeTuple", "Tabular"]
