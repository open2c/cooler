from __future__ import annotations

from typing import Callable, Dict, Iterable, Optional, Tuple, TypeVar, Union

import numpy as np
import pandas as pd

T = TypeVar('T')
U = TypeVar('U')
MapFunctor = Callable[[Callable[[T], U], Iterable[T]], Iterable[U]]
GenomicRangeSpecifier = Union[str , Tuple[str, Optional[int], Optional[int]]]
GenomicRangeTuple = Tuple[str, int, int]
Tabular = Union[pd.DataFrame, Dict[str, np.ndarray]]

__all__ = ["GenomicRangeSpecifier", "GenomicRangeTuple", "MapFunctor", "Tabular"]
