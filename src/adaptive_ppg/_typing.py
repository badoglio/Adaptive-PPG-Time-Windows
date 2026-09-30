"""Array type aliases shared by the package."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int_]
BoolArray = NDArray[np.bool_]

__all__ = ["BoolArray", "FloatArray", "IntArray"]
