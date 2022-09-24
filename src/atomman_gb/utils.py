from __future__ import annotations

from functools import reduce
from math import gcd

import numpy as np
from numpy.typing import NDArray


def is_integer_array(array: NDArray, atol=1e-8) -> bool:
    return np.allclose(array, np.around(array), atol=atol)


def gcd_on_list(nums: list[int]) -> int:
    return reduce(lambda x, y: gcd(x, y), nums, nums[0])


def extgcd(a, b):
    """
    Extended Euclidean algorithm for ax + by = gcd(a, b)
    Return (gcd(a, b), x, y)
    """
    if b == 0:
        return (a, 1, 0)
    else:
        g, xx, yy = extgcd(b, a % b)
        x = yy
        y = xx - yy * (a // b)
        return (g, x, y)


def find_smallest_multiplier(array: NDArray, max_denominator: int = 1000000):
    for a in range(1, max_denominator + 1):
        if is_integer_array(array * a):
            return a
    raise ValueError("Failed to make it integer!")
