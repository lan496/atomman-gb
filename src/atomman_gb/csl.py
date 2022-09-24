"""Calculate coincidence-site lattice (CSL)."""
from __future__ import annotations

import numpy as np
from hsnf.lattice import compute_dual, compute_union
from numpy.typing import NDArray

from atomman_gb.utils import is_integer_array


def get_csl(R: NDArray, sigma: int) -> NDArray:
    """Return CSL from rotation matrix ``R`` which transform lattice-1 to lattice-2.

    Parameters
    ----------
    R: array, (3, 3)
        Rotation matrix in fractional coordinates
    sigma: int
        Least positive integer such that ``sigma * R`` and ``sigma * np.linalg.inv(R)`` are integer matrices.

    Returns
    -------
    csl: array, (3, 3)
        Intersection of lattice-1 and lattice-2
    """
    # Let lattice-1 be L1 = Z^3 and lattice-2 be L2 = R Z^3.
    # CSL = L1 cap L2
    #     = dual(dual(L[I]) + dual(L[R]))
    #     = dual( L([I|R]) )
    dbasis = compute_union(
        sigma * np.eye(3).astype(int), np.around(sigma * R).astype(int), row_wise=False
    )
    # csl[:, i] is the i-th basis vector of CSL
    csl = compute_dual(dbasis, row_wise=False) * sigma

    # sanity check
    assert is_integer_array(csl)
    assert is_integer_array(np.dot(np.linalg.inv(R), csl))
    assert np.linalg.matrix_rank(csl) == 3

    csl = np.around(csl)

    return csl
