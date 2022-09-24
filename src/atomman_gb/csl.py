"""Calculate coincidence-site lattice (CSL)."""
from __future__ import annotations

import numpy as np
from hsnf import column_style_hermite_normal_form
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
    #     = dual( L(I|R) )
    # Ref: https://cseweb.ucsd.edu/classes/wi10/cse206a/lec2.pdf
    dunion = np.concatenate([np.eye(3), R], axis=1)  # (3, 6)
    dbasis_hnf, _ = column_style_hermite_normal_form(np.around(dunion * sigma).astype(int))
    dbasis = dbasis_hnf[:, :3].astype(np.float_) / sigma  # (3, 3)
    # rbasis[:, i] is the i-th basis vector of CSL
    csl = dbasis @ np.linalg.inv(dbasis.T @ dbasis)

    # sanity check
    assert is_integer_array(csl)
    assert is_integer_array(np.dot(np.linalg.inv(R), csl))
    assert np.linalg.matrix_rank(csl) == 3

    csl = np.around(csl)

    return csl
