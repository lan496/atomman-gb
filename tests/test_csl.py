from __future__ import annotations

import numpy as np

from atomman_gb.csl import get_csl
from atomman_gb.utils import is_integer_array


def test_csl():
    # Sigma=5, [001]
    R = (
        np.array(
            [
                [3, -4, 0],
                [4, 3, 0],
                [0, 0, 5],
            ]
        )
        / 5
    )
    sigma = 5
    csl = get_csl(R, sigma)

    assert is_integer_array(csl)
    assert is_integer_array(np.dot(np.linalg.inv(R), csl))
