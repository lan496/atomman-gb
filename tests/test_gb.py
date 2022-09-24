from __future__ import annotations

from collections import Counter

import numpy as np
from atomman import System

from atomman_gb.gb import CubicGBGenerator, CubicGBInfo


def test_cubic_macro():
    # Adapted from Aimsgb
    # list of sigmas (<= 30) whose angles are less than 90 degree
    counts_expect = {
        (0, 0, 1): [
            5,
            13,
            17,
            25,
            29,
        ],
        (1, 1, 0): [
            3,
            9,
            11,
            17,
            19,
            27,
        ],
        (1, 1, 1): [
            3,
            7,
            13,
            19,
            21,
        ],
    }

    for uvw in [(0, 0, 1), (1, 1, 0), (1, 1, 1)]:
        gbinfo = CubicGBInfo(uvw=uvw, max_sigma=30)
        c = Counter([d["sigma"] for d in gbinfo.datum if d["theta"] < np.pi / 2])
        assert set(c.keys()) == set(counts_expect[uvw])


def test_rotation():
    uvw = (0, 0, 1)
    gbinfo = CubicGBInfo(uvw=uvw, max_sigma=5)
    data = gbinfo.datum[1]  # sigma5, 53.1 deg

    rotation_expect = np.array(
        [
            [3 / 5, -4 / 5, 0],
            [4 / 5, 3 / 5, 0],
            [0, 0, 1],
        ]
    )
    assert np.allclose(data["rotation"], rotation_expect)


def test(system_fcc: System):
    uvw = (0, 0, 1)
    gbinfo = CubicGBInfo(uvw=uvw, max_sigma=5)

    data = gbinfo.datum[0]  # sigma5, 39.97 deg

    gbgen = CubicGBGenerator.make_tilt(
        system_fcc, data["rotation"], uvw, data["csl"], data["symmetric_tilt"][0]["plane_mid"]
    )
    gbs = gbgen.generate()  # noqa
