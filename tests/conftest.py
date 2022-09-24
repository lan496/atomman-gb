from __future__ import annotations

import atomman.unitconvert as uc
import numpy as np
import pytest
from atomman import Atoms, Box, System


def get_fcc(a: float) -> System:
    """
    A1 structure
    """
    vects = a * np.eye(3)
    box = Box(vects=vects)
    atype = 1
    pos = [
        [0, 0, 0],
        [0, 0.5, 0.5],
        [0.5, 0, 0.5],
        [0.5, 0.5, 0],
    ]
    atoms = Atoms(atype=atype, pos=pos)

    # set "scale=True" for fractional coordinates!
    system = System(atoms=atoms, box=box, scale=True)
    return system


@pytest.fixture
def system_fcc() -> System:
    a = uc.set_in_units(4.05, "angstrom")
    system = get_fcc(a)
    system.symbols = ("Al",)
    return system
