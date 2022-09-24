from __future__ import annotations

from itertools import product

import numpy as np
from numpy.typing import NDArray
from atomman import Atoms, Box, System

from hsnf import smith_normal_form


def make_supercell(system: System, matrix: NDArray):
    """Transform basis vectors of system as ``(a' b' c') = system.box.vects.T @ matrix``.
    """
    D, P, _ = smith_normal_form(matrix)  # D = P @ matrix @ Q
    Pinv = np.around(np.linalg.inv(P)).astype(int)
    invariants = tuple(D.diagonal())

    # column-wise lattice vectors
    A = system.box.vects.T

    # distinct lattice points in sublattice corresponding to `matrix`
    points = []
    for factor in product(*[range(f) for f in invariants]):
        p = A @ Pinv @ np.array(factor)
        points.append(p)

    new_vects = (A @ matrix).T
    new_box = Box(vects=new_vects)

    new_natoms = system.atoms.natoms * len(points)
    new_pos = []
    new_atype = []
    for i in range(system.atoms.natoms):
        pos_i = system.atoms.pos[i]
        atype_i = system.atoms.atype[i]

        new_atype.extend([atype_i] * len(points))
        for p in points:
            new_pos.append(pos_i + p)

    new_atoms = Atoms(
        natoms=new_natoms,
        atype=new_atype,
        pos=new_pos,
    )

    new_system = System(
        atoms=new_atoms,
        box=new_box,
        scale=False,
        symbols=system.symbols,
    )
    return new_system


def rotate_system(system: System, R: NDArray):
    """Rotate cell such that ``A := system.box.vects.T`` -> ``A @ R``."""
    # column-wise lattice vectors
    A = system.box.vects.T

    new_vects = (A @ R).T
    new_box = Box(vects=new_vects)

    new_pos = []
    for i in range(system.atoms.natoms):
        pos_i = R @ system.atoms.pos[i]
        new_pos.append(pos_i)

    new_atoms = Atoms(
        natoms=system.atoms.natoms,
        atype=system.atoms.atype,
        pos=new_pos,
    )

    new_system = System(
        atoms=new_atoms,
        box=new_box,
        scale=False,
        symbols=system.symbols,
    )
    return new_system
