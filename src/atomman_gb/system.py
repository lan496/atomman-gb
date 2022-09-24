"""Manipulate Atomman's system."""
from __future__ import annotations

from itertools import product

import numpy as np
from atomman import Atoms, Box, System
from hsnf import smith_normal_form
from numpy.typing import NDArray


def make_supercell(system: System, matrix: NDArray):
    """Transform basis vectors of system as ``(a' b' c') = system.box.vects.T @ matrix``."""
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
    """Rotate cell such that A := system.box.vects.T -> R @ A."""
    new_vects = system.box.vects @ R.T
    new_box = Box(vects=new_vects)

    frac_coords = system.scale(system.atoms.pos)

    new_atoms = Atoms(
        natoms=system.atoms.natoms,
        atype=system.atoms.atype,
        pos=frac_coords,
    )
    new_system = System(
        atoms=new_atoms,
        box=new_box,
        scale=True,
        symbols=system.symbols,
    )
    return new_system


def normalize(system: System, rotation=None):
    """Normalize system w.r.t. LAMMPS convention."""
    # QR decomposition of lattice basis, gb_vects.T := (a, b, c) =: QR, where Q is orthogonal
    # matrix and R is upper triangular matrix. Then,
    # vects @ Q = ( Q^T @ (a,b,c) )^T = R^{T} is lower triangular, which is required for LAMMPS inputs.
    if rotation is None:
        vects = system.box.vects
        q, r = np.linalg.qr(vects.T)
        # change sign such that make diagonal terms in R^{T} positive
        q = q @ np.diag([(-1 if r[i, i] < 0 else 1) for i in range(3)])
        rotation = q.T
    new_system = rotate_system(system, rotation)

    # wrap atoms within unitcell
    frac_coords = np.remainder(new_system.scale(new_system.atoms.pos), 1)
    wrapped_pos = frac_coords @ new_system.box.vects
    new_system.atoms.pos = wrapped_pos

    return new_system, rotation
