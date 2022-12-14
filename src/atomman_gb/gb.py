"""Grain-boundary generator for cubic system."""
from __future__ import annotations

from itertools import permutations
from math import gcd
from typing import Any, Tuple

import numpy as np
from atomman import Atoms, Box, System
from numpy.typing import NDArray
from scipy.spatial import distance_matrix
from typing_extensions import TypeAlias  # for Python<3.10

from atomman_gb.csl import get_csl
from atomman_gb.system import make_supercell, normalize, rotate_system
from atomman_gb.utils import (
    extgcd,
    find_smallest_multiplier,
    gcd_on_list,
    is_integer_array,
)

Axis: TypeAlias = Tuple[int, int, int]
Plane: TypeAlias = Tuple[int, int, int]


class CubicGBInfo:
    """Grain-boundary information for cubic system."""

    def __init__(self, uvw: Axis, max_sigma: int):
        """Grain-boundary information for cubic system.

        Parameters
        ----------
        uvw: (int, int, int)
            Axis direction
        max_sigma: int
            Maximum number to search for sigma
        """
        self._uvw = uvw
        self._max_sigma = max_sigma
        assert gcd_on_list(self.uvw) == 1
        assert self.max_sigma > 1

        basis1, basis2 = self._get_planer_basis()
        self._planer_basis = np.array([basis1, basis2, self.uvw], dtype=np.float_).T

        angle_list = self._enumerate_cubic_csl_angle()
        datum = []
        for R, theta, sigma, _, _ in angle_list:
            csl = get_csl(R, sigma)

            # enumerate symmetric tilt planes
            symmetric_tilt_planes = self._enumerate_symmetric_tilt_plane(theta)

            datum.append(
                {
                    "sigma": sigma,
                    "theta": theta,
                    "rotation": R,
                    "csl": np.around(csl),  # transformation matrix from lattice-1 to CSL
                    "symmetric_tilt": symmetric_tilt_planes,  # [(plane_1, plane_2, plane_mid)]
                }
            )

        self._datum = sorted(datum, key=lambda d: (d["sigma"], d["theta"]))

    @property
    def uvw(self):
        """Return axis direction."""
        return self._uvw

    @property
    def max_sigma(self):
        """Maximum number to search for sigma."""
        return self._max_sigma

    @property
    def planer_basis(self) -> NDArray:
        """Return two basis vectors orthogonal to uvw and uvw."""
        return self._planer_basis

    @property
    def datum(self) -> list[dict[str, Any]]:
        """Return list of dictionary of grain-boundary information.

        Each dictionary has the following keys:
            - sigma
            - theta: angle of GB
            - rotation: rotation matrix in fractional coordinates
            - csl: transformation matrix from lattice-1 to CSL
            - csl_uvw: transformation matrix from lattice-1 to supercell of CSL parallel ot uvw
            - symmetric_tilt: [(plane_1, plane_2, plane_mid)]
        """
        return self._datum

    def _get_planer_basis(self):
        """
        Find basis vectors orthogonal to uvw.

        Ref: A. D. Banadaki and S. Patala, J. Appl. Cryst. 48, 2, (2015).
        """
        uvw = list(self.uvw[:])
        # swap axis such that l != 0
        swap = None
        if uvw[2] == 0:
            if uvw[0] != 0:
                swap = (0, 2)
            elif uvw[1] != 0:
                swap = (1, 2)
            else:
                raise ValueError("unreachable!")
        if swap:
            uvw[swap[0]], uvw[swap[1]] = uvw[swap[1]], uvw[swap[0]]

        gkl = gcd(abs(uvw[1]), abs(uvw[2]))
        basis1 = [0, -uvw[2] // gkl, uvw[1] // gkl]
        basis2 = []

        # let basis2 = [x, y, -(hx+ky)/l], volume of [uvw, basis1, basis2] is
        # |x| * (h**2 + k**2 + l**2) / gkl. To primitive cell, we search for x from one to greater.
        for x in range(1, gkl + 1):
            # Solve Diophantine equation, ky + lz = -hx for y and z
            if abs(uvw[0] * x) % gkl != 0:
                continue
            _, yy, zz = extgcd(uvw[1], uvw[2])
            y = -yy * uvw[0] * x // gkl
            basis2 = [x, y, -(uvw[0] * x + uvw[1] * y) // uvw[2]]
            break

        if swap:
            # swap back
            basis1[swap[0]], basis1[swap[1]] = basis1[swap[1]], basis1[swap[0]]
            basis2[swap[0]], basis2[swap[1]] = basis2[swap[1]], basis2[swap[0]]
            # reverse order for handedness
            basis1, basis2 = basis2, basis1

        # TODO: lattice reduction
        return [basis1, basis2]

    def _enumerate_cubic_csl_angle(self):
        """Enumerate CSL angles.

        O(sigma^3)

        Returns
        -------
        ret: list of (R, theta, sigma, m, n) where
            tan (theta/2) = n / m * sqrt(h ** 2 + k ** 2 + l ** 2)
        """
        norm2 = sum([self.uvw[i] ** 2 for i in range(3)])

        ret = []
        # skip sigma=1 because then CSL is equivalent to the original lattice!
        for sigma in range(2, self.max_sigma + 1):
            twosqrtsigma = int(np.ceil(2 * np.sqrt(sigma) + 0.5))
            for m in range(1, twosqrtsigma + 1):
                for n in range(1, twosqrtsigma + 1):
                    # m and n should be coprime
                    if gcd(m, n) != 1:
                        continue

                    # sigma is odd factor of S
                    S = m * m + norm2 * n * n
                    while S % 2 == 0:
                        S //= 2
                    if S != sigma or S > self.max_sigma:
                        continue

                    tan2 = n / m * np.sqrt(norm2)
                    theta = 2 * np.arctan(tan2)
                    R = self._get_rotation_matrix(theta)

                    # sigma is the least positive interger s.t. sigma * R and sigma * R^-1 are integer matrics.
                    assert is_integer_array(R * sigma)
                    assert is_integer_array(np.linalg.inv(R) * sigma)
                    g1 = gcd_on_list(np.around(R * sigma).reshape(-1).astype(int).tolist())
                    g2 = gcd_on_list(
                        np.around(np.linalg.inv(R) * sigma).reshape(-1).astype(int).tolist()
                    )
                    assert gcd(g1, g2) == 1

                    ret.append((R, theta, sigma, m, n))

        return ret

    def _get_rotation_matrix(self, theta: float) -> NDArray:
        rho = np.array(self.uvw) / np.linalg.norm(self.uvw)
        K = np.array(
            [
                [0, -rho[2], rho[1]],
                [rho[2], 0, -rho[0]],
                [-rho[1], rho[0], 0],
            ]
        )
        R = np.eye(3) + K * np.sin(theta) + np.dot(K, K) * (1 - np.cos(theta))
        # We do not need to transform to crystallographic coordinates for primitive cubic systems!
        return R

    def _enumerate_symmetric_tilt_plane(self, theta):
        # For symmetric tilt plane, (hkl) and R(hkl) are transformed to each other by mirror operations.
        # For primitive cubic systems, mirror planes are {100} and {110}
        mirror_planes = np.array(
            [
                [1, 0, 0],
                [0, 1, 0],
                [0, 0, 1],
                [1, 1, 0],
                [1, -1, 0],
                [1, 0, 1],
                [1, 0, -1],
                [0, 1, 1],
                [0, 1, -1],
            ]
        )

        def purify_hkl(hkl):
            hkl = hkl / np.min(np.abs(hkl[np.nonzero(hkl)]))
            mult = find_smallest_multiplier(hkl)
            hkl = np.around(hkl * mult).astype(int)
            return hkl

        ret = []
        used = set()
        for plane in mirror_planes:
            if tuple(plane) in used:
                continue
            if np.sum(np.array(self.uvw) * np.array(plane)) != 0:
                continue

            Rhalf = self._get_rotation_matrix(theta / 2)
            # Be careful this is contravariance!
            plane_1 = purify_vect(np.dot(Rhalf, plane))  # rotate -theta/2
            plane_2 = purify_vect(np.dot(Rhalf.T, plane))  # rotate theta/2

            ret.append(
                {
                    "plane_1": plane_1,
                    "plane_2": plane_2,
                    "plane_mid": plane,
                }
            )

            # visit equivalent planes
            for p in permutations(plane, r=3):
                used.add(tuple(p))

        return ret


class CubicGBGenerator:
    """Grain-boundary generator for cubic system."""

    def __init__(self, system: System, transform: NDArray, rotation: NDArray) -> None:
        """Grain-boundary generator for cubic system.

        Parameters
        ----------
        system
        transform:
            transformation matrix of given `system` to monoclinic grain boundary.
            Boundary plane is assumed to be perpendicular to z-axis.
        rotation:
            rotation matrix from lattice-1 to lattice-2
        """
        self._initial_system = system
        self._transform = transform
        self._rotation = rotation

        # system1: expand by `transform`
        self._system1 = make_supercell(self.initial_system, self._transform)

        # Transformation matrix: original system -> rotated -> system-2
        transform2 = np.linalg.solve(self._rotation, self._transform)
        assert is_integer_array(transform2)
        transform2 = np.around(transform2)

        # system2: Rotate initial system by theta and expand by `transform2`
        system2 = make_supercell(self.initial_system, transform2)
        self._system2 = rotate_system(system2, self._rotation)
        assert np.allclose(self._system1.box.vects, self._system2.box.vects)  # sublattice of CSL

        # Rotate axes for LAMMPS
        self._system1, rot = normalize(self._system1)
        self._system2, _ = normalize(self._system2, rot)

        # TODO: handle multiple inequivalent shifts for complex system!
        z1 = unique_zcoords(self._system1.atoms.pos[:, 2])
        self._system1_shift = np.array([0, 0, (z1[0] + z1[1]) / 2])
        z2 = unique_zcoords(self._system2.atoms.pos[:, 2])
        self._system2_shift = np.array([0, 0, (z2[0] + z2[1]) / 2])

    @property
    def initial_system(self) -> System:
        """Return initial system."""
        return self._initial_system

    def generate(
        self,
        expand_times: int = 1,
        ab_shifts: NDArray | None = None,
        c_thickness: float = 0,
        dmin: float | None = None,
    ) -> list[System]:
        """Generate grain boundary.

        Parameters
        ----------
        expand_times: int
        ab_shifts: array, (-1, 2)
        c_thickness: float
            width between system 1 and 2
        dmin: float
            delete one of too close atoms less than ``dmin``
        """
        if ab_shifts is None:
            ab_shifts = np.zeros((1, 2))

        # expand `2 * expand_times` times along c-axis
        expand = np.diag([1, 1, 2 * expand_times])
        system1 = make_supercell(self._system1, expand)
        system2 = make_supercell(self._system2, expand)

        # choose z <= 0.5 for system1
        eps = 1e-8
        frac_coords1 = np.remainder(
            system1.scale(system1.atoms.pos) + self._system1_shift[None, :], 1
        )
        select1 = frac_coords1[:, 2] < 0.5 + eps
        # choose z >= 0.5 for system2
        frac_coords2 = np.remainder(
            system2.scale(system2.atoms.pos) + self._system2_shift[None, :], 1
        )
        select2 = frac_coords2[:, 2] > 0.5 - eps

        new_vects = system1.box.vects[:]
        new_vects[2, 2] += c_thickness
        new_box = Box(vects=new_vects)

        pos1 = np.dot(frac_coords1[select1], system1.box.vects)
        select1_on_plane = np.logical_and(
            frac_coords1[:, 2] < 0.5 + eps, frac_coords1[:, 2] > 0.5 - eps
        )
        pos1_on_plane = np.dot(frac_coords1[select1_on_plane], system1.box.vects)
        atype1 = system1.atoms.atype[select1]

        new_symbols = system1.symbols + system2.symbols

        ret = []
        for xa, xb in ab_shifts:
            frac_coords2_shifted = frac_coords2[select2]
            frac_coords2_shifted[:, 0] += xa
            frac_coords2_shifted[:, 1] += xb
            pos2 = np.dot(frac_coords2_shifted, system2.box.vects)
            pos2[:, 2] += c_thickness

            # remove atoms too close pairs from `pos2`
            if dmin and (pos1_on_plane.shape[0] > 0):
                dists = distance_matrix(pos1_on_plane, pos2)
                mask = np.min(dists, axis=0) > dmin
                pos2 = pos2[mask, :]
            else:
                mask = np.ones(pos2.shape[0], dtype=bool)

            new_natoms = len(pos1) + len(pos2)
            new_pos = np.concatenate([pos1, pos2])

            atype2 = system2.atoms.atype[select2][mask] + len(system1.symbols)
            new_atype = np.array(atype1.tolist() + atype2.tolist())

            new_atoms = Atoms(
                natoms=new_natoms,
                atype=new_atype,
                pos=new_pos,
            )

            new_system = System(
                atoms=new_atoms,
                box=new_box,
                scale=False,
                symbols=new_symbols,
            )

            ret.append(new_system)

        return ret

    @classmethod
    def make_tilt(cls, system: System, rotation: NDArray, uvw: Axis, csl: NDArray, plane: Plane):
        """Construct GB generator for symmetric tilt GB."""
        axis_b = np.array(uvw)
        axis_c = np.array(plane)
        axis_a = np.cross(axis_b, axis_c)
        axis_a = axis_a // gcd_on_list(axis_a)
        B = np.transpose([axis_a, axis_b, axis_c])

        # Transformation matrix: original system -> system-1
        # We need to find transformation matrix M such that
        # csl @ M = B @ np.diag([l, m, n])
        Mtmp = np.linalg.solve(csl, B)
        M = Mtmp @ np.diag([find_smallest_multiplier(Mtmp[:, i]) for i in range(3)])
        transform = csl @ M

        return cls(system, transform=transform, rotation=rotation)

    @classmethod
    def make_twist(cls, system: System, rotation: NDArray, uvw: Axis, csl: NDArray):
        """Construct GB generator for twist GB."""
        axis_c = np.array(uvw)
        csl_uvw = find_csl_parallel_uvw(csl, uvw)
        axis_b = csl_uvw[:, 1]
        axis_a = np.cross(axis_b, axis_c)
        axis_a = axis_a // gcd_on_list(axis_a)
        B = np.transpose([axis_a, axis_b, axis_c])

        Mtmp = np.linalg.solve(csl, B)
        M = Mtmp @ np.diag([find_smallest_multiplier(Mtmp[:, i]) for i in range(3)])
        transform = csl @ M

        return cls(system, transform=transform, rotation=rotation)


def accommodate_vector_in_lattice(matrix, v):
    """Find integer vector ``x`` s.t. ``matrix @ x = v``."""
    x = np.linalg.solve(matrix, v)
    mult = 1
    if not is_integer_array(x):
        # expand cell
        mult = find_smallest_multiplier(x)
        x = np.around(x * mult)
    return x


def purify_vect(hkl):
    """Round vector to integers."""
    hkl = hkl / np.min(np.abs(hkl[np.nonzero(hkl)]))
    mult = find_smallest_multiplier(hkl)
    hkl = np.around(hkl * mult).astype(int)
    return hkl


def find_csl_parallel_uvw(csl: NDArray, uvw: Axis, atol: float = 1e-8):
    """Make basis of CSL parallel to uvw.

    An obtained basis may expand supercell!
    """

    def _is_zero(array) -> bool:
        return np.allclose(array, 0, atol=atol)

    # substitute the last basis vector with uvw
    c = accommodate_vector_in_lattice(csl, uvw)
    to_parallel = np.zeros((3, 3))
    to_parallel[:, 2] = c
    # choose two vectors from `csl` such that to_parallel is invertible
    for i in range(3):
        if _is_zero(c[i]):
            continue
        to_parallel[(i + 1) % 3, 0] = 1
        to_parallel[(i + 2) % 3, 1] = 1
        break

    csl_uvw = csl @ to_parallel
    assert not _is_zero(np.linalg.det(csl_uvw))

    return csl_uvw


def unique_zcoords(zcoords, decimals=8):
    """Make one-axis coordinates unique."""
    _, indices = np.unique(np.around(zcoords, decimals), return_index=True)
    return zcoords[indices]
