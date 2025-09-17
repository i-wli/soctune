"""Primitive cells and reciprocal-space points for the periodic pilot."""

from math import cos, pi, sin, sqrt

import numpy as np
from ase import Atoms
from ase.build import bulk, mx2
from ase.io import read


def cubic_cell(formula: str, structure: str, a: float) -> Atoms:
    return bulk(formula, structure, a=a)


def tmdc_cell(formula: str, a: float, thickness: float, vacuum: float) -> Atoms:
    return mx2(formula, kind="2H", a=a, thickness=thickness, vacuum=vacuum)


def methyl_bilayer(
    element: str, a: float, buckling: float, xc: float, ch: float, vacuum: float
) -> Atoms:
    """Build a buckled honeycomb with inversion-related methyl groups."""
    cell = [(a, 0, 0), (-a / 2, sqrt(3) * a / 2, 0), (0, 0, 0)]
    symbols = []
    positions = []
    radius = ch * sqrt(8 / 9)
    for x, y, side in [(0, 0, 1), (a / 2, sqrt(3) * a / 6, -1)]:
        z = side * buckling / 2
        carbon_z = z + side * xc
        symbols.extend((element, "C", "H", "H", "H"))
        positions.extend(((x, y, z), (x, y, carbon_z)))
        for i in range(3):
            angle = i * 2 * pi / 3
            positions.append(
                (x + side * radius * cos(angle),
                 y + side * radius * sin(angle),
                 carbon_z + side * ch / 3)
            )
    atoms = Atoms(symbols, positions=positions, cell=cell, pbc=(True, True, False))
    atoms.center(vacuum=vacuum, axis=2)
    return atoms


def special_kpoints(atoms: Atoms, labels: tuple[str, ...]) -> dict[str, np.ndarray]:
    """Return ASE special points in the input cell's scaled coordinates."""
    points = atoms.cell.bandpath(pbc=atoms.pbc).special_points
    return {label: points[label] for label in labels}


def strained(atoms: Atoms, eps: float) -> Atoms:
    """Apply hydrostatic strain along the periodic axes at fixed fractional coordinates."""
    sample = atoms.copy()
    cell = np.asarray(atoms.cell).copy()
    cell[atoms.pbc] *= 1 + eps
    sample.set_cell(cell, scale_atoms=True)
    return sample


def read_structure(path) -> Atoms:
    atoms = read(path, index=0)
    if not atoms.pbc.any() or atoms.cell.rank != 3:
        raise ValueError(f"Expected a periodic structure with a complete cell in {path}")
    return atoms
