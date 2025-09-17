"""Band splittings from ordered scalar and spinor eigenvalues in eV."""

from dataclasses import dataclass

import numpy as np


class InvalidManifold(ValueError):
    """The requested scalar manifold or SOC multiplicity is absent."""

    def __init__(self, message, pattern=()):
        super().__init__(message)
        self.pattern = pattern


@dataclass(frozen=True)
class Observable:
    name: str
    kpoint: str
    kind: str
    spatial_size: int


@dataclass(frozen=True)
class Measurement:
    value_ev: float
    pattern: tuple[int, ...]
    neighbour_gap_ev: float


def manifold_pattern(groups, occupied, manifold):
    """List scalar multiplicities from the gap outwards to the target."""
    start, stop = manifold
    if start >= occupied:
        return tuple(b - a for a, b in groups if occupied <= a and b <= stop)
    return (stop - start,)


def _manifolds(energies: np.ndarray, tol: float) -> list[tuple[int, int]]:
    """Group levels whose full energy span is within the tolerance."""
    groups = []
    start = 0
    for stop in range(1, len(energies)):
        if energies[stop] - energies[start] > tol:
            groups.append((start, stop))
            start = stop
    return [*groups, (start, len(energies))]


def splitting(
    scalar_ev: np.ndarray,
    soc_ev: np.ndarray,
    electrons: int,
    observable: Observable,
    tolerance_ev: float,
) -> Measurement:
    """Measure the splitting and identify the target from the scalar gap."""
    scalar_ev = np.asarray(scalar_ev)
    soc_ev = np.asarray(soc_ev)
    occupied = electrons // 2
    if (scalar_ev.ndim != 1 or soc_ev.shape != (2 * len(scalar_ev),)
            or electrons % 2 or not 0 < occupied < len(scalar_ev)):
        raise ValueError("Expected one scalar spectrum, its doubled spinor spectrum, and an even electron count")
    if any(not np.isfinite(e).all() or np.any(np.diff(e) < 0) for e in (scalar_ev, soc_ev)):
        raise ValueError("Eigenvalues must be finite and sorted in ascending order")

    groups = _manifolds(scalar_ev, tolerance_ev)
    size = observable.spatial_size
    if observable.kind == "valence":
        manifold = next(((a, b) for a, b in groups if b == occupied), None)
    elif observable.kind == "conduction":
        manifold = next(
            ((a, b) for a, b in groups if a >= occupied and (size == 1 or b - a > 1)),
            None,
        )
        if any(a < occupied < b for a, b in groups):
            manifold = None
    elif observable.kind == "dirac":
        manifold = next(((a, b) for a, b in groups if a < occupied < b), None)
    else:
        raise ValueError(f"Unknown observable kind {observable.kind!r}")

    pattern = manifold_pattern(groups, occupied, manifold) if manifold is not None else ()
    if manifold is None or manifold[1] - manifold[0] != size:
        actual = "absent" if manifold is None else f"{manifold[1] - manifold[0]}-fold"
        reason = f"scalar {observable.kind} manifold is {actual}, expected {size}-fold"
        if observable.kind == "dirac":
            gap = scalar_ev[occupied] - scalar_ev[occupied - 1]
            reason += f" (scalar frontier gap {gap:.4f} eV)"
        raise InvalidManifold(reason, pattern)

    start, stop = (2 * i for i in manifold)
    levels = soc_ev[start:stop]
    blocks = _manifolds(levels, tolerance_ev)
    multiplicities = sorted(b - a for a, b in blocks)
    expected = {1: [1, 1], 2: [2, 2], 3: [2, 4]}[size]
    if multiplicities != expected:
        raise InvalidManifold(
            f"SOC multiplicities {multiplicities}, expected {expected} "
            f"(energy span {levels[-1] - levels[0]:.4f} eV)", pattern,
        )
    first, second = (levels[a:b].mean() for a, b in blocks)
    a, b = manifold
    below = scalar_ev[a] - scalar_ev[a - 1] if a else np.inf
    above = scalar_ev[b] - scalar_ev[b - 1] if b < len(scalar_ev) else np.inf
    return Measurement(float(second - first), pattern, float(min(below, above)))


def gap_change(scalar_ev: np.ndarray, soc_ev: np.ndarray, electrons: int) -> float:
    """Return the change in the signed band gap on the sampled k-points."""
    nocc = electrons // 2
    scalar_gap = scalar_ev[:, nocc].min() - scalar_ev[:, nocc - 1].max()
    soc_gap = soc_ev[:, electrons].min() - soc_ev[:, electrons - 1].max()
    return float(soc_gap - scalar_gap)
