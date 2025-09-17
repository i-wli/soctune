"""Geometry descriptors and weighted ridge models for L1 log scales."""

from dataclasses import dataclass

import numpy as np
from ase.neighborlist import neighbor_list

from .periodic_fit import LOWER, UPPER


@dataclass(frozen=True)
class Settings:
    alpha: float = 1.0

    def __post_init__(self):
        if not self.alpha > 0:
            raise ValueError("Ridge alpha must be positive")


def mean_nearest_neighbor(atoms, element):
    """Mean over atoms of one species of the distance to their closest neighbour."""
    # A periodic image along the shortest lattice vector bounds every nearest distance.
    cutoff = atoms.cell.lengths()[atoms.pbc].min() * 1.01
    indices, distances = neighbor_list("id", atoms, cutoff)
    nearest = np.full(len(atoms), np.inf)
    np.minimum.at(nearest, indices, distances)
    selected = np.asarray(atoms.get_chemical_symbols()) == element
    return float(nearest[selected].mean())


@dataclass(frozen=True)
class Ridge:
    mean: np.ndarray
    std: np.ndarray
    intercept: float
    coefficients: np.ndarray

    def predict(self, descriptors):
        x = (np.asarray(descriptors) - self.mean) / self.std
        log_scale = self.intercept + x @ self.coefficients
        return np.exp(np.clip(log_scale, np.log(LOWER), np.log(UPPER)))


def fit_ridge(descriptors, scales, stderr, alpha):
    """Fit log scales weighted by inverse delta-method variance, with a free intercept."""
    descriptors = np.asarray(descriptors)
    scales = np.asarray(scales)
    stderr = np.asarray(stderr)
    mean = descriptors.mean(axis=0)
    std = descriptors.std(axis=0)
    std[std == 0] = 1.0
    x = (descriptors - mean) / std
    design = np.column_stack((np.ones(len(x)), x))
    sqrt_weight = scales / stderr
    weighted = design * sqrt_weight[:, None]
    penalty = np.eye(design.shape[1]) * alpha
    penalty[0, 0] = 0.0
    coef = np.linalg.solve(weighted.T @ weighted + penalty, weighted.T @ (np.log(scales) * sqrt_weight))
    return Ridge(mean, std, float(coef[0]), coef[1:])


def leave_one_out(samples):
    """Yield each held-out sample and its disjoint training set, in input order."""
    for i, held_out in enumerate(samples):
        yield held_out, samples[:i] + samples[i + 1:]
