"""Bounded SOC scale fitting with a Gaussian prior and numerical sensitivities."""

from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares

LOWER, UPPER = 0.5, 1.5


@dataclass(frozen=True)
class Settings:
    tau: float = 0.25
    relative_error: float = 0.1
    fd_step: float = 0.02
    per_sample: bool = True

    def __post_init__(self):
        if not self.tau > 0 or not 0 < self.fd_step <= 0.5:
            raise ValueError("Fit requires a positive prior width and a finite-difference step in (0, 0.5]")
        if not np.isfinite(self.relative_error) or self.relative_error <= 0:
            raise ValueError("Fit requires a finite positive relative model error")


@dataclass(frozen=True)
class FitResult:
    scales: np.ndarray
    stderr: np.ndarray
    identifiable: np.ndarray
    at_bound: np.ndarray
    max_effect_ev: np.ndarray
    chi2_reduced: float | None


def fit_scales(forward, n_parameters, values_ev, sigma_ev, settings, *, baseline=None):
    """Fit scales with an injected forward function returning eligible splittings in eV."""
    values = np.asarray(values_ev, dtype=float)
    sigma = np.asarray(sigma_ev, dtype=float)
    denominator = settings.relative_error * np.maximum(np.abs(values), sigma)
    if values.size == 0 or np.any(denominator <= 0):
        raise ValueError("Fit requires compatible reference targets with positive residual denominators")

    scales = np.ones(n_parameters)
    base = forward(scales) if baseline is None else np.asarray(baseline)
    parameters = np.arange(n_parameters)

    def derivatives(point, prediction, indices):
        columns = np.empty((values.size, len(indices)))
        for column, index in enumerate(indices):
            step = settings.fd_step if point[index] + settings.fd_step <= UPPER else -settings.fd_step
            perturbed = point.copy()
            perturbed[index] += step
            columns[:, column] = (forward(perturbed) - prediction) / step
        return columns

    initial = derivatives(scales, base, parameters)
    effects = np.abs(initial) * (UPPER - LOWER)
    identifiable = np.any(effects > 3 * sigma[:, None], axis=0)
    active = parameters[identifiable]
    stderr = np.full(n_parameters, np.nan)
    at_bound = np.zeros(n_parameters, dtype=bool)
    data_residual = (base - values) / denominator

    if active.size:
        prior_jac = np.eye(len(active)) / settings.tau
        prediction = base

        def residual(x):
            nonlocal prediction
            point = scales.copy()
            point[active] = x
            prediction = base if np.all(x == 1) else forward(point)
            return np.r_[(prediction - values) / denominator, (x - 1) / settings.tau]

        def jacobian(x):
            # TRF calls the Jacobian at the point whose residual it just evaluated.
            point = scales.copy()
            point[active] = x
            jac = initial[:, active] if np.all(x == 1) else derivatives(point, prediction, active)
            return np.vstack((jac / denominator[:, None], prior_jac))

        fit = least_squares(
            residual, scales[active], jac=jacobian, method="trf", bounds=(LOWER, UPPER),
        )
        if not fit.success:
            raise RuntimeError(f"SOC scale fit failed: {fit.message}")
        scales[active] = fit.x
        data_residual = fit.fun[:values.size]
        stderr[active] = np.sqrt(np.diag(np.linalg.inv(fit.jac.T @ fit.jac)))
        at_bound[active] = (np.isclose(fit.x, LOWER, atol=1e-6, rtol=0)
                            | np.isclose(fit.x, UPPER, atol=1e-6, rtol=0))

    dof = values.size - active.size
    chi2_reduced = float(data_residual @ data_residual / dof) if dof > 0 else None
    return FitResult(scales, stderr, identifiable, at_bound, effects.max(axis=0), chi2_reduced)
