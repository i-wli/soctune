import numpy as np
import pytest

from ml_soc_xtb.periodic_fit import Settings, fit_scales


def test_fit_recovers_known_scale_with_sub_sigma_noise():
    sensitivities = np.linspace(0.18, 0.22, 64)

    def forward(scales):
        return 0.01 + sensitivities * (scales[0] - 1)

    noise = np.random.default_rng(14).normal(0, 1e-5, 64)
    sigma = np.full(64, 0.0001)
    result = fit_scales(forward, 1, forward([1.2]) + noise, sigma, Settings())
    assert result.identifiable.tolist() == [True]
    assert result.scales[0] == pytest.approx(1.2, abs=0.004)
    assert result.stderr[0] > 0
    assert result.at_bound.tolist() == [False]


def test_fit_leaves_insensitive_scale_at_one():
    def forward(scales):
        assert scales[1] in (1.0, 1.02)
        return np.array([0.04 * scales[0]])

    result = fit_scales(forward, 2, [0.05], [0.001], Settings())
    assert result.identifiable.tolist() == [True, False]
    assert result.scales[1] == 1.0
    assert np.isnan(result.stderr[1])
    assert result.chi2_reduced is None


def test_reduced_chi_squared_uses_data_only_and_subtracts_active_parameters():
    def forward(scales):
        return np.array([0.01 * scales[0], 0.02])

    settings = Settings(relative_error=0.1, tau=0.25)
    values = np.array([0.012, 0.025])
    sigma = np.array([0.013, 0.001])
    result = fit_scales(forward, 1, values, sigma, settings)

    assert not result.identifiable[0]
    assert result.chi2_reduced == pytest.approx(((0.002 / 0.0013) ** 2 + 2**2) / 2)
    # Lower the numerical σ so the first row's derivative exceeds 3σ.
    sigma[0] = 0.001
    result = fit_scales(forward, 1, values, sigma, settings)
    assert result.identifiable[0]
    residuals = (forward(result.scales) - values) / (0.1 * values)
    assert result.chi2_reduced == pytest.approx(residuals @ residuals)


@pytest.mark.parametrize(("target", "bound"), [(0.01, 0.5), (0.20, 1.5)])
def test_fit_flags_a_bound(target, bound):
    def forward(scales):
        return np.full(100, 0.04 * scales[0])

    result = fit_scales(forward, 1, np.full(100, target), np.full(100, 0.001), Settings())
    assert result.scales[0] == pytest.approx(bound, abs=1e-6)
    assert result.at_bound.tolist() == [True]
