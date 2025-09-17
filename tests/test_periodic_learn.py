import numpy as np

from ml_soc_xtb.periodic_learn import fit_ridge, leave_one_out


def test_weighted_ridge_leave_one_out_uses_training_descriptors_only():
    x = np.linspace(-1, 1, 9)[:, None]
    exact = np.exp(0.04 + 0.15 * x[:, 0])
    noisy = exact.copy()
    noisy[4] = 1.35
    stderr = np.full(9, 0.002)
    stderr[4] = 0.9
    predictions = []
    for held_out, training in leave_one_out(tuple(range(len(x)))):
        indices = list(training)
        assert held_out not in indices
        model = fit_ridge(x[indices], noisy[indices], stderr[indices], alpha=1.0)
        np.testing.assert_allclose(model.mean, x[indices].mean(axis=0))
        predictions.append(float(model.predict(x[held_out])))
    np.testing.assert_allclose(predictions, exact, atol=0.001)
