import numpy as np
import pytest

from ml_soc_xtb.periodic_observables import InvalidManifold, Observable, gap_change, splitting


TOL = 0.0005


@pytest.mark.parametrize(
    "scalar,soc,kind,size,pattern",
    [
        ([-4, -1, -1, -1, 2], [-4, -4, -1.2, -1.2, -1, -1, -1, -1, 2, 2], "valence", 3, (3,)),
        ([-4, -1, -1, -1, 2], [-4, -4, -1, -1, -1, -1, -0.8, -0.8, 2, 2], "valence", 3, (3,)),
        ([-4, -3, -2, -1, 1, 4, 4], [-4, -4, -3, -3, -2, -2, -1, -1, 1, 1, 3.9, 3.9, 4.1, 4.1],
         "conduction", 2, (1, 2)),
        ([-6, -5, -4, -3, 1, 4, 4], [-6, -6, -5, -5, -4, -4, -3, -3, 0.9, 1.1, 4, 4, 4, 4],
         "conduction", 1, (1,)),
        ([-6, -5, -4, -1, 2], [-6, -6, -5, -5, -4, -4, -1.1, -0.9, 2, 2], "valence", 1, (1,)),
        ([-6, -5, -4, 0, 0, 3], [-6, -6, -5, -5, -4, -4, -0.1, -0.1, 0.1, 0.1, 3, 3],
         "dirac", 2, (2,)),
    ],
)
def test_splitting_scalar_manifold_to_spinor_multiplets(scalar, soc, kind, size, pattern):
    observable = Observable("test", "G", kind, size)
    measurement = splitting(scalar, soc, 8, observable, TOL)
    assert measurement.value_ev == pytest.approx(0.2)
    assert measurement.pattern == pattern


@pytest.mark.parametrize("below,above", [(-0.03, 1.0), (-1.0, 0.03)])
def test_splitting_reports_the_closer_scalar_neighbour(below, above):
    scalar = [-2, below, 0, 0, above]
    soc = [-2, -2, below, below, -0.005, -0.005, 0.005, 0.005, above, above]
    measurement = splitting(scalar, soc, 8, Observable("l_vb", "L", "valence", 2), TOL)
    assert measurement.neighbour_gap_ev == pytest.approx(0.03)


def test_conduction_requires_the_first_degenerate_manifold_to_have_the_expected_size():
    scalar = [-4, -3, -2, -1, 1, 2, 2, 3, 3, 3]
    observable = Observable("gamma_cb", "G", "conduction", 3)
    with pytest.raises(InvalidManifold, match="conduction manifold is 2-fold, expected 3-fold"):
        splitting(scalar, np.repeat(scalar, 2), 8, observable, TOL)


def test_split_valence_doublet_must_have_two_degenerate_pairs():
    scalar = [-4, -3, -1, -1, 1]
    soc = [-4, -4, -3, -3, -1.2, -1.1, -1, -1, 1, 1]
    with pytest.raises(InvalidManifold, match=r"SOC multiplicities \[1, 1, 2\]"):
        splitting(scalar, soc, 8, Observable("l_vb", "L", "valence", 2), TOL)


def test_dirac_observable_flags_a_gap_already_present_without_soc():
    scalar = [-4, -3, -2, -0.1, 0.1, 1]
    with pytest.raises(InvalidManifold, match="scalar frontier gap 0.2000 eV"):
        splitting(scalar, np.repeat(scalar, 2), 8, Observable("k_dirac", "K", "dirac", 2), TOL)


def test_gap_diagnostic_uses_indirect_extrema_and_preserves_signed_gap():
    scalar = np.array([[0, 0.9], [1, 1.4]])
    soc = np.array([[-0.1, 0.1, 0.8, 1], [0.8, 1.1, 1.8, 2]])
    assert gap_change(scalar, soc, 2) == pytest.approx(-0.2)
