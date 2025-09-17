from pathlib import Path

import numpy as np
import pytest
from ase.build import mx2

from ml_soc_xtb.periodic_observables import Measurement
from ml_soc_xtb.periodic_reference import Target
from ml_soc_xtb.periodic_structures import strained
from ml_soc_xtb.periodic_workflow import common_observables, compare_target, load_study, splitting_trends


CONFIG = Path(__file__).resolve().parents[1] / "configs" / "si_strain_demo.py"


def test_config_snapshot_rejects_changed_bytes(tmp_path):
    config = tmp_path / "config.py"
    config.write_bytes(CONFIG.read_bytes())
    root = tmp_path / "out"
    load_study(config, root)
    assert (root / "config.py").read_bytes() == config.read_bytes()
    config.write_bytes(config.read_bytes() + b"\n# edited\n")
    with pytest.raises(ValueError, match="Config differs"):
        load_study(config, root)


def test_pattern_mismatch_excludes_observable_from_residual():
    row = dict(status="ok", value_ev=0.009, pattern="2", neighbour_gap_ev=0.02)
    target = Target(Measurement(0.015, (1, 2), 2.0), 0.001, (0.0, 0.0, 0.001), "test")
    compare_target(row, target)
    assert row["status"] == "manifold_mismatch"
    assert row["reference_ev"] == 0.015
    assert "residual_ev" not in row


@pytest.mark.parametrize("gap, status", [(2.5, "ok"), (0.8, "unisolated_manifold")])
def test_isolation_mismatch_uses_reference_splitting_for_both_ratios(gap, status):
    target = Target(Measurement(0.2, (2,), 1.8), 0.001, (0.0, 0.0, 0.001), "test")
    for value in (0.02, 0.35):
        row = dict(status="ok", value_ev=value, pattern="2", neighbour_gap_ev=gap)
        compare_target(row, target)
        assert row["status"] == status
        if status == "ok":
            assert row["residual_ev"] == pytest.approx(value - 0.2)
        else:
            assert "xTB 0.250, reference 0.111" in row["error"]
            assert "residual_ev" not in row


def test_fixed_observables_and_trends_use_every_sample_of_each_material():
    materials = {"s0": "Si", "s1": "Si", "g0": "Ge", "g1": "Ge"}
    descriptors = [
        dict(sample=name, material="Si", element="Si", shell="p", distance_angstrom=distance)
        for name, distance in (("s1", 2.5), ("s0", 2.0))
    ]
    eligible = [
        dict(sample=name, material=materials[name], observable=obs)
        for name, obs in (("s0", "vb"), ("s1", "vb"), ("s1", "cb"), ("g0", "vb"), ("g1", "cb"))
    ]
    assert common_observables(eligible, materials) == {("Si", "vb")}
    loo = []
    for row in eligible:
        x = int(row["sample"].endswith("1"))
        for level, sign in (("L-1", 1), ("L0", 2), ("L0 common", 3), ("L2", -3)):
            loo.append(dict(row, level=level, value_ev=0.02 + sign * x / 1000, reference_ev=0.02 - x / 1000))
    assert list(splitting_trends(loo, descriptors, materials)) == [
        ("Si", "Si p", "vb", "-2.000", "2.000", "4.000", "6.000", "-6.000"),
    ]


def test_strain_only_scales_periodic_lattice_vectors():
    atoms = mx2("MoS2", vacuum=10)
    result = strained(atoms, 0.03)
    np.testing.assert_allclose(result.cell[:2], atoms.cell[:2] * 1.03)
    np.testing.assert_allclose(result.cell[2], atoms.cell[2])
    np.testing.assert_allclose(result.get_scaled_positions(), atoms.get_scaled_positions())
