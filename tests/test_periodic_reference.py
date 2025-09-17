import csv

import pytest

from ml_soc_xtb.periodic_observables import Measurement, Observable, splitting
from ml_soc_xtb.periodic_reference import REFERENCE_FIELDS, Resolution, Target, gpaw_spectra, load_reference
from ml_soc_xtb.periodic_structures import cubic_cell, special_kpoints


@pytest.mark.parametrize("field", ["sigma_ev", "pattern", "neighbour_gap_ev"])
def test_reference_csv_rejects_missing_column(tmp_path, field):
    path = tmp_path / "reference.csv"
    path.write_text(",".join(name for name in REFERENCE_FIELDS if name != field) + "\n")
    with pytest.raises(ValueError, match=field):
        load_reference(path)


@pytest.mark.parametrize("field, value", [
    ("sigma_ev", ""), ("pattern", ""), ("neighbour_gap_ev", ""), ("neighbour_gap_ev", "0"),
])
def test_reference_csv_rejects_invalid_value(tmp_path, field, value):
    path = tmp_path / "reference.csv"
    target = Target(Measurement(0.03, (1, 2), 1.4), 0.001, (0.0, 0.0, 0.001), "test")
    row = target.row("sample", "l_cb")
    row[field] = value
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=REFERENCE_FIELDS)
        writer.writeheader()
        writer.writerow(row)
    with pytest.raises(ValueError, match="invalid reference"):
        load_reference(path)


@pytest.mark.dft_integration
def test_gpaw_si_gamma_valence_smoke(tmp_path):
    atoms = cubic_cell("Si", "diamond", 5.576263382)
    points = special_kpoints(atoms, ("G",))
    count, scalar, spinor = gpaw_spectra(atoms, points, Resolution(300, (4, 4, 4), 12), 0.025, tmp_path)
    obs = Observable("gamma_vb", "G", "valence", 3)
    measured = splitting(scalar[0], spinor[0], count, obs, 0.0005)
    assert measured.pattern == (3,)
    assert measured.value_ev > 0
    assert (tmp_path / "eigenvalues.npz").is_file()
