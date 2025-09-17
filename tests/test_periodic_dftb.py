import re
from subprocess import CalledProcessError

import numpy as np
import pytest
from ase.units import Hartree

from configs import pilot_periodic as pilot
from configs.si_strain_demo import STUDY as si_demo
from ml_soc_xtb import periodic_dftb as dftb
from ml_soc_xtb.constants import default_soc_constants_path, load_soc_constants
from ml_soc_xtb.dftbplus_parser import parse_band_out_file
from ml_soc_xtb.periodic_observables import splitting
from ml_soc_xtb.periodic_structures import cubic_cell, special_kpoints


def test_relax_constraints_and_band_restart_hsd():
    mesh = pilot.K_MESHES["zincblende"]
    cubic = dftb.input_hsd("relax", pilot.DFTB_SETTINGS, mesh=mesh, cubic=True)
    slab = dftb.input_hsd("relax", pilot.DFTB_SETTINGS, mesh=pilot.K_MESHES["ti_bilayer"])
    band = dftb.input_hsd("band", pilot.DFTB_SETTINGS, kpoints=[(0.0, 0.0, 0.0), (0.5, 0.5, 0.5)])
    restart = dftb.input_hsd("scc", pilot.DFTB_SETTINGS, mesh=mesh, read_charges=True)
    assert "Driver = ConjugateGradient" in cubic
    assert "Isotropic = Yes" in cubic
    assert "MaxForceComponent = 0.0001" in cubic
    assert "SupercellFolding" in cubic and "0.5 0.5 0.5" in cubic
    assert "SpinOrbit" not in cubic
    assert "FixAngles = Yes" in slab and "FixLengths = { No No Yes }" in slab
    assert "ReadInitialCharges = Yes" in band
    assert "MaxSCCIterations = 1" in band and "ConvergentSCCOnly = No" in band
    assert "0.5 0.5 0.5 1" in band and "SupercellFolding" not in band
    assert "ReadInitialCharges = Yes" in restart and "ConvergentSCCOnly = Yes" in restart
    assert "Options { WriteResultsTag = Yes }" in band and "WriteResultsTag = Yes" in restart


def test_tag_energies_preserve_kpoint_order_and_precision(tmp_path):
    (tmp_path / "results.tag").write_text(
        "eigenvalues         :real:3:3,2,1\n"
        "   -1.419519338712345E+00       -1.4329871234567E-02\n"
        "    0.1234567890123E-01       -1.150596430617468E+00\n"
        "   -1.816110264482706E-02      3.017163754893200E-01\n"
        "filling             :real:3:3,2,1\n"
        "    2.000000000000000E+00\n"
    )
    expected = [[-1.419519338712345, -1.4329871234567e-2, 0.1234567890123e-1],
                [-1.150596430617468, -1.816110264482706e-2, 3.017163754893200e-1]]
    np.testing.assert_allclose(dftb._energies(tmp_path), np.asarray(expected) * Hartree, rtol=1e-14)


def test_spinorbit_basis_order_and_semicore_omission():
    table = load_soc_constants(default_soc_constants_path())
    constants = dftb.spinorbit_constants(("Mo", "H", "Zn", "S"), table, pilot.SHELL_ORDER, pilot.FITTED_SHELL)
    assert constants["Mo"] == pytest.approx((0.00335420435128631, 0, 0))
    assert constants["H"] == constants["Zn"] == (0, 0)
    assert constants["S"] == pytest.approx((0, 0.00291963168457041, 0))
    scaled = dftb.spinorbit_constants(("Mo", "S"), table, pilot.SHELL_ORDER, pilot.FITTED_SHELL, {("Mo", "d"): 1.2})
    assert scaled["Mo"] == pytest.approx(np.asarray(constants["Mo"]) * 1.2)
    assert scaled["S"] == constants["S"]
    text = dftb.input_hsd("scc", pilot.DFTB_SETTINGS, mesh=pilot.K_MESHES["tmdc_monolayer"], soc=constants)
    assert "Dual = Yes" in text
    slots = re.search(r"Mo = \{ ([^}]+) \}", text)[1]
    assert tuple(map(float, slots.split())) == pytest.approx(constants["Mo"])


@pytest.mark.parametrize("error", [
    "-> SCC is NOT converged, maximal SCC iterations exceeded",
    "-> Failure in LAPACK routine zpotrf, unable to complete Cholesky factorization of B, info: 23",
])
def test_dftb_scientific_failure_preserves_stdout_and_flags_without_retry(tmp_path, monkeypatch, error):
    def runner(run_dir, executable):
        raise CalledProcessError(1, [str(executable)], output=f"ERROR!\n{error}\n", stderr="ERROR STOP")

    monkeypatch.setattr(dftb, "run_dftbplus", runner)
    with pytest.raises(dftb.OutOfDomain) as flag:
        dftb.run_stage(tmp_path, pilot.MATERIALS[0].build(), "dftb+", "")
    assert str(flag.value) == error
    assert (tmp_path / "stdout.log").read_text() == f"ERROR!\n{error}\n"


def test_dftb_input_errors_fail_fast(tmp_path, monkeypatch):
    def runner(run_dir, executable):
        raise CalledProcessError(1, [str(executable)], output="-> Invalid character(s).", stderr="ERROR STOP")

    monkeypatch.setattr(dftb, "run_dftbplus", runner)
    with pytest.raises(RuntimeError, match="Invalid character") as failure:
        dftb.run_stage(tmp_path, pilot.MATERIALS[0].build(), "dftb+", "")
    assert not isinstance(failure.value, dftb.OutOfDomain)


@pytest.mark.dft_integration
def test_periodic_si_fixed_geometry_regression(tmp_path, dftbplus_executable, monkeypatch):
    monkeypatch.setenv("OMP_NUM_THREADS", "1")
    monkeypatch.setenv("OPENBLAS_NUM_THREADS", "1")
    atoms = cubic_cell("Si", "diamond", 5.576263382)
    labels = ("G", "L")
    kpoints = special_kpoints(atoms, labels)
    np.testing.assert_allclose(kpoints["L"], [0.5, 0.5, 0.5])
    args = (dftbplus_executable, si_demo.dftb, si_demo.meshes["zincblende"], kpoints.values())
    table = load_soc_constants(default_soc_constants_path())
    constants = dftb.spinorbit_constants(("Si",), table, si_demo.shell_order, si_demo.fitted_shell)
    scalar = dftb.spectrum(atoms, tmp_path / "scalar", *args)
    soc = dftb.spectrum(atoms, tmp_path / "soc", *args, soc=constants)
    for mode, spectrum in (("scalar", scalar), ("soc", soc)):
        for stage, values in (("scc", spectrum.mesh_ev), ("band", spectrum.points_ev)):
            rounded = parse_band_out_file(tmp_path / mode / stage / "band.out")["energy_ev"].to_numpy()
            # Compare using DFTB+ 24.1's conversion and 0.001 eV print precision.
            dftb_ev = values / Hartree * 27.2113845
            np.testing.assert_allclose(dftb_ev, rounded.reshape(values.shape), atol=0.0005, rtol=0)

    expected = ((0.0330, (3,)), (0.0261, (2,)), (0.0244, (3,)), (0.0088, (1, 2)))
    for obs, (value, pattern) in zip(si_demo.observables["zincblende"], expected, strict=True):
        i = labels.index(obs.kpoint)
        result = splitting(scalar.points_ev[i], soc.points_ev[i], scalar.electrons, obs, si_demo.tolerance_ev)
        assert result.value_ev == pytest.approx(value, abs=0.0002)
        assert result.pattern == pattern
        if obs.name == "l_cb":
            assert result.neighbour_gap_ev == pytest.approx(0.0215, abs=0.0002)
