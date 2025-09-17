"""Scalar relaxation and self-consistent SOC mesh and band runs with DFTB+."""

import re
from dataclasses import dataclass
from pathlib import Path
from shutil import copyfile
from subprocess import CalledProcessError

import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.units import Hartree

from .constants import get_xi_ref
from .dftbplus_runner import run_dftbplus


class OutOfDomain(RuntimeError):
    """DFTB+ reported SCC failure or a non-positive-definite overlap."""


class RelaxationFailed(RuntimeError):
    """DFTB+ reached the geometry step limit without convergence."""


@dataclass(frozen=True)
class Settings:
    temperature_ev: float
    scc_tolerance: float
    max_scc_iterations: int
    max_geometry_steps: int
    max_force_ha_bohr: float


def spinorbit_constants(elements, table, shell_order, fitted_shell, scales=None):
    """Scale the selected Jha–Heine values in DFTB+ basis order."""
    scales = {} if scales is None else scales
    return {
        element: tuple(
            get_xi_ref(element, shell, table) * scales.get((element, shell), 1.0)
            if shell == fitted_shell[element] else 0.0
            for shell in shell_order[element]
        )
        for element in dict.fromkeys(elements)
    }


def _spinorbit_block(constants) -> list[str]:
    lines = ["  SpinOrbit = {", "    Dual = Yes"]
    for element, values in constants.items():
        slots = " ".join(f"{value:.16g}" for value in values)
        lines.append(f"    {element} = {{ {slots} }}")
    return [*lines, "  }"]


def input_hsd(
    stage: str, settings: Settings, *, mesh=None, kpoints=None, cubic=False, soc=None, read_charges=False
) -> str:
    """Write HSD for one of the three periodic stages."""
    lines = ['Geometry = GenFormat {', '  <<< "geometry.gen"', '}']
    if stage == "relax":
        constraint = ["  Isotropic = Yes"] if cubic else [
            "  FixAngles = Yes", "  FixLengths = { No No Yes }",
        ]
        lines.extend([
            "Driver = ConjugateGradient {",
            "  LatticeOpt = Yes",
            *constraint,
            f"  MaxSteps = {settings.max_geometry_steps}",
            f"  MaxForceComponent = {settings.max_force_ha_bohr:.16g}",
            '  OutputPrefix = "relaxed"',
            "}",
        ])
    lines.extend([
        "Hamiltonian = xTB {",
        '  Method = "GFN1-xTB"',
        "  SCC = Yes",
        "  ShellResolvedSCC = Yes",
        f"  SCCTolerance = {settings.scc_tolerance:.16g}",
        f"  Filling = Fermi {{ Temperature [eV] = {settings.temperature_ev:.16g} }}",
    ])
    if stage in ("relax", "scc"):
        (nx, ny, nz), (sx, sy, sz) = mesh
        if read_charges:
            lines.append("  ReadInitialCharges = Yes")
        lines.extend([
            f"  MaxSCCIterations = {settings.max_scc_iterations}",
            "  ConvergentSCCOnly = Yes",
            "  KPointsAndWeights = SupercellFolding {",
            f"    {nx} 0 0",
            f"    0 {ny} 0",
            f"    0 0 {nz}",
            f"    {sx} {sy} {sz}",
            "  }",
        ])
    elif stage == "band":
        lines.extend([
            "  ReadInitialCharges = Yes",
            "  MaxSCCIterations = 1",
            "  ConvergentSCCOnly = No",
            "  KPointsAndWeights = {",
        ])
        lines.extend("    " + " ".join(f"{x:.16g}" for x in point) + " 1" for point in kpoints)
        lines.append("  }")
    else:
        raise ValueError(f"Unknown DFTB+ stage {stage!r}")
    if soc is not None:
        if stage == "relax":
            raise ValueError("Geometry relaxation requires a scalar calculation")
        lines.extend(_spinorbit_block(soc))
    return "\n".join([
        *lines, "}", "Options { WriteResultsTag = Yes }", "Analysis { WriteBandOut = Yes }", "",
    ])


def run_stage(run_dir: Path, atoms: Atoms, executable: Path, hsd: str) -> str:
    run_dir.mkdir(parents=True, exist_ok=True)
    periodic = atoms.copy()
    periodic.pbc = True  # Treat slabs as 3D supercells with explicit vacuum.
    write(run_dir / "geometry.gen", periodic)
    (run_dir / "dftb_in.hsd").write_text(hsd)
    try:
        result = run_dftbplus(run_dir, executable)
    except CalledProcessError as exc:
        (run_dir / "stdout.log").write_text(exc.stdout)
        (run_dir / "stderr.log").write_text(exc.stderr)
        for line in exc.stdout.splitlines():
            normalized = line.lower()
            if "scc is not converged" in normalized or "cholesky factorization of b" in normalized:
                raise OutOfDomain(line.strip()) from exc
        raise RuntimeError(f"DFTB+ failed in {run_dir}\n{exc.stdout}\n{exc.stderr}") from exc
    (run_dir / "stdout.log").write_text(result.stdout)
    (run_dir / "stderr.log").write_text(result.stderr)
    return result.stdout


def relax(atoms, run_dir, executable, settings, mesh, *, cubic):
    stdout = run_stage(
        run_dir, atoms, executable, input_hsd("relax", settings, mesh=mesh, cubic=cubic)
    )
    if "Geometry converged" not in stdout:
        raise RelaxationFailed("Geometry did not converge. See " + str(run_dir / "stdout.log"))
    relaxed = read(run_dir / "relaxed.gen")
    relaxed.pbc = atoms.pbc
    initial = np.asarray(atoms.cell)
    final = np.asarray(relaxed.cell)
    scales = relaxed.cell.lengths() / atoms.cell.lengths()
    if cubic:
        scales[:] = scales[0]
    else:
        scales[2] = 1.0
    if not np.allclose(final, initial * scales[:, None], rtol=0.0, atol=1e-7):
        raise RuntimeError(f"DFTB+ violated the lattice constraints in {run_dir}.")
    return relaxed


@dataclass(frozen=True)
class Spectrum:
    electrons: int
    mesh_ev: np.ndarray
    points_ev: np.ndarray


def _energies(run_dir: Path) -> np.ndarray:
    path = run_dir / "results.tag"
    with path.open() as stream:
        for line in stream:
            if not line.startswith("eigenvalues "):
                continue
            header = re.fullmatch(r"eigenvalues\s+:real:3:(\d+),(\d+),1\s*", line)
            if header is None:
                raise ValueError(f"Expected one scalar or spinor channel in {path}: {line.strip()}")
            bands, points = map(int, header.groups())
            values = []
            for row in stream:
                if ":" in row:
                    break
                values.extend(map(float, row.split()))
            return np.asarray(values).reshape(points, bands) * Hartree
    raise ValueError(f"No eigenvalues in {path}")


def spectrum(atoms, run_dir, executable, settings, mesh, kpoints, *, soc=None, initial_charges=None) -> Spectrum:
    mesh_dir = run_dir / "scc"
    band_dir = run_dir / "band"
    if initial_charges is not None:
        mesh_dir.mkdir(parents=True, exist_ok=True)
        copyfile(initial_charges, mesh_dir / "charges.bin")
    run_stage(mesh_dir, atoms, executable, input_hsd(
        "scc", settings, mesh=mesh, soc=soc, read_charges=initial_charges is not None,
    ))
    band_dir.mkdir(parents=True, exist_ok=True)
    copyfile(mesh_dir / "charges.bin", band_dir / "charges.bin")
    stdout = run_stage(
        band_dir, atoms, executable, input_hsd("band", settings, kpoints=kpoints, soc=soc)
    )
    counts = re.findall(r"Nr\. of (?:up |down )?electrons:\s+(\S+)", stdout)
    count = sum(map(float, counts))
    if not counts or not count.is_integer():
        raise ValueError(f"Expected an integer electron count in {band_dir / 'stdout.log'}")
    return Spectrum(int(count), _energies(mesh_dir), _energies(band_dir))
