"""Provisional GPAW references and a common CSV schema for external targets."""

import csv
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

from .periodic_observables import InvalidManifold, Measurement, splitting


KNOBS = ("ecut", "kmesh", "nbands")
DELTA_FIELDS = tuple(f"delta_{knob}_ev" for knob in KNOBS)
REFERENCE_FIELDS = (
    "sample", "observable", "value_ev", "sigma_ev", *DELTA_FIELDS, "pattern", "neighbour_gap_ev", "source",
)


@dataclass(frozen=True)
class Resolution:
    ecut: float
    kmesh: tuple[int, int, int]
    nbands: int


@dataclass(frozen=True)
class Settings:
    base: Resolution
    tight: Resolution
    temperature_ev: float = 0.025


@dataclass(frozen=True)
class Target:
    measurement: Measurement
    sigma_ev: float
    deltas_ev: tuple[float, ...]
    source: str

    def row(self, sample, observable):
        return dict(
            sample=sample, observable=observable, value_ev=self.measurement.value_ev,
            sigma_ev=self.sigma_ev, pattern=format_pattern(self.measurement.pattern),
            neighbour_gap_ev=self.measurement.neighbour_gap_ev, source=self.source,
            **dict(zip(DELTA_FIELDS, self.deltas_ev, strict=True)),
        )


def format_pattern(pattern):
    return ",".join(map(str, pattern))


def load_reference(path: Path) -> dict[tuple[str, str], Target]:
    """Load uniquely keyed targets with scalar patterns and adjacent gaps."""
    targets = {}
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        missing = set(REFERENCE_FIELDS) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path} is missing reference columns: {', '.join(sorted(missing))}")
        for row in reader:
            key = row["sample"].strip(), row["observable"].strip()
            try:
                pattern = tuple(int(size) for size in row["pattern"].split(","))
                value = float(row["value_ev"])
                sigma = float(row["sigma_ev"])
                gap = float(row["neighbour_gap_ev"])
                deltas = tuple(float(row[field]) for field in DELTA_FIELDS)
                source = row["source"].strip()
            except (AttributeError, TypeError, ValueError) as exc:
                raise ValueError(f"{path}:{reader.line_num} has an invalid reference row for {key}") from exc
            if (not all(key) or not source or min(pattern) < 1 or not np.isfinite([value, sigma, gap, *deltas]).all()
                    or gap <= 0 or min(value, sigma, *deltas) < 0 or max(value, sigma) == 0):
                raise ValueError(f"{path}:{reader.line_num} has an invalid reference row for {key}")
            if key in targets:
                raise ValueError(f"{path}:{reader.line_num} repeats reference target {key}")
            targets[key] = Target(Measurement(value, pattern, gap), sigma, deltas, source)
    if not targets:
        raise ValueError(f"{path} has no reference targets")
    return targets


def gpaw_spectra(atoms, points, resolution, temperature_ev, run_dir):
    """Return scalar and post-SCF spinor eigenvalues at explicit scaled k-points."""
    from gpaw import PW, FermiDirac
    from gpaw.calculator import GPAW
    from gpaw.spinorbit import soc_eigenstates

    run_dir.mkdir(parents=True, exist_ok=True)
    cell = atoms.copy()
    cell.pbc = True
    calc = GPAW(
        mode=PW(resolution.ecut), xc="PBE", kpts=resolution.kmesh,
        occupations=FermiDirac(temperature_ev), spinpol=False,
        txt=str(run_dir / "scf.log"),
    )
    cell.calc = calc
    cell.get_potential_energy()
    bands = calc.fixed_density(
        kpts=list(points.values()), symmetry="off", nbands=resolution.nbands,
        convergence={"bands": resolution.nbands}, txt=str(run_dir / "bands.log"),
    )
    scalar = np.stack([bands.get_eigenvalues(kpt=i) for i in range(len(points))])
    spinor = soc_eigenstates(bands).eigenvalues()
    electrons = bands.get_number_of_electrons()
    if electrons != int(electrons):
        raise ValueError(f"Expected an integer GPAW electron count in {run_dir}")
    np.savez(
        run_dir / "eigenvalues.npz", labels=list(points), kpoints=list(points.values()),
        electrons=int(electrons), scalar_ev=scalar, soc_ev=spinor,
    )
    return int(electrons), scalar, spinor


def converge_reference(atoms, points, observations, settings, tolerance_ev, run_dir):
    """Use the all-tight values and the largest single-knob change from base."""
    from gpaw import __version__

    runs = {
        "base": settings.base,
        **{knob: replace(settings.base, **{knob: getattr(settings.tight, knob)}) for knob in KNOBS},
        "tight": settings.tight,
    }
    indices = {label: i for i, label in enumerate(points)}
    measurements = {}
    for label, resolution in runs.items():
        electrons, scalar, spinor = gpaw_spectra(
            atoms, points, resolution, settings.temperature_ev, run_dir / label,
        )
        measurements[label] = {}
        for obs in observations:
            i = indices[obs.kpoint]
            try:
                measurements[label][obs.name] = splitting(scalar[i], spinor[i], electrons, obs, tolerance_ev)
            except InvalidManifold as exc:
                raise InvalidManifold(f"{run_dir}/{label}, {obs.name}: {exc}") from exc

    targets = {}
    for obs in observations:
        name = obs.name
        tight = measurements["tight"][name]
        if any(run[name].pattern != tight.pattern for run in measurements.values()):
            raise InvalidManifold(f"{run_dir}: scalar pattern changes across convergence runs for {name}")
        base = measurements["base"][name]
        deltas = tuple(abs(measurements[knob][name].value_ev - base.value_ev) for knob in KNOBS)
        targets[name] = Target(tight, max(deltas), deltas, f"GPAW {__version__} PW/PBE post-SCF SOC (provisional)")
    return targets
