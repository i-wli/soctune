"""Run a config-driven periodic SOC study."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from dataclasses import dataclass
from itertools import count
from pathlib import Path
from threading import Lock

import numpy as np
from ase import Atoms
from ase.io import write

from . import periodic_dftb as dftb
from . import periodic_fit as fitting
from . import periodic_learn as learning
from . import periodic_reference as reference
from .constants import default_soc_constants_path, load_soc_constants
from .dftbplus_runner import resolve_dftbplus_path
from .periodic_observables import InvalidManifold, Observable, gap_change, splitting
from .periodic_structures import cubic_cell, methyl_bilayer, read_structure, special_kpoints, strained, tmdc_cell


BUILDERS = {"zincblende": cubic_cell, "tmdc_monolayer": tmdc_cell, "ti_bilayer": methyl_bilayer}
ISOLATION_MISMATCH = 0.1
LOO_LEVELS = ("L-1", "L0", "L0 common", "L2")
SAMPLE_FIELDS = ("sample", "material", "family", "strain", "status", "geometry", "error")
SPECTRUM_FIELDS = (
    "sample", "material", "family", "strain", "observable", "kind", "status", "value_ev", "pattern", "neighbour_gap_ev",
    "reference_ev", "sigma_ev", "residual_ev", "relative_residual", "delta_eg_ev", "error",
)
PARAMETER_FIELDS = (
    "level", "sample", "element", "shell", "scale", "stderr", "identifiable", "at_bound", "max_effect_ev",
    "chi2_reduced",
)
DESCRIPTOR_FIELDS = ("sample", "material", "element", "shell", "distance_angstrom")
RIDGE_FIELDS = ("material", "element", "shell", "mean_angstrom", "std_angstrom", "intercept", "coefficient")


@dataclass(frozen=True)
class Material:
    id: str
    family: str
    builder_args: tuple | None = None
    strains: tuple[float, ...] = (0.0,)
    structure_path: Path | None = None

    def __post_init__(self):
        if (self.builder_args is None) == (self.structure_path is None):
            raise ValueError(f"{self.id} requires exactly one of builder_args or structure_path")
        if not self.strains or len(set(self.strains)) != len(self.strains):
            raise ValueError(f"{self.id} requires a non-empty set of distinct strains")
        if any(not np.isfinite(eps) or eps <= -1 for eps in self.strains):
            raise ValueError(f"{self.id} requires finite strains greater than -1")

    def build(self):
        return BUILDERS[self.family](*self.builder_args)


@dataclass(frozen=True)
class Study:
    materials: tuple[Material, ...]
    observables: dict[str, tuple[Observable, ...]]
    dftb: dftb.Settings
    executable: Path
    meshes: dict
    shell_order: dict[str, tuple[str, ...]]
    fitted_shell: dict[str, str | None]
    tolerance_ev: float
    reference: reference.Settings | Path | None = None
    fit: fitting.Settings | None = None
    learn: learning.Settings | None = None
    preserve_cubic_basis: bool = True

    def __post_init__(self):
        if not self.materials or len({m.id for m in self.materials}) != len(self.materials):
            raise ValueError("Study requires distinct material IDs")
        if self.fit is not None and self.reference is None:
            raise ValueError("Fitting requires reference settings or an external reference CSV")
        if self.learn is not None and (self.fit is None or not self.fit.per_sample):
            raise ValueError("Learning requires per-sample fits")


@dataclass(frozen=True)
class Sample:
    row: dict
    observations: tuple[Observable, ...]
    atoms: Atoms | None
    points: dict


def load_study(config: Path, root: Path) -> Study:
    """Load STUDY from the specified Python source and record its exact bytes."""
    source = config.read_bytes()
    snapshot = root / "config.py"
    if snapshot.exists() and snapshot.read_bytes() != source:
        raise ValueError(f"Config differs from the snapshot in {snapshot}")
    spec = importlib.util.spec_from_file_location("_periodic_study_config", config)
    if spec is None:
        raise ValueError(f"Expected a Python configuration file: {config}")
    module = importlib.util.module_from_spec(spec)
    exec(compile(source, config, "exec"), module.__dict__)
    study = module.STUDY
    if not isinstance(study, Study):
        raise TypeError(f"{config} must export STUDY as a Study")
    root.mkdir(parents=True, exist_ok=True)
    if not snapshot.exists():
        snapshot.write_bytes(source)
    return study


def write_csv(path, fields, rows, *, append=False):
    header = not append or not path.exists()
    with path.open("a" if append else "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        if header:
            writer.writeheader()
        writer.writerows(rows)


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def sample_id(material, eps):
    return f"{material.id}_{eps:+g}" if len(material.strains) > 1 else material.id


def relax_stage(study, root, config_dir, executable, jobs):
    def relax_material(material):
        run_dir = root / "relax" / material.id
        status, error = "ok", ""
        try:
            if material.structure_path is not None:
                parent = read_structure(config_dir / material.structure_path)
            else:
                initial = material.build()
                cubic = material.family == "zincblende"
                parent = dftb.relax(
                    initial, run_dir, executable, study.dftb, study.meshes[material.family], cubic=cubic,
                )
                if cubic and study.preserve_cubic_basis:
                    # Diamond and zincblende have a symmetry-fixed fractional basis.
                    parent.set_scaled_positions(initial.get_scaled_positions(wrap=False))
        except dftb.OutOfDomain as exc:
            status, error = "out_of_domain", str(exc)
        except dftb.RelaxationFailed as exc:
            status, error = "relaxation_failed", str(exc)

        rows = []
        for eps in material.strains:
            name = sample_id(material, eps)
            relative = Path("samples") / name / "geometry.extxyz"
            row = dict(sample=name, material=material.id, family=material.family, strain=eps,
                       status=status, error=error, geometry="")
            if status == "ok":
                path = root / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                write(path, strained(parent, eps))
                row["geometry"] = str(relative)
            rows.append(row)
        return rows

    with ThreadPoolExecutor(max_workers=jobs) as pool:
        rows = [row for subset in pool.map(relax_material, study.materials) for row in subset]
    write_csv(root / "samples.csv", SAMPLE_FIELDS, rows)


def read_samples(study, root):
    samples = {}
    for row in read_csv(root / "samples.csv"):
        observations = study.observables[row["family"]]
        labels = tuple(dict.fromkeys(obs.kpoint for obs in observations))
        atoms, points = None, {}
        if row["status"] == "ok":
            atoms = read_structure(root / row["geometry"])
            try:
                points = special_kpoints(atoms, labels)
            except KeyError as exc:
                row.update(status="invalid_cell", error=f"Cell has no ASE special point {exc}")
        samples[row["sample"]] = Sample(row, observations, atoms, points)
    return samples


def reference_stage(study, root, config_dir, jobs):
    samples = read_samples(study, root)
    expected = {
        (name, obs.name)
        for name, sample in samples.items() if sample.row["status"] == "ok"
        for obs in sample.observations
    }
    if isinstance(study.reference, Path):
        source = config_dir / study.reference
        targets = reference.load_reference(source)
        if targets.keys() != expected:
            raise ValueError(f"External reference missing {sorted(expected - targets.keys())} "
                             f"and has unexpected targets {sorted(targets.keys() - expected)}")
        (root / "reference.csv").write_bytes(source.read_bytes())
        return

    with ProcessPoolExecutor(max_workers=jobs) as pool:
        futures = {
            name: pool.submit(
                reference.converge_reference, sample.atoms, sample.points, sample.observations,
                study.reference, study.tolerance_ev, root / "reference" / name,
            )
            for name, sample in samples.items() if sample.row["status"] == "ok"
        }
        rows = []
        for sample, future in futures.items():
            targets = future.result()
            rows.extend(target.row(sample, obs) for obs, target in targets.items())
    write_csv(root / "reference.csv", reference.REFERENCE_FIELDS, rows)


def measurements(study, sample, scalar, soc):
    if scalar.electrons != soc.electrons:
        raise ValueError(f"{sample.row['sample']}: scalar and SOC electron counts differ")
    delta_gap = gap_change(scalar.mesh_ev, soc.mesh_ev, scalar.electrons)
    indices = {label: i for i, label in enumerate(sample.points)}
    rows = []
    for obs in sample.observations:
        i = indices[obs.kpoint]
        row = dict(delta_eg_ev=delta_gap)
        try:
            value = splitting(scalar.points_ev[i], soc.points_ev[i], scalar.electrons, obs, study.tolerance_ev)
        except InvalidManifold as exc:
            row.update(status="invalid_manifold", error=str(exc), pattern=reference.format_pattern(exc.pattern))
        else:
            row.update(status="ok", value_ev=value.value_ev, pattern=reference.format_pattern(value.pattern),
                       neighbour_gap_ev=value.neighbour_gap_ev)
        rows.append(row)
    return rows


def compare_target(row, target):
    if target is None:
        return
    measurement = target.measurement
    row.update(reference_ev=measurement.value_ev, sigma_ev=target.sigma_ev)
    expected = reference.format_pattern(measurement.pattern)
    if row.get("pattern") and row["pattern"] != expected:
        row.update(status="manifold_mismatch", error=f"xTB pattern {row['pattern']}; reference pattern {expected}")
    if row["status"] != "ok":
        return
    ref_value = abs(measurement.value_ev)
    xtb_ratio = ref_value / row["neighbour_gap_ev"]
    ref_ratio = ref_value / measurement.neighbour_gap_ev
    if abs(xtb_ratio - ref_ratio) > ISOLATION_MISMATCH:
        row.update(status="unisolated_manifold", error=f"Reference splitting / adjacent scalar gap: "
                   f"xTB {xtb_ratio:.3f}, reference {ref_ratio:.3f}")
        return
    residual = row["value_ev"] - measurement.value_ev
    row.update(residual_ev=residual, relative_residual=residual / max(ref_value, target.sigma_ev))


def result_rows(sample, results, targets):
    sample_fields = {key: sample.row[key] for key in SAMPLE_FIELDS[:4]}
    rows = []
    for obs, result in zip(sample.observations, results, strict=True):
        row = dict(**sample_fields, observable=obs.name, kind=obs.kind, **result)
        key = row["sample"], obs.name
        compare_target(row, targets.get(key))
        rows.append(row)
    return rows


class Forward:
    """Run SOC at a fixed sample and use its L-1 scalar spectrum."""

    def __init__(self, study, root, executable, workers):
        self.study = study
        self.root = root
        self.executable = executable
        self.workers = workers
        self.samples = read_samples(study, root)
        self.constants = load_soc_constants(default_soc_constants_path())
        self.targets = {}
        if study.reference is not None:
            self.targets = reference.load_reference(root / "reference.csv")

    def spectrum(self, sample, directory, soc=None, initial_charges=None):
        return dftb.spectrum(
            sample.atoms, directory, self.executable, self.study.dftb,
            self.study.meshes[sample.row["family"]], sample.points.values(), soc=soc,
            initial_charges=initial_charges,
        )

    def soc(self, sample, directory, scales, initial_charges=None):
        constants = dftb.spinorbit_constants(
            sample.atoms.get_chemical_symbols(), self.constants,
            self.study.shell_order, self.study.fitted_shell, scales,
        )
        return self.spectrum(sample, directory, soc=constants, initial_charges=initial_charges)

    def baseline(self):
        def run(name):
            sample = self.samples[name]
            status, error = sample.row["status"], sample.row["error"]
            if status == "ok":
                directory = self.root / "baseline" / name
                try:
                    scalar = self.spectrum(sample, directory / "scalar")
                    np.savez(directory / "scalar.npz", electrons=scalar.electrons,
                             mesh_ev=scalar.mesh_ev, points_ev=scalar.points_ev)
                    soc = self.soc(sample, directory / "soc", {})
                except dftb.OutOfDomain as exc:
                    status, error = "out_of_domain", str(exc)
                else:
                    return result_rows(sample, measurements(self.study, sample, scalar, soc), self.targets)
            results = [dict(status=status, error=error) for _ in sample.observations]
            return result_rows(sample, results, self.targets)

        rows = [row for result in self.workers.map(run, self.samples) for row in result]
        write_csv(self.root / "baseline.csv", SPECTRUM_FIELDS, rows)

    def predict(self, name, directory, scales):
        sample = self.samples[name]
        baseline = self.root / "baseline" / name
        with np.load(baseline / "scalar.npz") as data:
            scalar = dftb.Spectrum(int(data["electrons"]), data["mesh_ev"], data["points_ev"])
        soc = self.soc(sample, directory, scales, baseline / "soc" / "scc" / "charges.bin")
        return result_rows(sample, measurements(self.study, sample, scalar, soc), self.targets)

    def parameters(self, names):
        elements = dict.fromkeys(
            element for name in names for element in self.samples[name].atoms.get_chemical_symbols()
        )
        return tuple((element, self.study.fitted_shell[element]) for element in elements
                     if self.study.fitted_shell[element] is not None)


class Evaluation:
    """Evaluate one fit or prediction and retain its per-sample DFTB+ outputs."""

    def __init__(self, model, level):
        self.model = model
        self.directory = model.root / "fit" / level
        self.directory.mkdir(parents=True, exist_ok=True)
        last = max((int(path.name.removeprefix("eval_")) for path in self.directory.glob("eval_*")), default=0)
        self.serial = count(last + 1)
        self.lock = Lock()
        self.scale_fields = {key: "_".join(("s", *key)) for key in model.parameters(
            name for name, sample in model.samples.items() if sample.atoms is not None
        )}

    def run(self, scales_by_sample):
        with self.lock:
            run_id = f"eval_{next(self.serial)}"
        run_dir = self.directory / run_id
        run_dir.mkdir()

        def predict(item):
            name, scales = item
            theta = {field: scales.get(key, 1.0) for key, field in self.scale_fields.items()}
            rows = self.model.predict(name, run_dir / name, scales)
            return [dict(evaluation=run_id, **theta, **row) for row in rows]

        rows = []
        fields = ("evaluation", *self.scale_fields.values(), *SPECTRUM_FIELDS)
        for subset in self.model.workers.map(predict, scales_by_sample.items()):
            with self.lock:
                write_csv(self.directory / "evaluations.csv", fields, subset, append=True)
            rows.extend({field: row.get(field, "") for field in SPECTRUM_FIELDS} for row in subset)
        return rows


def eligible_rows(root):
    return [row for row in read_csv(root / "baseline.csv") if row["status"] == "ok"]


def train(evaluation, training_rows):
    names = tuple(dict.fromkeys(row["sample"] for row in training_rows))
    parameters = evaluation.model.parameters(names)

    def forward(theta):
        scales = dict(zip(parameters, theta, strict=True))
        rows = evaluation.run(dict.fromkeys(names, scales))
        predictions = {(row["sample"], row["observable"]): row for row in rows}
        values = []
        for row in training_rows:
            key = row["sample"], row["observable"]
            prediction = predictions[key]
            if prediction["status"] != "ok":
                raise InvalidManifold(f"{evaluation.directory}, {key}: {prediction['status']} ({prediction['error']})")
            values.append(prediction["value_ev"])
        return np.asarray(values)

    result = fitting.fit_scales(
        forward, len(parameters),
        [float(row["reference_ev"]) for row in training_rows],
        [float(row["sigma_ev"]) for row in training_rows], evaluation.model.study.fit,
        baseline=[float(row["value_ev"]) for row in training_rows],
    )
    return parameters, result


def fitted_parameters(level, sample, parameters, result):
    for i, (element, shell) in enumerate(parameters):
        yield dict(
            level=level, sample=sample, element=element, shell=shell, scale=result.scales[i],
            stderr=result.stderr[i] if result.identifiable[i] else "",
            identifiable=bool(result.identifiable[i]), at_bound=bool(result.at_bound[i]),
            max_effect_ev=result.max_effect_ev[i], chi2_reduced=result.chi2_reduced,
        )


def common_observables(rows, sample_materials):
    """Return material/observable pairs eligible for every sample of that material."""
    material_counts = Counter(sample_materials.values())
    matching = Counter((row["material"], row["observable"]) for row in rows)
    return {pair for pair, count in matching.items() if count == material_counts[pair[0]]}


def fit_stage(model, jobs):
    baseline = eligible_rows(model.root)
    names = tuple(dict.fromkeys(row["sample"] for row in baseline))

    def fit_level(task):
        label, evaluation, training_rows, predicted_names, sample = task
        parameters, result = train(evaluation, training_rows)
        scales = dict(zip(parameters, result.scales, strict=True))
        predictions = [dict(level=label, **row) for row in evaluation.run(dict.fromkeys(predicted_names, scales))]
        return predictions, list(fitted_parameters(label, sample, parameters, result))

    tasks = [("L0", Evaluation(model, "L0"), baseline, names, "")]
    if model.study.fit.per_sample:
        materials = {name: sample.row["material"] for name, sample in model.samples.items()}
        common = common_observables(baseline, materials)
        evaluator = Evaluation(model, "L1")
        for name in names:
            rows = [row for row in baseline
                    if row["sample"] == name and (row["material"], row["observable"]) in common]
            if not rows:
                raise ValueError(f"{name}: L1 requires an observable eligible across all {materials[name]} samples")
            tasks.append(("L1", evaluator, rows, (name,), name))
    valence = [row for row in baseline if row["kind"] == "valence"]
    if valence:
        tasks.append(("L0 valence-only", Evaluation(model, "valence_control"), valence, names, ""))

    predictions, results = [], []
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        for rows, parameters in pool.map(fit_level, tasks):
            predictions.extend(rows)
            results.extend(parameters)
    write_csv(model.root / "parameters.csv", PARAMETER_FIELDS, results)
    write_csv(model.root / "fit.csv", ("level", *SPECTRUM_FIELDS), predictions)


def descriptor_rows(model, params):
    for row in params:
        if row["level"] == "L1":
            sample = model.samples[row["sample"]]
            yield {key: row[key] for key in ("sample", "element", "shell")} | dict(
                material=sample.row["material"], distance_angstrom=learning.mean_nearest_neighbor(
                    sample.atoms, row["element"],
                ),
            )


def ridge_for(labels, descriptors, alpha):
    distances = [[descriptors[(r["sample"], r["element"], r["shell"])]] for r in labels]
    scales = [float(r["scale"]) for r in labels]
    errors = [float(r["stderr"]) for r in labels]
    return learning.fit_ridge(distances, scales, errors, alpha)


def learn_stage(model, jobs):
    root = model.root
    baseline = eligible_rows(root)
    names = tuple(dict.fromkeys(row["sample"] for row in baseline))
    materials = {name: sample.row["material"] for name, sample in model.samples.items()}
    counts = Counter(materials[name] for name in names)
    if not counts or min(counts.values()) < 2:
        raise ValueError("L2 leave-one-out evaluation requires at least two eligible samples per material")
    common = common_observables(baseline, materials)
    params = read_csv(root / "parameters.csv")
    inputs = list(descriptor_rows(model, params))
    descriptors = {(r["sample"], r["element"], r["shell"]): r["distance_angstrom"] for r in inputs}
    groups = {}
    for row in params:
        if row["level"] == "L1" and row["identifiable"] == "True":
            key = materials[row["sample"]], row["element"], row["shell"]
            groups.setdefault(key, []).append(row)

    def fit_fold(split):
        held_out, training = split
        held_rows = [row for row in baseline if row["sample"] == held_out]
        held_keys = {(row["sample"], row["observable"]) for row in held_rows}
        predictions = [dict(level="L-1", **row) for row in held_rows]
        train_rows = [row for row in baseline if row["sample"] in training]
        common_rows = [row for row in train_rows if (row["material"], row["observable"]) in common]
        fits, fold_parameters = [], []
        for level, directory, rows in (("L0", "L0", train_rows), ("L0 common", "L0_common", common_rows)):
            evaluation = Evaluation(model, f"loo/{held_out}/{directory}")
            parameters, fit = train(evaluation, rows)
            fits.append((level, evaluation, dict(zip(parameters, fit.scales, strict=True))))
            fold_parameters.extend(fitted_parameters(level, held_out, parameters, fit))
        l2 = {}
        for element, shell in model.parameters((held_out,)):
            group = materials[held_out], element, shell
            matching = [row for row in groups.get(group, ()) if row["sample"] in training]
            if matching:
                ridge = ridge_for(matching, descriptors, model.study.learn.alpha)
                x = descriptors[(held_out, element, shell)]
                l2[(element, shell)] = float(ridge.predict([x]))
            else:
                l2[(element, shell)] = 1.0
        fits.append(("L2", Evaluation(model, f"loo/{held_out}/L2"), l2))
        for level, runner, scales in fits:
            rows = runner.run({held_out: scales})
            for row in rows:
                key = row["sample"], row["observable"]
                if key in held_keys and row["status"] != "ok":
                    raise InvalidManifold(f"Held-out {level}, {key}: {row['status']} ({row['error']})")
                if key in held_keys:
                    predictions.append(dict(level=level, **row))
        return predictions, fold_parameters

    loo, loo_parameters = [], []
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        for rows, parameters in pool.map(fit_fold, learning.leave_one_out(names)):
            loo.extend(rows)
            loo_parameters.extend(parameters)
    models = []
    for (material, element, shell), labels in groups.items():
        ridge = ridge_for(labels, descriptors, model.study.learn.alpha)
        models.append(dict(
            material=material, element=element, shell=shell, mean_angstrom=ridge.mean[0],
            std_angstrom=ridge.std[0], intercept=ridge.intercept, coefficient=ridge.coefficients[0],
        ))
    write_csv(root / "descriptors.csv", DESCRIPTOR_FIELDS, inputs)
    write_csv(root / "ridge.csv", RIDGE_FIELDS, models)
    write_csv(root / "loo.csv", ("level", *SPECTRUM_FIELDS), loo)
    write_csv(root / "loo_l0_parameters.csv", PARAMETER_FIELDS, loo_parameters)


def table(lines, headers, rows):
    lines.extend(("", "| " + " | ".join(headers) + " |", "| " + " | ".join("---" for _ in headers) + " |"))
    lines.extend("| " + " | ".join(map(str, row)) + " |" for row in rows)


def mev(value):
    return f"{1000 * float(value):.3f}" if value != "" else "—"


def report_plot(root, params):
    from matplotlib.figure import Figure

    rows = read_csv(root / "descriptors.csv")
    descriptors = {(r["sample"], r["element"], r["shell"]): float(r["distance_angstrom"]) for r in rows}
    materials = {r["sample"]: r["material"] for r in rows}
    labels = [r for r in params if r["level"] == "L1" and r["identifiable"] == "True"]
    models = read_csv(root / "ridge.csv")
    if not models:
        return
    figure = Figure(figsize=(6, 3.8 * len(models)), layout="constrained")
    axes = figure.subplots(len(models), 1, squeeze=False)
    for ax, row in zip(axes[:, 0], models, strict=True):
        material, element, shell = row["material"], row["element"], row["shell"]
        selected = [r for r in labels if materials[r["sample"]] == material
                    and (r["element"], r["shell"]) == (element, shell)]
        global_fit = next(r for r in params if r["level"] == "L0" and (r["element"], r["shell"]) == (element, shell))
        distance = [descriptors[(r["sample"], element, shell)] for r in selected]
        ridge = learning.Ridge(
            np.array([float(row["mean_angstrom"])]), np.array([float(row["std_angstrom"])]),
            float(row["intercept"]), np.array([float(row["coefficient"])]),
        )
        x = np.linspace(min(distance), max(distance), 100)
        ax.errorbar(distance, [float(r["scale"]) for r in selected],
                    yerr=[float(r["stderr"]) for r in selected], fmt="o", color="#174b71",
                    capsize=3, label="L1 ± 1 SE")
        ax.axhline(float(global_fit["scale"]), color="#515763", linestyle="--", label="L0")
        ax.plot(x, ridge.predict(x[:, None]), color="#b7502b", label="L2")
        title = f"{element} {shell}" if material == element else f"{material} · {element} {shell}"
        ax.set(xlabel="Mean nearest-neighbour distance (Å)", ylabel="SOC scale", title=title)
        ax.grid(axis="y", alpha=0.25)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncols=3, frameon=False)
    figure.savefig(root / "scales.png", dpi=160)


def splitting_trends(loo, descriptors, sample_materials):
    """Fit OLS slopes of reference and held-out splittings against distance, in meV/Å."""
    common = common_observables([row for row in loo if row["level"] == "L-1"], sample_materials)
    spectra = {(r["sample"], r["observable"], r["level"]): r for r in loo}
    distances = {}
    for row in descriptors:
        group = row["material"], row["element"], row["shell"]
        distances.setdefault(group, {})[row["sample"]] = float(row["distance_angstrom"])
    for (material, element, shell), by_sample in distances.items():
        x = np.array(list(by_sample.values()))
        centered = x - x.mean()
        total = centered @ centered
        if total == 0:
            continue
        for obs in sorted(obs for mat, obs in common if mat == material):
            values = [[float(spectra[(name, obs, "L-1")]["reference_ev"]) for name in by_sample]]
            values.extend([float(spectra[(name, obs, level)]["value_ev"]) for name in by_sample]
                          for level in LOO_LEVELS)
            slopes = np.asarray(values) @ centered / total
            yield material, f"{element} {shell}", obs, *(mev(slope) for slope in slopes)


def report_stage(study, root):
    baseline = read_csv(root / "baseline.csv")
    lines = ["# Periodic SOC study"]
    if study.fit is not None:
        lines.extend(("", "Scales multiply the shell-resolved Jha–Heine constants."))
    if study.reference is None:
        lines.extend(("", "## L-1 splittings"))
        table(lines, ("Sample", "Observable", "Splitting (meV)", "Status"), (
            (r["sample"], r["observable"], mev(r["value_ev"]), r["status"]) for r in baseline
        ))
        fitted = []
    else:
        targets = reference.load_reference(root / "reference.csv")
        sources = ", ".join(dict.fromkeys(target.source for target in targets.values()))
        lines.extend(("", f"The reference source is {sources}."))
        fitted = read_csv(root / "fit.csv") if study.fit is not None else []
        levels = tuple(dict.fromkeys(r["level"] for r in fitted))
        lookup = {(r["sample"], r["observable"], r["level"]): r for r in fitted}
        for kind in dict.fromkeys(r["kind"] for r in baseline):
            lines.extend(("", f"## {kind.title()} residuals", "",
                          "Values and signed residuals (xTB − reference) are in meV."))
            table(lines, ("Sample", "Observable", "Reference", "σ", "L-1", *levels), (
                (r["sample"], r["observable"], mev(r["reference_ev"]), mev(r["sigma_ev"]), mev(r["residual_ev"]),
                 *(mev(lookup[(r["sample"], r["observable"], level)]["residual_ev"])
                   if (r["sample"], r["observable"], level) in lookup else "—" for level in levels))
                for r in baseline if r["kind"] == kind
            ))

    flags = [r for r in baseline if r["status"] != "ok"]
    if flags:
        lines.extend(("", "## Applicability flags"))
        table(lines, ("Sample", "Observable", "Status", "Detail"), (
            (r["sample"], r["observable"], r["status"], r["error"]) for r in flags
        ))

    lines.extend(("", "## Band-gap diagnostic",
                  "", "ΔE_g is the signed SOC-induced change in the band gap on the DFTB+ SCC mesh, in meV."))
    by_level = {(r["sample"], r["level"]): r["delta_eg_ev"] for r in fitted}
    by_sample = {r["sample"]: r["delta_eg_ev"] for r in baseline}
    diagnostic_levels = tuple(dict.fromkeys(r["level"] for r in fitted))
    table(lines, ("Sample", "L-1", *diagnostic_levels), (
        (name, mev(value), *(mev(by_level.get((name, level), "")) for level in diagnostic_levels))
        for name, value in by_sample.items()
    ))

    if study.fit is not None:
        params = read_csv(root / "parameters.csv")
        lines.extend(("", "## Fit scales and identifiability", "",
                      f"The fit assumes a relative model error of {100 * study.fit.relative_error:g}%. "
                      "Reduced χ² uses the data residuals and subtracts one degree of freedom per fitted scale.",
                      "",
                      "The effect column is the largest estimated splitting change across the full scale range. "
                      "Identifiability requires a change greater than 3σ for at least one eligible observable. "
                      "Standard errors include the Gaussian prior in the Jacobian."))
        if study.fit.per_sample:
            lines.extend(("", "L1 fits use the observables eligible across all samples of the same material."))
        table(lines, ("Fit", "Sample", "Shell", "Scale", "SE", "Identifiable", "At bound", "Max effect (meV)",
                      "Reduced χ²"), (
            (r["level"], r["sample"] or "all", f"{r['element']} {r['shell']}", f"{float(r['scale']):.4f}",
             f"{float(r['stderr']):.4f}" if r["stderr"] else "—",
             r["identifiable"], r["at_bound"], mev(r["max_effect_ev"]),
             f"{float(r['chi2_reduced']):.3f}" if r["chi2_reduced"] else "—")
            for r in params
        ))

    if study.learn is not None:
        loo = read_csv(root / "loo.csv")
        lines.extend(("", "## Leave-one-sample-out evaluation", "",
                      "Held-out mean absolute splitting errors in meV. "
                      "L0 common fits a constant scale using the same observables as L1. "
                      "The L2 folds holding out the shortest and longest distances extrapolate."))
        groups = tuple(dict.fromkeys((r["sample"], r["kind"]) for r in loo))

        def mae(level, sample=None, kind=None):
            return mev(np.mean([abs(float(r["residual_ev"])) for r in loo if r["level"] == level
                                and (sample is None or sample == r["sample"])
                                and (kind is None or kind == r["kind"])]))

        table(lines, ("Held-out sample", "Bands", *LOO_LEVELS), [
            (sample, kind, *(mae(level, sample, kind) for level in LOO_LEVELS)) for sample, kind in groups
        ] + [
            ("All", kind, *(mae(level, kind=kind) for level in LOO_LEVELS))
            for kind in dict.fromkeys(r["kind"] for r in loo)
        ] + [("All", "All", *(mae(level) for level in LOO_LEVELS))])
        materials = {r["sample"]: r["material"] for r in baseline}
        lines.extend(("", "## Splitting trends", "",
                      "Linear slopes in meV/Å against mean nearest-neighbour distance, "
                      "for observables eligible at every sample of the material. "
                      "Slopes for fitted models use held-out predictions."))
        table(lines, ("Material", "Shell", "Observable", "Reference", *LOO_LEVELS),
              splitting_trends(loo, read_csv(root / "descriptors.csv"), materials))
        report_plot(root, params)
        if (root / "scales.png").is_file():
            lines.extend(("", "![L1 scales with standard errors, L0 and all-sample L2](scales.png)"))

    (root / "report.md").write_text("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--dftbplus-executable", type=Path)
    parser.add_argument("--jobs", type=int, default=1, help="Number of samples calculated concurrently")
    parser.add_argument("stage", choices=("relax", "reference", "baseline", "fit", "learn", "report"))
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    config = args.config.resolve()
    root = (args.output_root or Path("examples/generated") / config.stem).resolve()
    study = load_study(config, root)
    stage_settings = {"reference": study.reference, "fit": study.fit, "learn": study.learn}
    if args.stage in stage_settings and stage_settings[args.stage] is None:
        parser.error(f"{args.stage} is not configured for {config.name}")
    if args.stage == "report":
        report_stage(study, root)
    elif args.stage == "reference":
        reference_stage(study, root, config.parent, args.jobs)
    else:
        executable = resolve_dftbplus_path(args.dftbplus_executable or study.executable)
        if args.stage == "relax":
            relax_stage(study, root, config.parent, executable, args.jobs)
        else:
            if args.stage == "baseline" and study.reference is None and not (root / "samples.csv").exists():
                relax_stage(study, root, config.parent, executable, args.jobs)
            with ThreadPoolExecutor(max_workers=args.jobs) as workers:
                model = Forward(study, root, executable, workers)
                if args.stage == "baseline":
                    model.baseline()
                elif args.stage == "fit":
                    fit_stage(model, args.jobs)
                else:
                    learn_stage(model, args.jobs)
    print(root)


if __name__ == "__main__":
    sys.modules["ml_soc_xtb.periodic_workflow"] = sys.modules[__name__]
    main()
