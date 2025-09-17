"""L-point splittings in scalar-relaxed GaAs and GaP."""

from pathlib import Path

from ml_soc_xtb.periodic_dftb import Settings as DftbSettings
from ml_soc_xtb.periodic_fit import Settings as FitSettings
from ml_soc_xtb.periodic_observables import Observable
from ml_soc_xtb.periodic_reference import Resolution, Settings as GpawSettings
from ml_soc_xtb.periodic_workflow import Material, Study


STUDY = Study(
    materials=(
        Material("GaAs", "zincblende", ("GaAs", "zincblende", 5.653)),
        Material("GaP", "zincblende", ("GaP", "zincblende", 5.451)),
    ),
    observables={"zincblende": (
        Observable("l_vb", "L", "valence", 2),
        Observable("l_cb", "L", "conduction", 2),
    )},
    executable=Path(".dftb-24.1/bin/dftb+"),
    dftb=DftbSettings(0.025, 1e-8, 200, 200, 1e-4),
    meshes={"zincblende": ((8, 8, 8), (0.5, 0.5, 0.5))},
    shell_order=dict.fromkeys(("Ga", "As", "P"), ("s", "p", "d")),
    fitted_shell=dict.fromkeys(("Ga", "As", "P"), "p"),
    tolerance_ev=0.0005,
    reference=GpawSettings(base=Resolution(400, (6, 6, 6), 16), tight=Resolution(600, (10, 10, 10), 32)),
    fit=FitSettings(tau=0.25, relative_error=0.1, fd_step=0.02, per_sample=False),
)
