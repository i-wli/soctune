"""Diamond Si at five hydrostatic strains."""

from pathlib import Path

from ml_soc_xtb.periodic_dftb import Settings as DftbSettings
from ml_soc_xtb.periodic_fit import Settings as FitSettings
from ml_soc_xtb.periodic_learn import Settings as LearnSettings
from ml_soc_xtb.periodic_observables import Observable
from ml_soc_xtb.periodic_reference import Resolution, Settings as GpawSettings
from ml_soc_xtb.periodic_workflow import Material, Study


STUDY = Study(
    materials=(Material("Si", "zincblende", ("Si", "diamond", 5.431),
                        strains=(-0.03, -0.015, 0.0, 0.015, 0.03)),),
    observables={"zincblende": (
        Observable("gamma_vb", "G", "valence", 3),
        Observable("l_vb", "L", "valence", 2),
        Observable("gamma_cb", "G", "conduction", 3),
        Observable("l_cb", "L", "conduction", 2),
    )},
    executable=Path(".dftb-24.1/bin/dftb+"),
    dftb=DftbSettings(0.025, 1e-8, 200, 200, 1e-4),
    meshes={"zincblende": ((8, 8, 8), (0.5, 0.5, 0.5))},
    shell_order={"Si": ("s", "p", "d")},
    fitted_shell={"Si": "p"},
    tolerance_ev=0.0005,
    reference=GpawSettings(base=Resolution(400, (6, 6, 6), 16), tight=Resolution(600, (10, 10, 10), 32)),
    fit=FitSettings(tau=0.25, relative_error=0.1, fd_step=0.02),
    learn=LearnSettings(alpha=1.0),
)
