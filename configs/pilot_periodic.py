"""Scientific settings for the GFN1-xTB periodic pilot screen."""

from pathlib import Path

from ml_soc_xtb.periodic_dftb import Settings
from ml_soc_xtb.periodic_observables import Observable
from ml_soc_xtb.periodic_workflow import Material, Study


DFTBPLUS_EXECUTABLE = Path(".dftb-24.1/bin/dftb+")
DFTB_SETTINGS = Settings(
    temperature_ev=0.025,
    scc_tolerance=1e-8,
    max_scc_iterations=200,
    max_geometry_steps=200,
    max_force_ha_bohr=1e-4,
)
# Energy-span tolerance for grouping degenerate levels.
DEGENERACY_TOLERANCE_EV = 0.0005


# Starting geometries in Å. The screen uses the scalar-relaxed cells.
VACUUM_ANGSTROM = 10.0  # on each side
CH_ANGSTROM = 1.09
MATERIALS = (
    Material("GaAs", "zincblende", ("GaAs", "zincblende", 5.653)),
    Material("InSb", "zincblende", ("InSb", "zincblende", 6.479)),
    Material("ZnSe", "zincblende", ("ZnSe", "zincblende", 5.668)),
    Material("Si", "zincblende", ("Si", "diamond", 5.431)),
    Material("Ge", "zincblende", ("Ge", "diamond", 5.658)),
    Material("GaSb", "zincblende", ("GaSb", "zincblende", 6.096)),
    Material("InAs", "zincblende", ("InAs", "zincblende", 6.058)),
    Material("InP", "zincblende", ("InP", "zincblende", 5.869)),
    Material("GaP", "zincblende", ("GaP", "zincblende", 5.451)),
    Material("AlSb", "zincblende", ("AlSb", "zincblende", 6.136)),
    Material("ZnTe", "zincblende", ("ZnTe", "zincblende", 6.104)),
    Material("CdTe", "zincblende", ("CdTe", "zincblende", 6.482)),
    Material("CdSe", "zincblende", ("CdSe", "zincblende", 6.050)),
    Material("alpha-Sn", "zincblende", ("Sn", "diamond", 6.489)),
    Material("MoS2", "tmdc_monolayer", ("MoS2", 3.18, 3.13, VACUUM_ANGSTROM)),
    Material("WS2", "tmdc_monolayer", ("WS2", 3.18, 3.15, VACUUM_ANGSTROM)),
    Material("MoSe2", "tmdc_monolayer", ("MoSe2", 3.32, 3.34, VACUUM_ANGSTROM)),
    Material("WSe2", "tmdc_monolayer", ("WSe2", 3.32, 3.35, VACUUM_ANGSTROM)),
    Material("MoTe2", "tmdc_monolayer", ("MoTe2", 3.52, 3.60, VACUUM_ANGSTROM)),
    Material("WTe2", "tmdc_monolayer", ("WTe2", 3.52, 3.60, VACUUM_ANGSTROM)),
    Material("Me-Bi", "ti_bilayer", ("Bi", 5.33, 0.9, 2.25, CH_ANGSTROM, VACUUM_ANGSTROM)),
    Material("Me-Sb", "ti_bilayer", ("Sb", 5.09, 0.9, 2.15, CH_ANGSTROM, VACUUM_ANGSTROM)),
    # Pb geometry seed: Table 2, DOI 10.1039/C9RA07531C.
    Material("Me-Pb", "ti_bilayer", ("Pb", 5.07407, 0.8863, 2.31, CH_ANGSTROM, VACUUM_ANGSTROM)),
)

OBSERVABLES = {
    "zincblende": (
        Observable("gamma_vb", "G", "valence", 3),
        Observable("l_vb", "L", "valence", 2),
        Observable("gamma_cb", "G", "conduction", 3),
        Observable("l_cb", "L", "conduction", 2),
    ),
    "tmdc_monolayer": (
        Observable("k_vb", "K", "valence", 1),
        Observable("k_cb", "K", "conduction", 1),
    ),
    "ti_bilayer": (Observable("k_dirac", "K", "dirac", 2),),
}

# SupercellFolding dimensions and shifts in units of the mesh step.
K_MESHES = {
    "zincblende": ((8, 8, 8), (0.5, 0.5, 0.5)),
    "tmdc_monolayer": ((12, 12, 1), (0, 0, 0)),
    "ti_bilayer": ((12, 12, 1), (0, 0, 0)),
}

# Shells from the installed tblite GFN1 export. Mo and W slot order also
# verified with isolated SOC perturbations in DFTB+ on MoS2 and WS2.
SHELL_ORDER = {
    "H": ("s", "s2"),
    "C": ("s", "p"),
    "Al": ("s", "p", "d"),
    "Si": ("s", "p", "d"),
    "P": ("s", "p", "d"),
    "S": ("s", "p", "d"),
    "Zn": ("s", "p"),
    "Ga": ("s", "p", "d"),
    "Ge": ("s", "p", "d"),
    "As": ("s", "p", "d"),
    "Se": ("s", "p", "d"),
    "Mo": ("d", "s", "p"),
    "Cd": ("s", "p"),
    "In": ("s", "p", "d"),
    "Sn": ("s", "p", "d"),
    "Sb": ("s", "p", "d"),
    "Te": ("s", "p", "d"),
    "W": ("d", "s", "p"),
    "Pb": ("s", "p"),
    "Bi": ("s", "p"),
}

# Apply the tabulated constants to these basis shells at scale 1.
# Group 1–12 tabulated p values are treated as semicore and excluded.
# Zn and Cd have no GFN1 d shell.
FITTED_SHELL = {
    "H": None, "C": "p", "Al": "p", "Si": "p", "P": "p", "S": "p",
    "Zn": None, "Ga": "p", "Ge": "p", "As": "p", "Se": "p", "Mo": "d",
    "Cd": None, "In": "p", "Sn": "p", "Sb": "p", "Te": "p", "W": "d",
    "Pb": "p", "Bi": "p",
}

STUDY = Study(
    materials=MATERIALS, observables=OBSERVABLES, executable=DFTBPLUS_EXECUTABLE,
    dftb=DFTB_SETTINGS, meshes=K_MESHES, shell_order=SHELL_ORDER,
    fitted_shell=FITTED_SHELL, tolerance_ev=DEGENERACY_TOLERANCE_EV,
    preserve_cubic_basis=False,
)
