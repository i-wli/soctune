"""SOC reference constants loading and validation."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

VALID_SHELLS = ("s", "p", "d", "f")
REQUIRED_COLUMNS = ("element", "shell", "xi_hartree")
VALID_ELEMENT_SYMBOLS = frozenset(
    (
        "H",
        "He",
        "Li",
        "Be",
        "B",
        "C",
        "N",
        "O",
        "F",
        "Ne",
        "Na",
        "Mg",
        "Al",
        "Si",
        "P",
        "S",
        "Cl",
        "Ar",
        "K",
        "Ca",
        "Sc",
        "Ti",
        "V",
        "Cr",
        "Mn",
        "Fe",
        "Co",
        "Ni",
        "Cu",
        "Zn",
        "Ga",
        "Ge",
        "As",
        "Se",
        "Br",
        "Kr",
        "Rb",
        "Sr",
        "Y",
        "Zr",
        "Nb",
        "Mo",
        "Tc",
        "Ru",
        "Rh",
        "Pd",
        "Ag",
        "Cd",
        "In",
        "Sn",
        "Sb",
        "Te",
        "I",
        "Xe",
        "Cs",
        "Ba",
        "La",
        "Ce",
        "Pr",
        "Nd",
        "Pm",
        "Sm",
        "Eu",
        "Gd",
        "Tb",
        "Dy",
        "Ho",
        "Er",
        "Tm",
        "Yb",
        "Lu",
        "Hf",
        "Ta",
        "W",
        "Re",
        "Os",
        "Ir",
        "Pt",
        "Au",
        "Hg",
        "Tl",
        "Pb",
        "Bi",
        "Po",
        "At",
        "Rn",
        "Fr",
        "Ra",
        "Ac",
        "Th",
        "Pa",
        "U",
        "Np",
        "Pu",
        "Am",
        "Cm",
        "Bk",
        "Cf",
        "Es",
        "Fm",
        "Md",
        "No",
        "Lr",
        "Rf",
        "Db",
        "Sg",
        "Bh",
        "Hs",
        "Mt",
        "Ds",
        "Rg",
        "Cn",
        "Nh",
        "Fl",
        "Mc",
        "Lv",
        "Ts",
        "Og",
    )
)


def default_soc_constants_path() -> Path:
    """Return the repository default constants table path."""

    return Path(__file__).resolve().parents[1] / "data" / "soc_parameters.csv"


def load_soc_constants(path: str | Path) -> pd.DataFrame:
    """Load and validate a Hartree-valued shell-resolved SOC constants table."""

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"SOC constants file does not exist: {path}")

    df = pd.read_csv(path)
    missing = [column for column in REQUIRED_COLUMNS if column not in df.columns]
    if missing:
        raise ValueError(f"SOC constants table is missing required columns: {missing}")
    if df.empty:
        raise ValueError("SOC constants table is empty")

    df = df.loc[:, REQUIRED_COLUMNS].copy()

    for idx, element in df["element"].items():
        if not isinstance(element, str) or element not in VALID_ELEMENT_SYMBOLS:
            raise ValueError(f"Invalid element symbol at row {idx}: {element!r}")

    for idx, shell in df["shell"].items():
        if shell not in VALID_SHELLS:
            raise ValueError(f"Invalid shell at row {idx}: {shell!r}; expected one of {VALID_SHELLS}")

    xi = pd.to_numeric(df["xi_hartree"], errors="coerce")
    if not np.isfinite(xi.to_numpy(dtype=float)).all():
        raise ValueError("SOC constants table contains non-finite xi_hartree values")
    if (xi < 0.0).any():
        raise ValueError("SOC constants table contains negative xi_hartree values")
    df["xi_hartree"] = xi.astype(float)

    duplicates = df.duplicated(["element", "shell"], keep=False)
    if duplicates.any():
        pairs = df.loc[duplicates, ["element", "shell"]].drop_duplicates().to_dict("records")
        raise ValueError(f"SOC constants table contains duplicate element-shell pairs: {pairs}")

    return df


def get_xi_ref(element: str, shell: str, df: pd.DataFrame) -> float:
    """Return a reference SOC constant in Hartree for one element-shell pair."""

    if element not in VALID_ELEMENT_SYMBOLS:
        raise ValueError(f"Invalid element symbol: {element!r}")
    if shell not in VALID_SHELLS:
        raise ValueError(f"Invalid shell: {shell!r}; expected one of {VALID_SHELLS}")

    matches = df[(df["element"] == element) & (df["shell"] == shell)]
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one SOC constant for {element}-{shell}, found {len(matches)}")
    value = float(matches.iloc[0]["xi_hartree"])
    if not np.isfinite(value):
        raise ValueError(f"Non-finite SOC constant for {element}-{shell}")
    return value
