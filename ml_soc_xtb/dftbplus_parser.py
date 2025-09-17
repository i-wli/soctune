"""Parsers for small DFTB+ output values."""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from .parser_patterns import FLOAT_RE

BAND_HEADER_RE = re.compile(r"^\s*KPT\s+(\d+)\s+SPIN\s+(\d+)\b")
BAND_ROW_RE = re.compile(
    rf"^\s*(\d+)\s+({FLOAT_RE})\s+({FLOAT_RE})\s*$",
    re.MULTILINE,
)


def _parse_float(text: str) -> float:
    return float(text.replace("D", "E").replace("d", "e"))


def parse_band_out(text: str) -> pd.DataFrame:
    """Parse DFTB+ band.out text into k-point, spin, state, energy, and occupation rows."""

    rows: list[dict[str, float | int]] = []
    current_kpt: int | None = None
    current_spin: int | None = None

    for line in text.splitlines():
        header_match = BAND_HEADER_RE.match(line)
        if header_match:
            current_kpt = int(header_match.group(1))
            current_spin = int(header_match.group(2))
            continue

        row_match = BAND_ROW_RE.match(line)
        if row_match and current_kpt is not None and current_spin is not None:
            rows.append(
                {
                    "kpt": current_kpt,
                    "spin": current_spin,
                    "state": int(row_match.group(1)),
                    "energy_ev": _parse_float(row_match.group(2)),
                    "occupation": _parse_float(row_match.group(3)),
                }
            )

    if not rows:
        raise ValueError("No valid DFTB+ band rows found in band.out")
    return pd.DataFrame(rows, columns=["kpt", "spin", "state", "energy_ev", "occupation"])


def parse_band_out_file(path: str | Path) -> pd.DataFrame:
    """Parse a DFTB+ band.out file."""

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"DFTB+ band output file does not exist: {path}")
    return parse_band_out(path.read_text(encoding="utf-8"))
