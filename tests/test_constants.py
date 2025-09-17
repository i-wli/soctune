from pathlib import Path

import pytest

from ml_soc_xtb.constants import default_soc_constants_path, get_xi_ref, load_soc_constants


def test_load_default_soc_constants():
    df = load_soc_constants(default_soc_constants_path())
    assert get_xi_ref("I", "p", df) == pytest.approx(0.0328189165248978)
    assert get_xi_ref("W", "d", df) == pytest.approx(0.0113279406280928)


def test_load_soc_constants_rejects_duplicate_pairs(tmp_path: Path):
    path = tmp_path / "bad.csv"
    path.write_text(
        "element,shell,xi_hartree\nC,p,0.001\nC,p,0.002\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="duplicate"):
        load_soc_constants(path)


def test_load_soc_constants_rejects_invalid_rows(tmp_path: Path):
    path = tmp_path / "bad.csv"
    path.write_text(
        "element,shell,xi_hartree\nXx,p,0.001\nC,g,0.002\nI,p,nan\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError):
        load_soc_constants(path)


def test_get_xi_ref_rejects_missing_and_bad_shell():
    df = load_soc_constants(default_soc_constants_path())
    with pytest.raises(ValueError, match="exactly one"):
        get_xi_ref("Pb", "d", df)
    with pytest.raises(ValueError, match="Invalid shell"):
        get_xi_ref("C", "g", df)


def test_load_soc_constants_rejects_negative_value(tmp_path: Path):
    path = tmp_path / "bad.csv"
    path.write_text("element,shell,xi_hartree\nC,p,-0.001\n", encoding="utf-8")
    with pytest.raises(ValueError, match="negative"):
        load_soc_constants(path)
