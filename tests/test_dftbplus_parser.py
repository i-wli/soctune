from pathlib import Path

import pytest

from ml_soc_xtb.dftbplus_parser import parse_band_out, parse_band_out_file


def test_parse_band_out_i2_format():
    text = """
 KPT            1  SPIN            1  KWEIGHT    1.0000000000000000
     1    -24.1822  1.00000
     2    -24.1822  1.00000
    14    -10.4246  1.00000
    15     -7.9611  0.00000
"""
    bands = parse_band_out(text)
    assert list(bands.columns) == ["kpt", "spin", "state", "energy_ev", "occupation"]
    assert bands.to_dict("records") == [
        {"kpt": 1, "spin": 1, "state": 1, "energy_ev": -24.1822, "occupation": 1.0},
        {"kpt": 1, "spin": 1, "state": 2, "energy_ev": -24.1822, "occupation": 1.0},
        {"kpt": 1, "spin": 1, "state": 14, "energy_ev": -10.4246, "occupation": 1.0},
        {"kpt": 1, "spin": 1, "state": 15, "energy_ev": -7.9611, "occupation": 0.0},
    ]


def test_parse_band_out_file(tmp_path: Path):
    path = tmp_path / "band.out"
    path.write_text(
        " KPT 1 SPIN 1 KWEIGHT 1.0\n"
        "     1    -1.0  1.00000\n"
        "     2     2.0  0.00000\n",
        encoding="utf-8",
    )
    bands = parse_band_out_file(path)
    assert bands["energy_ev"].tolist() == [-1.0, 2.0]


def test_parse_band_out_rejects_empty_or_malformed_text():
    with pytest.raises(ValueError, match="No valid DFTB\\+ band rows"):
        parse_band_out("")
    with pytest.raises(ValueError, match="No valid DFTB\\+ band rows"):
        parse_band_out("KPT 1 SPIN 1\nnot a band row\n")
