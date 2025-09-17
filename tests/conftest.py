from pathlib import Path

import pytest


def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption(
        "--dftbplus-executable",
        type=Path,
        help="Path to the DFTB+ executable for opt-in integration tests.",
    )


@pytest.fixture(scope="module")
def dftbplus_executable(pytestconfig: pytest.Config) -> Path:
    executable = pytestconfig.getoption("--dftbplus-executable")
    if executable is None:
        pytest.skip("Pass --dftbplus-executable to run DFTB+ integration tests.")
    return executable
