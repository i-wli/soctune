from pathlib import Path
from unittest.mock import Mock

import pytest

import ml_soc_xtb.dftbplus_runner as dftbplus_runner
from ml_soc_xtb.dftbplus_runner import _dftbplus_run_command


def test_run_resolves_relative_executable_before_changing_directory(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    executable = tmp_path / "dftb+"
    executable.touch()
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    run = Mock()
    monkeypatch.setattr(dftbplus_runner.subprocess, "run", run)

    dftbplus_runner.run_dftbplus(run_dir, executable="dftb+")

    assert run.call_args.args[0][-1] == str(executable)
    assert run.call_args.kwargs["cwd"] == run_dir


def test_invalid_executable_does_not_use_environment_or_path(tmp_path, monkeypatch):
    executable = tmp_path / "dftb+"
    executable.touch(mode=0o755)
    monkeypatch.setenv("PATH", str(tmp_path))
    monkeypatch.setenv("DFTBPLUS_EXECUTABLE", str(executable))
    monkeypatch.setenv("DFTBPLUS_BIN", str(executable))
    missing = tmp_path / "missing"

    with pytest.raises(FileNotFoundError, match="DFTB\\+ executable is not a file"):
        dftbplus_runner.run_dftbplus(tmp_path, executable=missing)


def test_posix_stack_limit_wrapper_uses_hard_limit_without_preexec(monkeypatch):
    monkeypatch.setattr(
        dftbplus_runner.resource,
        "getrlimit",
        lambda limit: (1024, 8192),
    )

    command = _dftbplus_run_command(Path("/tmp/dftb+"))

    assert command == [
        "/bin/sh",
        "-c",
        'ulimit -s "$1" || exit $?; exec "$2"',
        "dftbplus-stack-limit",
        "8",
        "/tmp/dftb+",
    ]
