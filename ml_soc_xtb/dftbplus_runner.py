"""External DFTB+ command execution."""

from __future__ import annotations

import subprocess
from pathlib import Path

import resource

def _dftbplus_run_command(path: Path) -> list[str]:
    """Wrap DFTB+ with a fail-fast child-only stack-limit increase on POSIX."""

    soft, hard = resource.getrlimit(resource.RLIMIT_STACK)
    if soft == hard:
        return [str(path)]
    if hard == resource.RLIM_INFINITY:
        shell_limit = "unlimited"
    else:
        shell_limit = str(int(hard) // 1024)
        if int(shell_limit) <= 0:
            raise RuntimeError(f"Invalid hard stack limit for DFTB+: {hard}")
    return [
        "/bin/sh",
        "-c",
        'ulimit -s "$1" || exit $?; exec "$2"',
        "dftbplus-stack-limit",
        shell_limit,
        str(path),
    ]


def resolve_dftbplus_path(executable: str | Path) -> Path:
    """Resolve the explicitly supplied DFTB+ executable path."""

    path = Path(executable).expanduser().absolute()
    if not path.is_file():
        raise FileNotFoundError(f"DFTB+ executable is not a file: {path}")
    return path


def run_dftbplus(run_dir: str | Path, executable: str | Path) -> subprocess.CompletedProcess[str]:
    """Run DFTB+ in a working directory."""

    path = resolve_dftbplus_path(executable)
    run_dir = Path(run_dir)
    if not run_dir.exists() or not run_dir.is_dir():
        raise FileNotFoundError(f"DFTB+ run directory does not exist: {run_dir}")
    return subprocess.run(
        _dftbplus_run_command(path),
        cwd=run_dir,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
