"""Compatibility gate for the installed ``vcc`` CLI (final-round packaging).

The pipeline uses exactly two CLI surfaces:

* ``vcc prep ... --dry-run --json``: it needs the options in ``REQUIRED_PREP_OPTIONS``,
  so every panel value is passed explicitly;
* ``vcc.prep.run_prep(...)`` under the CLI's own interpreter (``vcc_pack_panel.py``): it
  needs the keywords in ``REQUIRED_RUN_PREP_KWARGS`` and the ``read_h5ad`` hook.

:func:`probe` reads these facts from the installed CLI and :func:`check` compares them
with the requirements. :func:`require` raises :class:`VccCompatibilityError` with
actionable steps if anything is missing. **Nothing adapts silently.** A version other
than the tested one passes only when every required capability is present, and it is
reported as untested.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
from dataclasses import asdict, dataclass, field
from pathlib import Path

TESTED_VERSIONS = ("0.2.0",)
REQUIRED_PREP_OPTIONS = (
    "--genes",
    "--perts",
    "--contexts",
    "--cells-per-pert",
    "--expected-gene-dim",
    "--verify-targets",
    "--check-cell-counts",
    "--dry-run",
    "--json",
)
REQUIRED_RUN_PREP_KWARGS = (
    "input_path",
    "output_path",
    "genes_path",
    "perts_path",
    "cells_per_pert",
    "required_contexts",
    "expected_gene_dim",
    "dry_run",
)

_INTROSPECT = (
    "import inspect, json\n"
    "from importlib.metadata import version\n"
    "import vcc.prep as p\n"
    "print(json.dumps({'package_version': version('vcc-cli'),"
    " 'run_prep_params': list(inspect.signature(p.run_prep).parameters),"
    " 'has_read_h5ad': hasattr(p, 'read_h5ad')}))\n"
)


class VccCompatibilityError(RuntimeError):
    """The installed vcc CLI cannot run the final-round packaging as written."""


@dataclass
class CliInfo:
    executable: str | None = None
    interpreter: str | None = None
    version: str | None = None
    package_version: str | None = None
    prep_help: str = ""
    run_prep_params: list[str] = field(default_factory=list)
    has_read_h5ad: bool = False
    probe_errors: list[str] = field(default_factory=list)


@dataclass
class CompatReport:
    compatible: bool
    version: str | None
    tested_version: bool
    problems: list[str]
    warnings: list[str]
    info: dict

    def as_dict(self) -> dict:
        return asdict(self)


def _version(text: str) -> str | None:
    m = re.search(r"\b(\d+\.\d+\.\d+)\b", text)
    return m.group(1) if m else None


def probe(command: str = "vcc") -> CliInfo:
    """Read version, ``prep`` options and the ``run_prep`` signature from the installed CLI."""
    info = CliInfo()
    exe = shutil.which(command)
    if not exe:
        info.probe_errors.append(f"{command!r} is not on PATH")
        return info
    info.executable = exe
    for args in ([exe, "--version"], [exe, "version"]):
        r = subprocess.run(args, capture_output=True, text=True, check=False)
        if r.returncode == 0 and (v := _version(r.stdout + r.stderr)):
            info.version = v
            break
    r = subprocess.run([exe, "prep", "--help"], capture_output=True, text=True, check=False)
    info.prep_help = r.stdout
    if r.returncode != 0:
        info.probe_errors.append(f"`vcc prep --help` exited {r.returncode}")
    first = Path(exe).read_bytes()[:512].split(b"\n", 1)[0].decode(errors="replace")
    if first.startswith("#!"):
        info.interpreter = first[2:].strip().split()[0]
        r = subprocess.run(
            [info.interpreter, "-c", _INTROSPECT], capture_output=True, text=True, check=False
        )
        if r.returncode == 0:
            meta = json.loads(r.stdout.strip().splitlines()[-1])
            info.package_version = meta["package_version"]
            info.run_prep_params = list(meta["run_prep_params"])
            info.has_read_h5ad = bool(meta["has_read_h5ad"])
        else:
            info.probe_errors.append(
                f"cannot import vcc.prep under {info.interpreter}: {r.stderr.strip()[-300:]}"
            )
    else:
        info.probe_errors.append(f"{exe} is not a Python entry-point script (no shebang)")
    return info


def check(info: CliInfo) -> CompatReport:
    """Compare a probe with the pipeline's requirements. Pure; no side effects."""
    problems = list(info.probe_errors)
    warnings = []
    if info.version is None:
        problems.append("could not read a version from `vcc --version`")
    if info.package_version and info.version and info.package_version != info.version:
        problems.append(
            f"`vcc --version` says {info.version} but the importable vcc-cli is "
            f"{info.package_version} (two installations?)"
        )
    missing = [
        o
        for o in REQUIRED_PREP_OPTIONS
        if not re.search(rf"(^|\s|,){re.escape(o)}\b", info.prep_help)
    ]
    if missing:
        problems.append(f"`vcc prep` lacks required options: {missing}")
    lacking = [k for k in REQUIRED_RUN_PREP_KWARGS if k not in info.run_prep_params]
    if lacking:
        problems.append(f"vcc.prep.run_prep lacks required keywords: {lacking}")
    if not info.has_read_h5ad:
        problems.append("vcc.prep.read_h5ad (the disk-mapped reader hook) is missing")
    tested = info.version in TESTED_VERSIONS
    if not tested and not problems:
        warnings.append(
            f"vcc {info.version} is untested (tested: {', '.join(TESTED_VERSIONS)}); every "
            "required capability is present, and behaviour is unchanged"
        )
    return CompatReport(
        not problems,
        info.version,
        tested,
        problems,
        warnings,
        asdict(info) | {"prep_help": f"<{len(info.prep_help)} chars>"},
    )


def require(info: CliInfo | None = None) -> CompatReport:
    """Probe (unless given) and raise with actionable steps if incompatible."""
    report = check(probe() if info is None else info)
    if not report.compatible:
        raise VccCompatibilityError(
            f"installed vcc CLI (version {report.version or 'unknown'}) is not compatible with "
            "the final-round packaging:\n  - "
            + "\n  - ".join(report.problems)
            + "\nWhat to do: install the tested CLI (`uv tool install vcc-cli==0.2.0`), or read "
            "the new CLI's changelog and update `vcc_pack_panel.py` / `package_final_panel.py` "
            "(with tests) before packaging. The prediction itself does not depend on the CLI."
        )
    return report
