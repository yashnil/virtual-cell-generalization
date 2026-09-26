"""Figure-source files, provenance records and saving.

Every figure is drawn from a small table in ``data/figure_sources/`` that is
extracted deterministically from frozen output artifacts
(:mod:`virtual_cell.visualization.sources`). Each table has a sidecar
``<name>.provenance.json`` recording the artifacts it was read from (with
SHA-256 and the freeze manifests that pin them), the report it illustrates, the
extraction and figure scripts, the date, and the git HEAD.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from collections.abc import Iterable
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
SOURCES_DIR = ROOT / "data" / "figure_sources"
FIGURES_DIR = ROOT / "reports" / "figures"
FREEZE_DIR = ROOT / "data" / "provenance" / "scperteval"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def git_head() -> dict[str, object]:
    def run(*args: str) -> str:
        return subprocess.run(
            ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
        ).stdout.strip()

    return {"commit": run("rev-parse", "HEAD"), "dirty": bool(run("status", "--porcelain"))}


def freezes_pinning(rel_path: str) -> list[str]:
    """Names of the freeze manifests that record a digest for ``rel_path``."""
    names = []
    for manifest in sorted(FREEZE_DIR.glob("*_freeze.txt")):
        for line in manifest.read_text().splitlines():
            parts = line.split()
            if len(parts) == 2 and len(parts[0]) == 64 and parts[1] == rel_path:
                names.append(manifest.name)
                break
    return names


def source_record(paths: Iterable[str]) -> list[dict[str, object]]:
    out = []
    for rel in paths:
        p = ROOT / rel
        out.append({"path": rel, "sha256": sha256(p), "freezes": freezes_pinning(rel)})
    return out


def write_source(
    name: str,
    table: pd.DataFrame,
    *,
    sources: Iterable[str],
    report: str,
    figure_script: str,
    extraction_script: str = "scripts/figures/extract_figure_sources.py",
    notes: str = "",
) -> Path:
    """Write ``data/figure_sources/<name>.csv`` plus its provenance sidecar."""
    SOURCES_DIR.mkdir(parents=True, exist_ok=True)
    csv = SOURCES_DIR / f"{name}.csv"
    table.to_csv(csv, index=False, float_format="%.10g")
    record = {
        "figure_source": f"data/figure_sources/{name}.csv",
        "source_files": source_record(sources),
        "report": report,
        "extraction_script": extraction_script,
        "figure_script": figure_script,
        "generated": date.today().isoformat(),
        "git_head": git_head(),
        "notes": notes,
    }
    (SOURCES_DIR / f"{name}.provenance.json").write_text(json.dumps(record, indent=2) + "\n")
    return csv


def load_source(name: str) -> pd.DataFrame:
    return pd.read_csv(SOURCES_DIR / f"{name}.csv")


def load_provenance(name: str) -> dict[str, object]:
    return json.loads((SOURCES_DIR / f"{name}.provenance.json").read_text())


def numeric_is_finite(table: pd.DataFrame, allow_nan: Iterable[str] = ()) -> bool:
    """True when every numeric value is finite, except NaN in ``allow_nan``
    columns (intentionally missing, e.g. a slope that is undefined)."""
    allowed = set(allow_nan)
    for col in table.select_dtypes(include=[np.number]).columns:
        vals = table[col].to_numpy(dtype=np.float64)
        if np.isinf(vals).any():
            return False
        if col not in allowed and np.isnan(vals).any():
            return False
    return True


def save_figure(fig, name: str) -> list[Path]:
    """Save ``reports/figures/<name>.png`` (300 dpi) and ``.svg``."""
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    out = []
    for ext in ("png", "svg"):
        path = FIGURES_DIR / f"{name}.{ext}"
        fig.savefig(path, metadata={"Date": None} if ext == "svg" else None)
        out.append(path)
    return out
