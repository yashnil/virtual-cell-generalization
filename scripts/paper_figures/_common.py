# ruff: noqa: E501 (figure labels and provenance notes are kept as single strings)
"""Provenance, integrity checks and saving for the paper-figure pipeline.

* Every frozen input is hashed. If the input is listed in its phase's output manifest
  (``outputs/<phase>/manifest_sha256.txt``) or a freeze file under ``data/provenance``, the hash must match,
  otherwise the build stops. Nothing upstream is ever written.
* Every figure-source table is written to ``data/figure_sources/paper/<name>.csv`` with a
  ``<name>.provenance.json`` sidecar.
* Every figure is written to ``reports/paper_figures/<name>.{svg,png,pdf}`` with a
  ``<name>.manifest.json`` sidecar listing the source tables (with SHA-256) it was drawn from.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from collections.abc import Iterable
from datetime import date
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = ROOT / "data" / "figure_sources" / "paper"
OUT_DIR = ROOT / "reports" / "paper_figures"


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


def _manifest_entries() -> dict[str, tuple[str, str]]:
    """rel path (from repo root) -> (sha256, manifest that pins it)."""
    out: dict[str, tuple[str, str]] = {}
    for man in sorted((ROOT / "outputs").glob("*/manifest_sha256.txt")):
        base = man.parent.relative_to(ROOT)
        for line in man.read_text().splitlines():
            parts = line.split()
            if len(parts) == 2 and len(parts[0]) == 64:
                rel = (base / parts[1].removeprefix("./")).as_posix()
                out[rel] = (parts[0], man.relative_to(ROOT).as_posix())
    for man in sorted((ROOT / "data" / "provenance").rglob("*.txt")):
        for line in man.read_text(errors="ignore").splitlines():
            parts = line.split()
            if len(parts) == 2 and len(parts[0]) == 64 and parts[1] not in out:
                out[parts[1].removeprefix("./")] = (parts[0], man.relative_to(ROOT).as_posix())
    return out


_MANIFEST: dict[str, tuple[str, str]] | None = None


def frozen(rel: str) -> Path:
    """Return the path of a frozen input after checking it against any manifest that pins it."""
    global _MANIFEST
    if _MANIFEST is None:
        _MANIFEST = _manifest_entries()
    p = ROOT / rel
    if not p.exists():
        raise FileNotFoundError(rel)
    if rel in _MANIFEST:
        want, man = _MANIFEST[rel]
        got = sha256(p)
        if got != want:
            raise RuntimeError(f"{rel} does not match {man}: {got} != {want}")
    return p


def source_record(paths: Iterable[str]) -> list[dict[str, object]]:
    global _MANIFEST
    if _MANIFEST is None:
        _MANIFEST = _manifest_entries()
    rec = []
    for rel in paths:
        entry = {"path": rel, "sha256": sha256(ROOT / rel)}
        if rel in _MANIFEST:
            entry["verified_against"] = _MANIFEST[rel][1]
        rec.append(entry)
    return rec


def write_source(
    name: str,
    table: pd.DataFrame,
    *,
    sources: Iterable[str],
    reports: Iterable[str],
    build_script: str,
    plot_script: str,
    notes: str,
    report_values: dict | None = None,
) -> Path:
    """Write ``data/figure_sources/paper/<name>.csv`` plus its provenance sidecar.

    ``report_values`` records numbers transcribed from a frozen report (never computed here), with the report
    section they come from; used only where no machine-readable frozen artefact holds them.
    """
    SRC_DIR.mkdir(parents=True, exist_ok=True)
    csv = SRC_DIR / f"{name}.csv"
    table.to_csv(csv, index=False, float_format="%.10g")
    rec = {
        "figure_source": csv.relative_to(ROOT).as_posix(),
        "sha256": sha256(csv),
        "frozen_inputs": source_record(sources),
        "reports": list(reports),
        "build_script": build_script,
        "plot_script": plot_script,
        "generated": date.today().isoformat(),
        "git_head": git_head(),
        "notes": notes,
    }
    if report_values:
        rec["transcribed_from_reports"] = report_values
    (SRC_DIR / f"{name}.provenance.json").write_text(
        json.dumps(rec, indent=2, ensure_ascii=False) + "\n"
    )
    return csv


def load(name: str) -> pd.DataFrame:
    return pd.read_csv(SRC_DIR / f"{name}.csv")


def save(
    fig, name: str, *, sources: Iterable[str], plot_script: str, upstream: Iterable[str]
) -> list[Path]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    paths = []
    for ext in ("svg", "png", "pdf"):
        p = OUT_DIR / f"{name}.{ext}"
        meta = (
            {"Date": None} if ext == "svg" else ({"CreationDate": None} if ext == "pdf" else None)
        )
        fig.savefig(p, metadata=meta, dpi=600 if ext == "png" else None)
        paths.append(p)
    srcs = [f"data/figure_sources/paper/{s}.csv" for s in sources]
    man = {
        "figure": name,
        "outputs": {
            p.suffix[1:]: {"path": p.relative_to(ROOT).as_posix(), "sha256": sha256(p)}
            for p in paths
        },
        "figure_sources": source_record(srcs),
        "plot_script": plot_script,
        "frozen_upstream": list(upstream),
        "generated": date.today().isoformat(),
        "git_head": git_head(),
    }
    (OUT_DIR / f"{name}.manifest.json").write_text(
        json.dumps(man, indent=2, ensure_ascii=False) + "\n"
    )
    return paths
