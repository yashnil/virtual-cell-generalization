"""Validate and record the C0_ATLASSHIFT_REPRODUCTION bundle (section 11). NOT submitted.

Applies the thirteen local checks used for v1 (``scripts/audit_arc_submission_candidate.py``,
via the frozen ``inspect_bundle`` path) to the uncompacted prediction, confirms the
target set equals ``pert_counts.csv``, runs ``vcc prep --dry-run --json`` on the
compacted file, and writes every checksum needed to reproduce the bundle:

* prediction / compact / ``.vcc`` SHA-256;
* upstream code commit and every vendored file's SHA-256, plus this repo's git HEAD;
* every raw source checksum (from the verified download log) and every prepared
  statistics file SHA-256.

Output: ``outputs/competition_v2/atlasshift_c0/c0_manifest.json``. It sets
``"submitted": false`` and nothing in this repository submits it.

Reproduce: ``uv run python scripts/competition_v2/validate_c0_bundle.py``
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import h5py
import pandas as pd

from virtual_cell.arc import bundle, generate
from virtual_cell.competition_v2 import ATLASSHIFT_COMMIT, XATLAS_REVISION

ROOT = Path(__file__).resolve().parents[2]
C0 = ROOT / "outputs" / "competition_v2" / "atlasshift_c0"
CONTROLS = ROOT / "data" / "raw" / "arc2026" / "controls"
UPSTREAM = ROOT / "third_party" / "atlasshift"
LOG = ROOT / "data" / "provenance" / "competition_v2" / "download_log.jsonl"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(16 * 1024**2), b""):
            h.update(block)
    return h.hexdigest()


def main() -> int:
    genes = pd.read_csv(CONTROLS / "gene_names.csv").gene_name.astype(str).tolist()
    panel = set(pd.read_csv(CONTROLS / "pert_counts.csv").target_gene.astype(str))
    prediction = C0 / "prediction.h5ad"
    compact = C0 / "prediction_compact.h5ad"
    package = C0 / "c0_atlasshift_reproduction_val.vcc"

    report = bundle.inspect_bundle(prediction, genes)
    sizes = report.cells_per_group
    checks = {
        "n_cells == 360000": report.n_cells == 360_000,
        "n_genes == 18533": report.n_genes == 18_533,
        "gene order matches gene_names.csv": report.gene_order_matches,
        "raw integer counts": report.integer_valued,
        "non-negative": report.non_negative,
        "finite": report.finite,
        "exactly 400 cells per (context, perturbation)": bool((sizes == 400).all()),
        "300 perturbations in every context": all(
            v == 300 for v in report.n_perturbations_per_context.values()
        ),
        "contexts are exactly A, B, C": report.contexts == ("A", "B", "C"),
        "no control cells emitted": "non-targeting" not in set(sizes.index.get_level_values(1)),
        "under max_counts_per_cell (1,000,000)": report.max_counts_per_cell
        < generate.MAX_COUNTS_PER_CELL,
        "under max_nnz (4,750,000,000)": report.nnz < 4_750_000_000,
        "under max_cell_dim (400,000)": report.n_cells <= 400_000,
    }
    targets = set(sizes.index.get_level_values(1))
    checks["target set equals pert_counts.csv"] = targets == panel
    with h5py.File(prediction, "r") as f:
        zeros = 0
        data = f["X/data"]
        for left in range(0, data.shape[0], 1 << 26):
            zeros += int((data[left : left + (1 << 26)] == 0).sum())
    checks["no explicit stored zeros"] = zeros == 0

    dry = subprocess.run(
        [
            "vcc",
            "prep",
            str(compact),
            "-g",
            str(CONTROLS / "gene_names.csv"),
            "--perts",
            str(CONTROLS / "pert_counts.csv"),
            "--dry-run",
            "--json",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    try:
        dry_report = json.loads(dry.stdout)
    except json.JSONDecodeError:
        dry_report = {"stdout": dry.stdout[-2000:], "stderr": dry.stderr[-2000:]}
    checks["vcc prep --dry-run exit 0"] = dry.returncode == 0
    (C0 / "vcc_prep_dry_run.json").write_text(json.dumps(dry_report, indent=2) + "\n")

    downloads = {}
    for line in LOG.read_text().splitlines():
        rec = json.loads(line)
        if str(rec.get("status", "")).startswith(("downloaded", "present")):
            downloads[rec["source"]] = {
                k: rec[k] for k in ("sha256", "bytes", "revision") if k in rec
            }
    head = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True
    ).stdout.strip()
    manifest = {
        "candidate": "C0_ATLASSHIFT_REPRODUCTION",
        "submitted": False,
        "provenance": f"external public baseline AtlasShift @ {ATLASSHIFT_COMMIT} (MIT), "
        "run unmodified",
        "partition_panel": "val / vcc2026-val-1",
        "prediction": {
            "path": str(prediction.relative_to(ROOT)),
            "sha256": sha256(prediction),
            "bytes": prediction.stat().st_size,
        },
        "compact": {
            "path": str(compact.relative_to(ROOT)),
            "sha256": sha256(compact),
            "bytes": compact.stat().st_size,
        },
        "vcc_package": {
            "path": str(package.relative_to(ROOT)),
            "sha256": sha256(package),
            "bytes": package.stat().st_size,
        }
        if package.exists()
        else None,
        "code": {
            "upstream_commit": ATLASSHIFT_COMMIT,
            "upstream_files_sha256": {
                p.name: sha256(p) for p in sorted(UPSTREAM.glob("*")) if p.is_file()
            },
            "repo_git_head": head,
            "upstream_python": "3.13 (third_party/atlasshift/.venv, requirements.txt pins)",
            "upstream_constants": {
                "weights": "K562 2, HCT116 1, HEK293T 1, H1 2; CD4 0.5",
                "amplitude": "log2fc 0.6, bulk_delta 0.3",
                "clip": 3,
                "prior_counts": 100000,
                "minimum_cells": 20,
                "promoter_fraction": 0.15,
                "seed": 20260910,
            },
        },
        "source_data": {
            "raw_downloads": downloads,
            "xatlas_revision": XATLAS_REVISION,
            "prepared_statistics_sha256": {
                p.name: sha256(p) for p in sorted((C0 / "data").glob("*.npz"))
            },
            "official_pairs_sha256": sha256(C0 / "data" / "official_pairs.csv"),
        },
        "bundle": {**report.as_dict(), "explicit_zeros": zeros},
        "checks": checks,
        "all_checks_pass": all(checks.values()),
        "vcc_prep_dry_run": dry_report,
    }
    (C0 / "c0_manifest.json").write_text(json.dumps(manifest, indent=2, default=float) + "\n")
    for k, v in checks.items():
        print(f"  [{'PASS' if v else 'FAIL'}] {k}")
    print(json.dumps({k: manifest[k] for k in ("prediction", "compact", "vcc_package")}, indent=2))
    return 0 if manifest["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
