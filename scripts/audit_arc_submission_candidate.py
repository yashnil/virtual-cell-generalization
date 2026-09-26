"""Pre-submission audit of the frozen Arc dry-run prediction (read-only).

Locks ``outputs/arc_dry_run_v1/arc_dry_run_v1.h5ad`` as submission candidate v1
without regenerating it:

1. verifies the official control bundle (checksums, 300 targets, 18,533 genes,
   contexts A/B/C, 400 cells per perturbation);
2. records file identity (size, mtime, SHA-256) and full streaming statistics
   of the prediction (counts, library sizes, genes detected, label checksums);
3. re-runs the thirteen frozen local checks through the frozen
   ``virtual_cell.arc.bundle.inspect_bundle`` code path;
4. re-runs ``vcc prep --dry-run`` and requires its report to equal the frozen
   ``outputs/arc_dry_run_v1/vcc_prep_dry_run.json`` exactly.

Nothing is fitted, generated, packaged or submitted. Writes
``outputs/arc_submission_v1/candidate_audit.json`` and exits non-zero on any
failure.

Reproduce: ``uv run python scripts/audit_arc_submission_candidate.py``
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from virtual_cell.arc import bundle, generate
from virtual_cell.data import arc2026

ROOT = Path(__file__).resolve().parents[1]
CONTROLS = ROOT / "data" / "raw" / "arc2026" / "controls"
CANDIDATE = ROOT / "outputs" / "arc_dry_run_v1" / "arc_dry_run_v1.h5ad"
FROZEN_DRY_RUN = ROOT / "outputs" / "arc_dry_run_v1" / "vcc_prep_dry_run.json"
FROZEN_SUMMARY = ROOT / "outputs" / "arc_dry_run_v1" / "summary.json"
CONTROL_SHA = ROOT / "data" / "provenance" / "arc2026_controls_sha256.txt"
OUTDIR = ROOT / "outputs" / "arc_submission_v1"
CHUNK = 5_000


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 24), b""):
            h.update(block)
    return h.hexdigest()


def sha256_lines(values) -> str:
    return hashlib.sha256("\n".join(map(str, values)).encode()).hexdigest()


def rule(title: str) -> None:
    print("\n" + "=" * 78 + f"\n{title}\n" + "=" * 78, flush=True)


def read_obs_column(f: h5py.File, name: str) -> np.ndarray:
    n = f["obs"][name]
    if isinstance(n, h5py.Group) and "categories" in n:
        cats = np.array(
            [c.decode() if isinstance(c, bytes) else str(c) for c in n["categories"][:]],
            dtype=object,
        )
        return cats[n["codes"][:]]
    raw = n["values"] if isinstance(n, h5py.Group) else n
    return np.array([v.decode() if isinstance(v, bytes) else str(v) for v in raw[:]], dtype=object)


def describe(x: np.ndarray) -> dict[str, float]:
    return {
        "min": float(x.min()),
        "median": float(np.median(x)),
        "mean": float(x.mean()),
        "max": float(x.max()),
    }


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    failures: list[str] = []

    rule("D. OFFICIAL INPUT BUNDLE")
    control_hashes = {}
    for line in CONTROL_SHA.read_text().splitlines():
        parts = line.split()
        if len(parts) == 2 and len(parts[0]) == 64:
            actual = sha256_file(ROOT / parts[1])
            control_hashes[parts[1]] = {
                "recorded": parts[0],
                "actual": actual,
                "ok": actual == parts[0],
            }
            print(f"  [{'OK' if actual == parts[0] else 'MISMATCH'}] {parts[1]}")
            if actual != parts[0]:
                failures.append(f"control checksum mismatch: {parts[1]}")
    manifest = json.loads((CONTROLS / "manifest.json").read_text())
    genes = list(arc2026.load_gene_names(CONTROLS))
    targets = list(arc2026.load_pert_counts(CONTROLS)["target_gene"])
    inputs = {
        "partition": manifest["partition"],
        "panel_id": manifest["panel_id"],
        "contexts": manifest["contexts"],
        "n_genes": len(genes),
        "n_targets": len(targets),
        "cells_per_pert": manifest["cells_per_pert"],
        "gene_order_sha256": sha256_lines(genes),
        "target_list_sha256": sha256_lines(targets),
    }
    expect = {
        "contexts": ["A", "B", "C"],
        "n_genes": 18_533,
        "n_targets": 300,
        "cells_per_pert": 400,
    }
    for k, v in expect.items():
        ok = inputs[k] == v
        print(f"  [{'OK' if ok else 'FAIL'}] {k} = {inputs[k]}")
        if not ok:
            failures.append(f"input {k} = {inputs[k]} != {v}")
    if len(set(targets)) != 300 or "non-targeting" in targets:
        failures.append("target list is not 300 unique non-control targets")

    rule("E. CANDIDATE IDENTITY AND STATISTICS")
    st = CANDIDATE.stat()
    identity = {
        "path": str(CANDIDATE),
        "bytes": st.st_size,
        "mtime_utc": datetime.fromtimestamp(st.st_mtime, UTC).isoformat(),
        "sha256": sha256_file(CANDIDATE),
    }
    print(json.dumps(identity, indent=2))
    frozen_summary = json.loads(FROZEN_SUMMARY.read_text())
    if st.st_size != frozen_summary["bundle_bytes"]:
        failures.append("candidate size differs from the frozen summary")

    with h5py.File(CANDIDATE, "r") as f:
        shape = tuple(int(v) for v in f["X"].attrs["shape"])
        encoding = f["X"].attrs.get("encoding-type")
        dtype = str(f["X/data"].dtype)
        index_dtype = str(f["X/indices"].dtype)
        indptr_dtype = str(f["X/indptr"].dtype)
        nnz = int(f["X/data"].shape[0])
        var = f["var"]
        node = var[var.attrs.get("_index", "_index")]
        values = node["values"] if isinstance(node, h5py.Group) else node
        cand_genes = [v.decode() if isinstance(v, bytes) else str(v) for v in values[:]]
        perts = read_obs_column(f, "target_gene")
        contexts = read_obs_column(f, "context")

    n_min = np.inf
    n_max = -np.inf
    n_nan = n_inf = n_neg = n_frac = n_zero = 0
    lib, det = [], []
    for _start, chunk in arc2026.stream_row_chunks(CANDIDATE, chunk_size=CHUNK):
        d = chunk.data
        if d.size:
            n_nan += int(np.isnan(d).sum())
            n_inf += int(np.isinf(d).sum())
            fin = d[np.isfinite(d)]
            n_neg += int((fin < 0).sum())
            n_frac += int((fin != np.floor(fin)).sum())
            n_zero += int((fin == 0).sum())
            n_min = min(n_min, float(fin.min()))
            n_max = max(n_max, float(fin.max()))
        lib.append(np.asarray(chunk.sum(axis=1)).ravel())
        det.append(np.diff(chunk.indptr))
    lib = np.concatenate(lib)
    det = np.concatenate(det)

    frame = pd.DataFrame({"context": contexts, "pert": perts})
    by_group = frame.groupby(["context", "pert"]).size()
    per_ctx_perts = frame.groupby("context")["pert"].apply(set)
    missing = {c: sorted(set(targets) - s) for c, s in per_ctx_perts.items()}
    extra = {c: sorted(s - set(targets)) for c, s in per_ctx_perts.items()}
    stats = {
        "shape": list(shape),
        "encoding": encoding,
        "dtype": dtype,
        "index_dtype": index_dtype,
        "indptr_dtype": indptr_dtype,
        "nnz": nnz,
        "density": nnz / (shape[0] * shape[1]),
        "sparsity": 1 - nnz / (shape[0] * shape[1]),
        "stored_value_min": n_min,
        "stored_value_max": n_max,
        "n_nan": n_nan,
        "n_inf": n_inf,
        "n_negative": n_neg,
        "n_fractional": n_frac,
        "n_explicit_zeros": n_zero,
        "library_size": describe(lib),
        "genes_detected": describe(det),
        "n_empty_cells": int((lib == 0).sum()),
        "cells_by_context": frame.context.value_counts().sort_index().to_dict(),
        "cells_per_perturbation": {
            "min": int(by_group.min()),
            "max": int(by_group.max()),
            "n_groups": int(len(by_group)),
        },
        "unique_perturbations_per_context": {c: len(s) for c, s in per_ctx_perts.items()},
        "missing_targets": missing,
        "extra_targets": extra,
        "non_targeting_cells": int((perts == "non-targeting").sum()),
        "n_genes": len(cand_genes),
        "gene_order_sha256": sha256_lines(cand_genes),
        "gene_order_matches_official": cand_genes == genes,
        "perturbation_labels_sha256": sha256_lines(perts),
        "context_labels_sha256": sha256_lines(contexts),
        "nnz_cap": 4_750_000_000,
        "nnz_fraction_of_cap": nnz / 4_750_000_000,
        "nnz_fraction_of_int32": nnz / 2**31,
    }
    print(
        json.dumps(
            {k: v for k, v in stats.items() if k not in ("missing_targets", "extra_targets")},
            indent=2,
            default=str,
        )
    )
    hard = {
        "no NaN": n_nan == 0,
        "no inf": n_inf == 0,
        "no negative": n_neg == 0,
        "no fractional": n_frac == 0,
        "no explicit stored zeros": n_zero == 0,
        "no empty cells": stats["n_empty_cells"] == 0,
        "no non-targeting cells": stats["non_targeting_cells"] == 0,
        "no missing target in any context": all(not v for v in missing.values()),
        "no extra target in any context": all(not v for v in extra.values()),
        "contexts exactly A/B/C, 120,000 each": stats["cells_by_context"]
        == {"A": 120_000, "B": 120_000, "C": 120_000},
        "gene order equals gene_names.csv": stats["gene_order_matches_official"],
        "nnz below CLI cap": nnz < 4_750_000_000,
    }
    for label, ok in hard.items():
        print(f"  [{'PASS' if ok else 'FAIL'}] {label}")
        if not ok:
            failures.append(label)

    rule("F1. THE THIRTEEN FROZEN LOCAL CHECKS (frozen inspect_bundle path)")
    report = bundle.inspect_bundle(CANDIDATE, genes)
    group_sizes = report.cells_per_group
    checks = {
        "n_cells == 360000": report.n_cells == 360_000,
        "n_genes == 18533": report.n_genes == 18_533,
        "gene order matches gene_names.csv": report.gene_order_matches,
        "raw integer counts": report.integer_valued,
        "non-negative": report.non_negative,
        "finite": report.finite,
        "exactly 400 cells per (context, perturbation)": bool((group_sizes == 400).all()),
        "300 perturbations in every context": all(
            v == 300 for v in report.n_perturbations_per_context.values()
        ),
        "contexts are exactly A, B, C": report.contexts == ("A", "B", "C"),
        "no control cells emitted": "non-targeting"
        not in set(group_sizes.index.get_level_values(1)),
        "under max_counts_per_cell (1,000,000)": report.max_counts_per_cell
        < generate.MAX_COUNTS_PER_CELL,
        "under max_nnz (4,750,000,000)": report.nnz < 4_750_000_000,
        "under max_cell_dim (400,000)": report.n_cells <= 400_000,
    }
    for label, ok in checks.items():
        print(f"  [{'PASS' if ok else 'FAIL'}] {label}")
        if not ok:
            failures.append(f"local check: {label}")
    frozen_checks = frozen_summary["checks"]
    if checks != frozen_checks:
        failures.append("local checks differ from the frozen dry-run checks")

    rule("F2. vcc prep --dry-run (must equal the frozen report)")
    cmd = [
        "vcc",
        "prep",
        str(CANDIDATE),
        "-g",
        str(CONTROLS / "gene_names.csv"),
        "--perts",
        str(CONTROLS / "pert_counts.csv"),
        "--dry-run",
        "--json",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    print("  $ " + " ".join(cmd))
    print(f"  exit code: {proc.returncode}")
    try:
        dry = json.loads(proc.stdout)
    except json.JSONDecodeError:
        dry = None
        failures.append("vcc prep --dry-run did not return JSON")
    frozen_dry = json.loads(FROZEN_DRY_RUN.read_text())
    same = dry == frozen_dry
    print(f"  identical to frozen dry-run report: {same}")
    if proc.returncode != 0:
        failures.append(f"vcc prep --dry-run exit code {proc.returncode}")
    if not same:
        failures.append("vcc prep --dry-run report differs from the frozen report")
    (OUTDIR / "vcc_prep_dry_run_recheck.json").write_text(proc.stdout)
    version = subprocess.run(["vcc", "--version"], capture_output=True, text=True).stdout.strip()

    audit = {
        "generated_utc": datetime.now(UTC).isoformat(),
        "vcc_version": version,
        "inputs": inputs,
        "control_checksums": control_hashes,
        "candidate": identity,
        "statistics": stats,
        "hard_checks": hard,
        "local_checks_13": checks,
        "local_checks_equal_frozen": checks == frozen_checks,
        "vcc_prep_dry_run": {
            "command": cmd,
            "exit_code": proc.returncode,
            "report": dry,
            "identical_to_frozen": same,
        },
        "failures": failures,
    }
    (OUTDIR / "candidate_audit.json").write_text(json.dumps(audit, indent=2, default=str) + "\n")
    rule("RESULT")
    print("  PASS" if not failures else "  FAIL:\n    " + "\n    ".join(failures))
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
