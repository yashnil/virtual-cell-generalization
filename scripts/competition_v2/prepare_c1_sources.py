"""Prepare the GREEN per-source statistics for C1 from raw public files (our code).

Jobs (each writes under ``outputs/competition_v2/c1_license_clean/sources/``):

* ``k562``  — K562 GWPS statistics + the public K562 fold's evaluation cells;
* ``h1``    — VCC 2025 H1 statistics (train + validation + test);
* ``h1cells`` — the public H1 fold's evaluation cells (train split; same seeded draw
  as ``scripts/competition_v2/run_atlasshift_public_h1.py``);
* ``cd4``   — CD4 genome-wide DE statistics.

X-Atlas is BLOCKED_PENDING_PERMISSION and has no job here. After each statistics job
the result is compared with the upstream AtlasShift file of the same name in
``outputs/competition_v2/atlasshift_c0/data`` (read only) and the comparison is written
to ``sources/equality_<job>.json``.

    uv run python scripts/competition_v2/prepare_c1_sources.py --job k562
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from virtual_cell.competition_v2 import sources  # noqa: E402

RAW = ROOT / "data" / "raw" / "competition_v2"
ARC = ROOT / "data" / "raw" / "arc2026" / "controls"
UPSTREAM = ROOT / "outputs" / "competition_v2" / "atlasshift_c0" / "data"
OUT = ROOT / "outputs" / "competition_v2" / "c1_license_clean" / "sources"
SEED = 20260910


def retained() -> list[str]:
    official = pd.read_csv(ARC / "pert_counts.csv").target_gene.astype(str).tolist()
    h1 = pd.read_csv(RAW / "pert_counts_Training.csv").target_gene.astype(str).tolist()
    return sources.retained_targets(official, h1)


def compare(ours: dict, upstream_path: Path) -> dict:
    """Array-by-array comparison with the upstream prepared file."""
    report = {"upstream": str(upstream_path.relative_to(ROOT)), "arrays": {}}
    with np.load(upstream_path, allow_pickle=False) as up:
        for key in up.files:
            a, b = np.asarray(ours[key]), up[key]
            entry = {"shape_equal": a.shape == b.shape}
            if a.shape == b.shape:
                if a.dtype.kind in "fc":
                    both = np.isfinite(a) & np.isfinite(b)
                    diff = np.abs(a.astype(np.float64) - b.astype(np.float64))[both]
                    scale = np.maximum(np.abs(b.astype(np.float64)[both]), 1e-30)
                    entry.update(
                        exact=bool(np.array_equal(a, b, equal_nan=True)),
                        nan_pattern_equal=bool(np.array_equal(np.isnan(a), np.isnan(b))),
                        max_abs=float(diff.max(initial=0)),
                        max_rel=float((diff / scale)[b[both] != 0].max(initial=0)),
                    )
                else:
                    entry["exact"] = bool(np.array_equal(a.astype(str), b.astype(str)))
            report["arrays"][key] = entry
    report["all_exact_or_close"] = all(
        e.get("shape_equal")
        and (e.get("exact") or e.get("max_rel", 1) < 1e-6 or e.get("max_abs", 1) < 1e-9)
        for e in report["arrays"].values()
    )
    return report


def save_cells(path: Path, counts: sparse.csr_matrix, labels, genes) -> None:
    counts = sparse.csr_matrix(counts, dtype=np.float32)
    np.savez_compressed(
        path,
        data=counts.data,
        indices=counts.indices,
        indptr=counts.indptr,
        shape=np.asarray(counts.shape),
        labels=np.asarray(labels).astype(str),
        genes=np.asarray(genes).astype(str),
    )


def job_k562() -> None:
    stats, cells = sources.prepare_k562(
        RAW / "K562_gwps_raw_singlecell_01.h5ad",
        retained(),
        seed=SEED,
        log=lambda m: print(m, flush=True),
    )
    np.savez_compressed(OUT / "K562_GWPS_CPM_full_statistics.npz", **stats)
    save_cells(OUT / "fold_cells_K562.npz", cells["counts"], cells["labels"], cells["genes"])
    eq = compare(stats, UPSTREAM / "K562_GWPS_CPM_full_statistics.npz")
    (OUT / "equality_k562.json").write_text(json.dumps(eq, indent=2))


def job_h1() -> None:
    files = {
        "train": "adata_Training.h5ad",
        "validation": "adata_Validation.h5ad",
        "test": "adata_Test.h5ad",
    }
    summaries = {
        split: sources.summarise_h1_split(RAW / f, log=lambda m: print(m, flush=True))
        for split, f in files.items()
    }
    stats = sources.assemble_h1(summaries)
    if len(stats["targets"]) != 300:
        raise ValueError("expected 300 public H1 targets")
    np.savez_compressed(OUT / "H1_2025_full_statistics.npz", **stats)
    eq = compare(stats, UPSTREAM / "H1_2025_full_statistics.npz")
    (OUT / "equality_h1.json").write_text(json.dumps(eq, indent=2))


def job_h1cells(max_real: int = 1000, n_controls: int = 8000) -> None:
    """Same draw as the C0 public-H1 harness: rng(SEED), per-label cap, sorted labels."""
    rng = np.random.default_rng(SEED)
    a = ad.read_h5ad(RAW / "adata_Training.h5ad", backed="r")
    labels = a.obs.target_gene.astype(str).to_numpy()
    picks = {}
    for t in np.unique(labels):
        rows = np.flatnonzero(labels == t)
        cap = n_controls if t == sources.CONTROL_LABEL else max_real
        if len(rows) > cap:
            rows = np.sort(rng.choice(rows, cap, replace=False))
        picks[t] = rows
    order = np.sort(np.concatenate(list(picks.values())))
    blocks = [sparse.csr_matrix(a.X[order[i : i + 20000]]) for i in range(0, len(order), 20000)]
    genes = a.var_names.astype(str).to_numpy()
    a.file.close()
    save_cells(OUT / "fold_cells_H1.npz", sparse.vstack(blocks, format="csr"), labels[order], genes)
    # the C0 harness continues this generator to split controls into reference / pool
    (OUT / "fold_cells_H1_rng_state.json").write_text(json.dumps(rng.bit_generator.state))


def job_cd4() -> None:
    stats = sources.prepare_cd4(RAW / "GWCD4i.DE_stats.h5ad", retained())
    np.savez_compressed(OUT / "CD4_DE_statistics.npz", **stats)
    eq = compare(stats, UPSTREAM / "CD4_DE_statistics.npz")
    (OUT / "equality_cd4.json").write_text(json.dumps(eq, indent=2))


def job_compare() -> None:
    """Re-run every equality comparison from the saved files."""
    for job, name in [
        ("k562", "K562_GWPS_CPM_full_statistics.npz"),
        ("h1", "H1_2025_full_statistics.npz"),
        ("cd4", "CD4_DE_statistics.npz"),
    ]:
        with np.load(OUT / name, allow_pickle=False) as d:
            eq = compare({k: d[k] for k in d.files}, UPSTREAM / name)
        (OUT / f"equality_{job}.json").write_text(json.dumps(eq, indent=2))
        print(job, eq["all_exact_or_close"], flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--job", choices=["k562", "h1", "h1cells", "cd4", "compare"], required=True)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    started = time.time()
    jobs = {
        "k562": job_k562,
        "h1": job_h1,
        "h1cells": job_h1cells,
        "cd4": job_cd4,
        "compare": job_compare,
    }
    jobs[args.job]()
    print(f"{args.job} done in {time.time() - started:.0f}s", flush=True)


if __name__ == "__main__":
    main()
