"""Build, validate and package the license-clean C1 Arc candidate. NEVER submits.

Stages (run in order; each is idempotent):

* ``emit``     — GREEN sources only (K562, H1 2025, CD4) + GENCODE promoter cap, the C1
  variant chosen by ``outputs/competition_v2/c1_license_clean/decision.json``;
  300 targets x 3 contexts x 400 cells -> ``prediction.h5ad``;
* ``package``  — lossless compaction and ``vcc.prep`` packaging with the vendored
  AtlasShift packaging utilities (``compact.py`` / ``pack.py``, run unmodified in their
  own pinned env; they reorder and compress only, no modelling);
* ``validate`` — the frozen local checks + ``vcc prep --dry-run --json`` + manifest.

    uv run python scripts/competition_v2/build_c1_candidate.py --stage emit
"""

from __future__ import annotations

import os

os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import argparse  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import multiprocessing as mp  # noqa: E402
import subprocess  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import anndata as ad  # noqa: E402
import h5py  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import sparse  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from virtual_cell.arc import bundle, generate  # noqa: E402
from virtual_cell.competition_v2 import fusion, generator, licensing, sources  # noqa: E402

BASE = ROOT / "outputs" / "competition_v2" / "c1_license_clean"
OUT = BASE  # candidate files live beside the fold results
SRC = BASE / "sources"
CONTROLS = ROOT / "data" / "raw" / "arc2026" / "controls"
RAW = ROOT / "data" / "raw" / "competition_v2"
UPSTREAM = ROOT / "third_party" / "atlasshift"
PRED = OUT / "prediction.h5ad"
COMPACT = OUT / "prediction_compact.h5ad"
PACKAGE = OUT / "c1_license_clean_val.vcc"
CELLS = 400
C1_SOURCES = ["K562", "H1", "CD4"]

_G: dict = {}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(16 * 1024**2), b""):
            h.update(block)
    return h.hexdigest()


def _emit_one(i: int):
    g = _G
    counts = generator.dual_moment_counts(
        g["template"],
        g["p_cpm"][i],
        g["p_bulk"][i],
        depths=g["depths"],
        seed=generator.seed_for(f"{g['ctx']}:{g['targets'][i]}"),
    )
    return sparse.csr_matrix(counts.astype(np.float32))


def stage_emit(jobs: int) -> None:
    if PRED.exists():
        raise FileExistsError(PRED)
    decision = json.loads((BASE / "decision.json").read_text())
    variant = decision["selected_variant"]
    licensing.assert_sources_allowed(C1_SOURCES + ["GENCODE"])
    targets = pd.read_csv(CONTROLS / "pert_counts.csv").target_gene.astype(str).to_numpy()
    genes = pd.read_csv(CONTROLS / "gene_names.csv").gene_name.astype(str).to_numpy()
    srcs = [
        sources.load_source(SRC / "K562_GWPS_CPM_full_statistics.npz"),
        sources.load_source(SRC / "H1_2025_full_statistics.npz"),
    ]
    with np.load(SRC / "CD4_DE_statistics.npz") as d:
        cd4 = {k: d[k] for k in d.files}
    pairs = generator.promoter_pairs(
        generator.gencode_tss(RAW / "gencode.v47.annotation.gtf.gz"), targets, genes
    )
    exclude = np.isin(genes, targets)
    agree, n_src = fusion.source_agreement(
        fusion.panel_source_vectors(srcs, cd4, targets, genes), exclude
    )
    if variant == "C1a":
        shrink = None
    elif variant.startswith("C1b_S2_f"):
        shrink = fusion.shrinkage_factor(
            "S2", len(targets), floor=float(variant[8:]), agreement=agree
        )
    else:
        raise ValueError(variant)
    pending = PRED.with_name("prediction.partial.h5ad")
    pending.unlink(missing_ok=True)
    writer = generator.CountWriter(pending, targets, genes, list("ABC"), CELLS)
    started = time.time()
    expected = {}
    for ci, ctx in enumerate("ABC"):
        ctrl = ad.read_h5ad(CONTROLS / f"context_{ctx}.h5ad")
        if not np.array_equal(ctrl.var_names.astype(str), genes):
            raise ValueError("control gene order differs from gene_names.csv")
        template, depths, mean, bulk = generator.control_template(
            ctrl.X, CELLS, generator.POOL_K, generator.SEED + ci
        )
        controls = {"log2fc": mean, "bulk_delta": bulk}
        effects = fusion.fused_effects(srcs, [1.0, 1.0], cd4, 1.0, targets, genes, controls)
        p_cpm, p_bulk = generator.expected_moments(
            effects, controls, targets, genes, shrink=shrink, pairs=pairs
        )
        expected[ctx] = (p_cpm, p_bulk, bulk)
        _G.clear()
        _G.update(
            template=template, depths=depths, p_cpm=p_cpm, p_bulk=p_bulk, ctx=ctx, targets=targets
        )
        with mp.get_context("fork").Pool(jobs) as pool:
            for k, block in enumerate(pool.imap(_emit_one, range(len(targets)), chunksize=2)):
                writer.append(block)
                if (k + 1) % 50 == 0:
                    print(f"{ctx}: {k + 1}/300 ({time.time() - started:.0f}s)", flush=True)
    writer.close()
    pending.rename(PRED)
    np.savez_compressed(
        OUT / "c1_expected_moments.npz",
        targets=targets,
        genes=genes,
        agreement=agree,
        n_sources=n_src,
        shrink=np.ones(len(targets)) if shrink is None else shrink,
        **{
            f"{c}_{k}": v
            for c, (a, b, bl) in expected.items()
            for k, v in [("p_cpm", a), ("p_bulk", b), ("ctrl_bulk", bl)]
        },
    )
    pairs.to_csv(OUT / "c1_promoter_pairs.csv", index=False)
    (OUT / "emit.json").write_text(
        json.dumps(
            {
                "variant": variant,
                "sources": C1_SOURCES,
                "weights": "equal (1 each, CD4 1)",
                "promoter_pairs": int(len(pairs)),
                "seconds": round(time.time() - started, 1),
                "nnz": int(writer.nnz),
            },
            indent=2,
        )
    )
    print(f"saved {PRED.name}: {writer.nobs} cells, {writer.nnz} nnz", flush=True)


def stage_package() -> None:
    py = UPSTREAM / ".venv" / "bin" / "python"
    scratch = OUT / "scratch"
    scratch.mkdir(exist_ok=True)
    if not COMPACT.exists():
        subprocess.run(
            [str(py), str(UPSTREAM / "compact.py"), str(PRED), str(COMPACT)], check=True, cwd=OUT
        )
    if not PACKAGE.exists():
        data = OUT / "pack_data"
        data.mkdir(exist_ok=True)
        for f in ["gene_names.csv", "pert_counts.csv"]:
            link = data / f
            if not link.exists():
                link.symlink_to(CONTROLS / f)
        subprocess.run(
            [
                str(py),
                str(UPSTREAM / "pack.py"),
                str(COMPACT),
                "--data-dir",
                str(data),
                "--output",
                str(PACKAGE),
                "--scratch-dir",
                str(scratch),
            ],
            check=True,
            cwd=OUT,
        )


def stage_validate() -> int:
    genes = pd.read_csv(CONTROLS / "gene_names.csv").gene_name.astype(str).tolist()
    panel = set(pd.read_csv(CONTROLS / "pert_counts.csv").target_gene.astype(str))
    report = bundle.inspect_bundle(PRED, genes)
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
        "target set equals pert_counts.csv": set(sizes.index.get_level_values(1)) == panel,
    }
    with h5py.File(PRED, "r") as f:
        data = f["X/data"]
        zeros = sum(
            int((data[i : i + (1 << 26)] == 0).sum()) for i in range(0, data.shape[0], 1 << 26)
        )
    checks["no explicit stored zeros"] = zeros == 0
    dry = subprocess.run(
        [
            "vcc",
            "prep",
            str(COMPACT),
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
    (OUT / "vcc_prep_dry_run.json").write_text(json.dumps(dry_report, indent=2) + "\n")
    decision = json.loads((BASE / "decision.json").read_text())
    head = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True
    ).stdout.strip()
    manifest = {
        "candidate": "C1_LICENSE_CLEAN",
        "selected_variant": decision["selected_variant"],
        "submitted": False,
        "provenance": "ours: reimplemented atlas backbone (AtlasShift concepts, MIT, "
        "attributed) on GREEN sources only",
        "partition_panel": "val / vcc2026-val-1",
        "sources": {n: licensing.STATUS[n] for n in C1_SOURCES + ["GENCODE"]},
        "excluded_sources": {n: licensing.STATUS[n] for n in ["HCT116", "HEK293T", "KADEN_RPE1"]},
        "prediction": {
            "path": str(PRED.relative_to(ROOT)),
            "sha256": sha256(PRED),
            "bytes": PRED.stat().st_size,
        },
        "compact": {
            "path": str(COMPACT.relative_to(ROOT)),
            "sha256": sha256(COMPACT),
            "bytes": COMPACT.stat().st_size,
        },
        "vcc_package": {
            "path": str(PACKAGE.relative_to(ROOT)),
            "sha256": sha256(PACKAGE),
            "bytes": PACKAGE.stat().st_size,
        }
        if PACKAGE.exists()
        else None,
        "source_statistics_sha256": {
            p.name: sha256(p) for p in sorted(SRC.glob("*_statistics.npz"))
        },
        "promoter_pairs_sha256": sha256(OUT / "c1_promoter_pairs.csv"),
        "repo_git_head": head,
        "constants": {
            "weights": "equal: K562 1, H1 1, CD4 1",
            "amplitude": "log2fc 0.6, bulk_delta 0.3",
            "clip": 3,
            "prior_counts": 100000,
            "minimum_cells": 20,
            "promoter_fraction": 0.15,
            "seed": generator.SEED,
        },
        "bundle": {**report.as_dict(), "explicit_zeros": zeros},
        "checks": checks,
        "all_checks_pass": all(checks.values()),
        "vcc_prep_dry_run": dry_report,
    }
    (OUT / "c1_manifest.json").write_text(json.dumps(manifest, indent=2, default=float) + "\n")
    for k, v in checks.items():
        print(f"  [{'PASS' if v else 'FAIL'}] {k}")
    return 0 if manifest["all_checks_pass"] else 1


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=["emit", "package", "validate"], required=True)
    parser.add_argument("--jobs", type=int, default=9)
    args = parser.parse_args()
    if args.stage == "emit":
        stage_emit(args.jobs)
    elif args.stage == "package":
        stage_package()
    else:
        return stage_validate()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
