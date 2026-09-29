"""Build, validate and package the C2 Arc candidate, only if C2 passed rule J. NEVER submits.

The mean predictor is C1's, unchanged: GREEN sources, equal weights and the promoter cap.
The candidate arm in ``outputs/competition_v2/c2_calibration/c2_decision.json`` fixes the
generator, the global amplitude ``a`` and, optionally, the per-target exponent ``beta``.

On Arc, the generator's control pool is the context's full 18,400-cell controls bundle
(there is no reference split: the hidden scorer holds the reference).

    uv run python scripts/competition_v2/build_c2_candidate.py --stage emit|package|validate
"""

from __future__ import annotations

import os

os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import argparse  # noqa: E402
import json  # noqa: E402
import multiprocessing as mp  # noqa: E402
import re  # noqa: E402
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
sys.path.insert(0, str(Path(__file__).resolve().parent))

import build_c1_candidate as b1  # noqa: E402
import run_c1_public_folds as c1f  # noqa: E402

from virtual_cell.arc import bundle, generate  # noqa: E402
from virtual_cell.competition_v2 import fusion, generator, licensing, realisation  # noqa: E402

CAL = ROOT / "outputs" / "competition_v2" / "c2_calibration"
OUT = ROOT / "outputs" / "competition_v2" / "c2_candidate"
PRED = OUT / "prediction.h5ad"
COMPACT = OUT / "prediction_compact.h5ad"
PACKAGE = OUT / "c2_val.vcc"
CONTROLS = b1.CONTROLS
CELLS = 400
N_DONORS = 400

_G: dict = {}


def parse_arm(arm: str) -> tuple[str, float, float | None]:
    m = re.fullmatch(r"(G[0-9a-z]+)_a([0-9.]+)(?:_beta([+-][0-9.]+))?", arm)
    if not m:
        raise ValueError(f"unexpected candidate arm {arm!r}")
    return m.group(1), float(m.group(2)), float(m.group(3)) if m.group(3) else None


def _emit_one(i: int):
    g = _G
    t = g["targets"][i]
    gen = g["gen"]
    if gen == "G0":
        counts = generator.dual_moment_counts(
            g["template"],
            g["p_cpm"][i],
            g["p_bulk"][i],
            depths=g["depths"],
            seed=generator.seed_for(f"{g['ctx']}:{t}"),
        )
    else:
        rng = np.random.default_rng(generator.seed_for(f"{g['ctx']}:{gen}:{t}"))
        if gen == "G1ci":
            drng = np.random.default_rng(generator.seed_for(f"{g['ctx']}:donors:{t}"))
            pick = np.sort(drng.choice(g["pool"].shape[0], N_DONORS, replace=False))
            donors = g["pool"][pick].toarray()
        else:
            donors = g["donors"]
        if gen in ("G1c", "G1ci"):
            ratio = realisation.effect_ratio(g["p_cpm"][i], g["c_cpm"])
            counts = realisation.g1_counts(donors, ratio, g["c_cpm"], rng=rng)
        elif gen == "G1b":
            ratio = realisation.effect_ratio(g["p_bulk"][i], g["c_bulk"])
            counts = realisation.g1_counts(donors, ratio, g["c_bulk"], rng=rng)
        elif gen == "G2":
            counts = realisation.g2_counts(donors.sum(1), g["p_bulk"][i], rng=rng)
        elif gen == "G3":
            counts = realisation.g3_counts(donors.sum(1), g["p_bulk"][i], g["phi"], rng=rng)
        else:
            raise ValueError(gen)
    return sparse.csr_matrix(counts.astype(np.float32))


def stage_emit(jobs: int) -> None:
    decision = json.loads((CAL / "c2_decision.json").read_text())
    if not decision["c2_pass"]:
        raise SystemExit("C2 did not pass rule J; nothing is built")
    arm = decision["c2_candidate"]
    gen, a, beta = parse_arm(arm)
    OUT.mkdir(parents=True, exist_ok=True)
    if PRED.exists():
        raise FileExistsError(PRED)
    licensing.assert_sources_allowed(b1.C1_SOURCES + ["GENCODE"])
    targets = pd.read_csv(CONTROLS / "pert_counts.csv").target_gene.astype(str).to_numpy()
    genes = pd.read_csv(CONTROLS / "gene_names.csv").gene_name.astype(str).to_numpy()
    green = c1f.load_green()
    srcs = [green["K562"], green["H1"]]
    cd4 = c1f.load_cd4()
    pairs = generator.promoter_pairs(
        generator.gencode_tss(b1.RAW / "gencode.v47.annotation.gtf.gz"), targets, genes
    )
    k = c1f.coverage(["K562", "H1"], green, cd4, targets, genes).sum(1).to_numpy().astype(float)
    if beta is None:
        shrink = None if a == 1.0 else np.full(len(targets), a)
    else:
        kbar = k[k > 0].mean()
        shrink = np.where(k > 0, a * (k / kbar) ** beta, a)
    pending = PRED.with_name("prediction.partial.h5ad")
    pending.unlink(missing_ok=True)
    writer = generator.CountWriter(pending, targets, genes, list("ABC"), CELLS)
    started = time.time()
    for ci, ctx in enumerate("ABC"):
        ctrl = ad.read_h5ad(CONTROLS / f"context_{ctx}.h5ad")
        if not np.array_equal(ctrl.var_names.astype(str), genes):
            raise ValueError("control gene order differs from gene_names.csv")
        pool = sparse.csr_matrix(ctrl.X)
        template, depths, mean, bulk = generator.control_template(
            pool, CELLS, generator.POOL_K, generator.SEED + ci
        )
        controls = {"log2fc": mean, "bulk_delta": bulk}
        effects = fusion.fused_effects(srcs, [1.0, 1.0], cd4, 1.0, targets, genes, controls)
        p_cpm, p_bulk = generator.expected_moments(
            effects, controls, targets, genes, shrink=shrink, pairs=pairs
        )
        rng = np.random.default_rng(generator.SEED + 7 + ci)
        donors = pool[np.sort(rng.choice(pool.shape[0], N_DONORS, replace=False))].toarray()
        phi = (
            realisation.fit_dispersion(
                pool[np.sort(rng.choice(pool.shape[0], 3000, replace=False))].toarray(), bulk
            )
            if gen == "G3"
            else None
        )
        _G.clear()
        _G.update(
            gen=gen,
            template=template,
            depths=depths,
            p_cpm=p_cpm / p_cpm.sum(1, keepdims=True),
            p_bulk=p_bulk / p_bulk.sum(1, keepdims=True),
            ctx=ctx,
            targets=targets,
            pool=pool,
            donors=donors,
            c_cpm=mean,
            c_bulk=bulk,
            phi=phi,
        )
        with mp.get_context("fork").Pool(jobs) as p:
            for j, block in enumerate(p.imap(_emit_one, range(len(targets)), chunksize=2)):
                writer.append(block)
                if (j + 1) % 50 == 0:
                    print(f"{ctx}: {j + 1}/300 ({time.time() - started:.0f}s)", flush=True)
    writer.close()
    pending.rename(PRED)
    (OUT / "emit.json").write_text(
        json.dumps(
            {
                "candidate_arm": arm,
                "generator": gen,
                "amplitude": a,
                "beta": beta,
                "sources": b1.C1_SOURCES,
                "promoter_pairs": int(len(pairs)),
                "seconds": round(time.time() - started, 1),
                "nnz": int(writer.nnz),
            },
            indent=2,
        )
    )


def stage_package() -> None:
    for name, value in {"OUT": OUT, "PRED": PRED, "COMPACT": COMPACT, "PACKAGE": PACKAGE}.items():
        setattr(b1, name, value)
    b1.stage_package()


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
    emit = json.loads((OUT / "emit.json").read_text())
    head = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True
    ).stdout.strip()
    manifest = {
        "candidate": "C2",
        "arm": emit["candidate_arm"],
        "submitted": False,
        "mean_predictor": "C1a, unchanged (GREEN sources only)",
        "sources": {n: licensing.STATUS[n] for n in b1.C1_SOURCES + ["GENCODE"]},
        "excluded_sources": {n: licensing.STATUS[n] for n in ["HCT116", "HEK293T", "KADEN_RPE1"]},
        "emit": emit,
        "prediction": {"path": str(PRED.relative_to(ROOT)), "sha256": b1.sha256(PRED)},
        "compact": {"path": str(COMPACT.relative_to(ROOT)), "sha256": b1.sha256(COMPACT)},
        "vcc_package": {
            "path": str(PACKAGE.relative_to(ROOT)),
            "sha256": b1.sha256(PACKAGE),
            "bytes": PACKAGE.stat().st_size,
        }
        if PACKAGE.exists()
        else None,
        "repo_git_head": head,
        "bundle": {**report.as_dict(), "explicit_zeros": zeros},
        "checks": checks,
        "all_checks_pass": all(checks.values()),
        "vcc_prep_dry_run": dry_report,
    }
    (OUT / "c2_manifest.json").write_text(json.dumps(manifest, indent=2, default=float) + "\n")
    for key, ok in checks.items():
        print(f"  [{'PASS' if ok else 'FAIL'}] {key}")
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
