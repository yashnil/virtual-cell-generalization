"""Public held-out benchmark of C0 (AtlasShift) on H1 2025, with the promoter ablation.

Sections 12E and 13. The VCC 2025 H1 training split (150 targets, Arc's own CRISPRi
Perturb-seq) is held out: its perturbed cells are read **only** by the scorer. C0 is
rebuilt from the remaining AtlasShift sources — K562 GWPS (w 2), X-Atlas HCT116 (1),
X-Atlas HEK293T (1) and CD4 DE (0.5) — with every upstream constant unchanged. The H1
source is dropped because it *is* the held-out context. Nothing is tuned.

The upstream functions are imported verbatim from ``third_party/atlasshift/model.py``
and ``predict.py``; this file only wires them to a public target context and scores the
result with the local cell-eval2 reimplementation (``virtual_cell.arc.metrics``).
Run it with the pinned AtlasShift interpreter so upstream code sees its own pins:

    PYTHONPATH=src:scripts third_party/atlasshift/.venv/bin/python \
        scripts/competition_v2/run_atlasshift_public_h1.py

Design (mirrors Arc): H1 control cells are split disjointly into a template pool
(the "control bundle" C0 is built from) and a scoring reference; 400 generated cells
per target; real cells capped at ``MAX_REAL`` per target for memory.

Arms: ``C0`` (promoter prior on), ``C0_no_promoter``, ``G0_control`` (C0 generator,
zero effect), ``ANCHOR_mean_response`` (ORACLE local b: the held-out panel's true mean
response through the same generator) and ``ANCHOR_split_half`` (ORACLE local r).
"""

from __future__ import annotations

import gzip
import json
import re
import sys
import tempfile
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[2]
UPSTREAM = ROOT / "third_party" / "atlasshift"
sys.path.insert(0, str(UPSTREAM))
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

import model as up  # noqa: E402  (upstream AtlasShift, verbatim)
import predict as up_predict  # noqa: E402
import prepare as up_prepare  # noqa: E402

import run_count_generator_benchmark as frozen  # noqa: E402  (scorer glue, unchanged)
from virtual_cell.arc import metrics as arc_metrics  # noqa: E402

RAW = ROOT / "data" / "raw" / "competition_v2"
PREPARED = ROOT / "outputs" / "competition_v2" / "atlasshift_c0" / "data"
OUT = ROOT / "outputs" / "competition_v2" / "atlasshift_public_h1"
OFFICIAL_GENES = ROOT / "data" / "raw" / "arc2026" / "controls" / "gene_names.csv"
SOURCES = [
    "K562_GWPS_CPM_full_statistics.npz",
    "HCT116_full_statistics.npz",
    "HEK293T_full_statistics.npz",
]
WEIGHTS = [2, 1, 1]  # upstream WEIGHTS with the held-out H1 entry removed
CELLS = up_predict.CELLS
SEED = up_predict.SEED
MAX_REAL = 1000
N_REFERENCE_CONTROLS = 4000
N_POOL_CONTROLS = 4000


def load_h1(path: Path, genes: np.ndarray, rng: np.random.Generator):
    """Real H1 cells restricted to ``genes``: controls and <= MAX_REAL per target."""
    a = ad.read_h5ad(path, backed="r")
    labels = a.obs.target_gene.astype(str).to_numpy()
    col = pd.Index(a.var_names.astype(str)).get_indexer(genes)
    if (col < 0).any():
        raise ValueError("evaluation genes missing from H1")
    picks = {}
    for t in np.unique(labels):
        rows = np.flatnonzero(labels == t)
        cap = N_REFERENCE_CONTROLS + N_POOL_CONTROLS if t == "non-targeting" else MAX_REAL
        if len(rows) > cap:
            rows = np.sort(rng.choice(rows, cap, replace=False))
        picks[t] = rows
    order = np.sort(np.concatenate(list(picks.values())))
    blocks = []
    for left in range(0, len(order), 20000):
        x = sparse.csr_matrix(a.X[order[left : left + 20000]])
        blocks.append(x[:, col])
    x = sparse.vstack(blocks, format="csr")
    a.file.close()
    pos = {r: i for i, r in enumerate(order)}
    return {t: x[[pos[r] for r in rows]] for t, rows in picks.items()}


def per_pert(pred, real, tg):
    n = real.lfc.shape[0]
    out = {k: np.zeros(n) for k in ("n_real", "n_pred", "k")}
    for p in range(n):
        rs = arc_metrics._drop_target(real.significant[p], tg[p])
        ps = arc_metrics._drop_target(pred.significant[p] & real.adjudicable[p], tg[p])
        out["n_real"][p], out["n_pred"][p] = rs.sum(), ps.sum()
        out["k"][p] = (ps & (np.sign(pred.lfc[p]) == np.sign(real.lfc[p]))).sum()
    with np.errstate(invalid="ignore", divide="ignore"):
        out["precision"] = out["k"] / out["n_pred"]
        out["yield"] = np.minimum(1, out["n_pred"] / out["n_real"])
    out["fid"] = arc_metrics.direction_fidelity_yield(pred, real, target_gene=tg)
    out["reach"] = arc_metrics.direction_reach(pred, real, target_gene=tg)
    out["jaccard"] = arc_metrics.sig_jaccard(pred, real, target_gene=tg)
    out["nmae"] = arc_metrics.lfc_nmae(pred, real, target_gene=tg)
    return out


def main() -> None:
    global OUT, SOURCES, WEIGHTS
    smoke = "--smoke" in sys.argv  # wiring check only: 10 targets, no K562
    if smoke:
        OUT = OUT.with_name(OUT.name + "_smoke")
        SOURCES, WEIGHTS = SOURCES[1:], WEIGHTS[1:]
    started = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    h1_genes = pd.read_csv(RAW / "gene_names.csv").iloc[:, 0].astype(str).to_numpy()
    official = pd.read_csv(OFFICIAL_GENES).gene_name.astype(str).to_numpy()
    genes = official[np.isin(official, h1_genes)]
    cells = load_h1(RAW / "adata_Training.h5ad", genes, rng)
    targets = np.array(sorted(t for t in cells if t != "non-targeting"))
    if smoke:
        targets = targets[:10]

    ctrl = cells.pop("non-targeting")
    perm = rng.permutation(ctrl.shape[0])
    reference = ctrl[perm[:N_REFERENCE_CONTROLS]].toarray().astype(np.int64)
    pool = ctrl[perm[N_REFERENCE_CONTROLS:]]
    real_blocks = [cells[t].toarray().astype(np.int32) for t in targets]

    with tempfile.TemporaryDirectory(dir=OUT) as tmp:
        pool_path = Path(tmp) / "pool.h5ad"
        ad.AnnData(pool.astype(np.float32), var=pd.DataFrame(index=genes)).write_h5ad(pool_path)
        template, depths, mean, bulk = up_predict.control_template(pool_path, genes, CELLS, 4, SEED)

    # ---- C0 effects, upstream functions and constants -------------------------
    sources = [up.load_xatlas(PREPARED / s, prior_counts=100000) for s in SOURCES]
    cd4, available = up.aligned_cd4(
        PREPARED / "CD4_DE_statistics.npz", targets, genes, common_subtract=1, center_scope="source"
    )
    base = {
        space: up.fuse_source_centered(
            sources, WEIGHTS, targets, genes, space=space, common_subtract=1
        )
        for space in ["log2fc", "bulk_delta"]
    }
    coverage_by_source = {
        Path(s).stem: np.isin(targets, src.targets[src.n_cells >= 20])
        for s, src in zip(SOURCES, sources, strict=True)
    }
    coverage_by_source["CD4"] = available.any(axis=1)
    n_sources = np.sum(list(coverage_by_source.values()), axis=0)
    del sources

    gtf_rows = []
    with gzip.open(RAW / "gencode.v47.annotation.gtf.gz", "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] != "gene":
                continue
            attrs = dict(re.findall('(\\w+) "([^"]*)"', f[8]))
            if "gene_name" in attrs:
                gtf_rows.append(
                    {
                        "gene": attrs["gene_name"],
                        "chromosome": f[0],
                        "strand": f[6],
                        "tss": int(f[3] if f[6] == "+" else f[4]),
                    }
                )
    pairs_path = OUT / "h1_panel_pairs.csv"
    up_prepare.pairs(pd.DataFrame(gtf_rows), targets, genes).to_csv(pairs_path, index=False)

    def desired(promoter: bool, zero: bool = False):
        out = []
        for space, control, amplitude in [("log2fc", mean, 0.6), ("bulk_delta", bulk, 0.3)]:
            effect, cov = base[space]
            effect, _ = up.add_cd4_family(
                effect, cov, cd4, available, 0.5, space=space, control_probability=control
            )
            if zero:
                effect = np.zeros_like(effect)
            prob = up.desired_mean(control, effect, space=space, amplitude=amplitude, clip=3)
            changed = []
            if promoter:
                prob, changed = up.apply_promoter_prior(
                    prob, control, targets, genes, pairs_path, 0.15
                )
            out.append(prob)
        return out, changed

    # oracle mean response of the held-out panel, as a two-moment target
    real_cpm = np.stack([(b / b.sum(axis=1, keepdims=True)).mean(axis=0) for b in real_blocks])
    real_bulk = np.stack([b.sum(axis=0) / b.sum() for b in real_blocks])
    mr_moments = [
        np.tile(real_cpm.mean(axis=0), (len(targets), 1)),
        np.tile(real_bulk.mean(axis=0), (len(targets), 1)),
    ]

    # ---- reference side ------------------------------------------------------
    tg = np.array([np.flatnonzero(genes == t)[0] if t in genes else -1 for t in targets])
    exclude = np.isin(genes, targets)
    keep = ~exclude
    ctrl_profile = arc_metrics.bulk_profile(reference)
    ctrl_disp = arc_metrics.jackknife_dispersion(reference)
    real_profiles = np.stack([arc_metrics.bulk_profile(b) for b in real_blocks])
    real_disp = np.array([arc_metrics.jackknife_dispersion(b) for b in real_blocks])
    real_de = arc_metrics.de_table(real_blocks, reference)
    # DE tables live on the control-gated *tested* axis, not the full gene axis, so the
    # target index must be re-expressed there (the frozen v1 harness passes the full-axis
    # index to both; see reports/competition_v2/competitive_baseline_expansion_v1.md).
    tested_pos = np.cumsum(real_de.tested) - 1
    tg_de = np.array([tested_pos[i] if i >= 0 and real_de.tested[i] else -1 for i in tg])
    true_norm = np.linalg.norm((real_profiles - ctrl_profile)[:, keep], axis=1)

    rows, frames, structure = [], [], []

    def score(
        name,
        blocks,
        *,
        real=real_de,
        cprof=ctrl_profile,
        cdisp=ctrl_disp,
        rprof=real_profiles,
        rdisp=real_disp,
        ref_ctrl=reference,
    ):
        profiles = np.stack([arc_metrics.bulk_profile(b) for b in blocks])
        disp = np.array([arc_metrics.jackknife_dispersion(b) for b in blocks])
        de = frozen.merge_de(
            [arc_metrics.de_table([b], ref_ctrl, tested=real.tested) for b in blocks], real.tested
        )
        mse = arc_metrics.expr_mse_unbiased_capped(
            profiles,
            rprof,
            cprof,
            pred_dispersion=disp,
            real_dispersion=rdisp,
            ctrl_dispersion=cdisp,
            target_gene=tg,
        )
        sc = {
            "pds_cosine": float(
                np.mean(
                    arc_metrics.pds_cosine(
                        profiles - cprof[None, :], rprof - cprof[None, :], exclude=exclude
                    )
                )
            ),
            "expr_mse_unbiased_capped_norm": float(mse.value),
            "de_wilcoxon_direction_fidelity_yield_raw": float(
                np.nanmean(arc_metrics.direction_fidelity_yield(de, real, target_gene=tg_de))
            ),
            "de_wilcoxon_direction_reach_raw": float(
                np.nanmean(arc_metrics.direction_reach(de, real, target_gene=tg_de))
            ),
            "de_wilcoxon_sig_jaccard": float(
                np.nanmean(arc_metrics.sig_jaccard(de, real, target_gene=tg_de))
            ),
            "de_wilcoxon_lfc_nmae": float(
                np.nanmean(arc_metrics.lfc_nmae(de, real, target_gene=tg_de))
            ),
        }
        rows.append({"arm": name, **sc})
        f = pd.DataFrame(per_pert(de, real, tg_de))
        f.insert(0, "perturbation", targets)
        f.insert(1, "arm", name)
        f["n_sources"] = n_sources
        f["pred_norm"] = np.linalg.norm((profiles - cprof)[:, keep], axis=1)
        f["true_norm"] = true_norm
        frames.append(f)
        structure.append(frozen.structure_stats(blocks, name))
        print(
            f"  {name:22s} PDS={sc['pds_cosine']:.3f} "
            f"FID={sc['de_wilcoxon_direction_fidelity_yield_raw']:.3f} "
            f"({time.time() - started:.0f}s)",
            flush=True,
        )

    promoter_changes = {}
    for name, (moments, changed) in {
        "C0": desired(True),
        "C0_no_promoter": desired(False),
        "G0_control": desired(False, zero=True),
        "ANCHOR_mean_response": (mr_moments, []),
    }.items():
        promoter_changes[name] = len(changed)
        blocks = [
            up.dual_moment_counts(
                template,
                moments[0][i] / moments[0][i].sum(),
                moments[1][i] / moments[1][i].sum(),
                depths=depths,
                seed=up.seed_for("H1:" + t, SEED),
            )
            for i, t in enumerate(targets)
        ]
        score(name, blocks)

    half_rng = np.random.default_rng(SEED + 2)
    rp = half_rng.permutation(len(reference))
    ca, cb = reference[rp[: len(rp) // 2]], reference[rp[len(rp) // 2 :]]
    ha, hb = [], []
    for b in real_blocks:
        i = half_rng.permutation(len(b))
        ha.append(b[i[: len(b) // 2]])
        hb.append(b[i[len(b) // 2 :]])
    score(
        "ANCHOR_split_half",
        ha,
        real=arc_metrics.de_table(hb, cb),
        cprof=arc_metrics.bulk_profile(cb),
        cdisp=arc_metrics.jackknife_dispersion(cb),
        rprof=np.stack([arc_metrics.bulk_profile(b) for b in hb]),
        rdisp=np.array([arc_metrics.jackknife_dispersion(b) for b in hb]),
        ref_ctrl=ca,
    )
    structure.append(frozen.structure_stats(real_blocks, "real_perturbed"))
    structure.append(frozen.structure_stats([reference], "real_control_reference"))

    pd.DataFrame(rows).to_csv(OUT / "scores_raw.csv", index=False)
    pd.concat(frames).to_csv(OUT / "per_perturbation.csv", index=False)
    pd.DataFrame(structure).to_csv(OUT / "structure.csv", index=False)
    (OUT / "summary.json").write_text(
        json.dumps(
            {
                "held_out": "H1 2025 train split (150 targets)",
                "sources": SOURCES + ["CD4_DE_statistics.npz"],
                "weights": WEIGHTS + [0.5],
                "n_genes": int(len(genes)),
                "n_targets": int(len(targets)),
                "targets_with_any_source": int((n_sources > 0).sum()),
                "coverage_by_source": {k: int(v.sum()) for k, v in coverage_by_source.items()},
                "promoter_pairs_on_panel": int(len(pd.read_csv(pairs_path))),
                "promoter_entries_capped": promoter_changes,
                "real_cells_per_target_cap": MAX_REAL,
                "reference_controls": N_REFERENCE_CONTROLS,
                "elapsed_seconds": round(time.time() - started, 1),
            },
            indent=2,
        )
    )
    print(f"done in {time.time() - started:.0f}s")


if __name__ == "__main__":
    main()
