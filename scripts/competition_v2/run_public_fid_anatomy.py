"""Public FID / yield anatomy of the frozen V1 model, plus the predeclared amplitude ladder.

Sections 16 and 17 of the competitive-baseline-expansion phase. Diagnostic only:
nothing here is selected, tuned, or fed back into any frozen artifact.

For one held-out public context (leave-one-context-out over the frozen four-context
design), the frozen V1 predictor ``delta = m_hat + w[tier] beta`` (M3b_basal_shrunk,
w = 0.50 / 0.25 / 0) is generated into counts with the frozen G1 transport generator
and scored against the context's **real** held-out cells with the local cell-eval2
reimplementation. The data path, control split, panel selection and scorer are
imported unchanged from ``scripts/run_count_generator_benchmark.py``.

Arms (each 300 perturbations x 400 cells):

* ``T2`` / ``T1`` / ``T0`` — every panel perturbation predicted as that tier, a = 1
* ``MIX_a{a}`` — Arc prevalence (7 / 79 / 214, seeded assignment) with the whole
  predicted effect scaled by the predeclared ladder a in {0.5, 0.75, 1, 1.25, 1.5, 2, 3}
* ``ANCHOR_mean_response`` — ORACLE local baseline b: the held-out context's own mean
  perturbation response, identical for every perturbation (reads held-out truth; an
  anchor, never a prediction)
* ``ANCHOR_split_half`` — ORACLE local replicate r: one disjoint half of the real
  cells (with its own control half) scored against the other half (one split, not five)

Per perturbation it records n_pred, n_real, k, directional precision k / n_pred,
yield min(1, n_pred / n_real), FID, REACH, Jaccard, NMAE and effect norms.

Reproduce: ``uv run python scripts/competition_v2/run_public_fid_anatomy.py --context <ctx>``
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import run_count_generator_benchmark as frozen  # noqa: E402
from virtual_cell.arc import bundle, generate  # noqa: E402
from virtual_cell.arc import metrics as arc_metrics  # noqa: E402
from virtual_cell.arc.panel import build_panel_map  # noqa: E402
from virtual_cell.modelling import mean_response as mr  # noqa: E402

OUTDIR = ROOT / "outputs" / "competition_v2" / "public_fid_anatomy"
AMPLITUDES = (0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0)  # predeclared in the phase brief
TIER_WEIGHTS = {2: 0.5, 1: 0.25, 0: 0.0}  # frozen v1
ARC_PREVALENCE = {2: 7, 1: 79, 0: 214}
SEED = 20260925


def per_pert_de(pred: arc_metrics.DETable, real: arc_metrics.DETable, tg: np.ndarray) -> dict:
    """Every quantity of FID's factorisation, per perturbation."""
    n = real.lfc.shape[0]
    rows = {k: np.full(n, np.nan) for k in ("n_real", "n_pred", "k")}
    for p in range(n):
        real_sig = arc_metrics._drop_target(real.significant[p], tg[p])
        pred_sig = arc_metrics._drop_target(pred.significant[p] & real.adjudicable[p], tg[p])
        same = np.sign(pred.lfc[p]) == np.sign(real.lfc[p])
        rows["n_real"][p] = real_sig.sum()
        rows["n_pred"][p] = pred_sig.sum()
        rows["k"][p] = (pred_sig & same).sum()
    with np.errstate(invalid="ignore", divide="ignore"):
        rows["precision"] = rows["k"] / rows["n_pred"]
        rows["yield"] = np.minimum(1.0, rows["n_pred"] / rows["n_real"])
    rows["fid"] = arc_metrics.direction_fidelity_yield(pred, real, target_gene=tg)
    rows["reach"] = arc_metrics.direction_reach(pred, real, target_gene=tg)
    rows["jaccard"] = arc_metrics.sig_jaccard(pred, real, target_gene=tg)
    rows["nmae"] = arc_metrics.lfc_nmae(pred, real, target_gene=tg)
    return rows


def atlasshift_model():
    """AtlasShift's ``model.py``, vendored verbatim (ablation arm only)."""
    upstream = str(ROOT / "third_party" / "atlasshift")
    if upstream not in sys.path:
        sys.path.insert(0, upstream)
    import model

    return model


def atlasshift_template(pool: np.ndarray, rng: np.random.Generator):
    """AtlasShift ``control_template`` logic on an in-memory pool: 400 cells, each the
    mean composition of 4 depth-sorted donors, depth = rounded donor mean."""
    cells, pool_k = frozen.CELLS_PER_PERT, 4
    lib = pool.sum(axis=1)
    comp = pool / lib[:, None]
    mean = comp.mean(axis=0)
    bulk = pool.sum(axis=0) / pool.sum()
    chosen = rng.choice(len(pool), cells * pool_k, replace=False)
    chosen = chosen[np.argsort(lib[chosen], kind="stable")]
    template = comp[chosen].reshape(cells, pool_k, -1).mean(axis=1)
    depths = np.rint(lib[chosen].reshape(cells, pool_k).mean(axis=1)).astype(np.int64)
    return template, depths, mean, bulk


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--context", required=True)
    parser.add_argument(
        "--generator",
        choices=["g1", "dualmoment"],
        default="g1",
        help="g1 = frozen V1 transport; dualmoment = AtlasShift's generator (ablation only)",
    )
    args = parser.parse_args()
    started = time.time()
    out = OUTDIR / (args.context if args.generator == "g1" else f"{args.context}__dualmoment")
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    delta, control_means, contexts, perts, shared_genes = frozen.load_design()
    target = contexts.index(args.context)
    sources = tuple(i for i in range(len(contexts)) if i != target)

    cell_counts = np.load(frozen.CANONICAL / "cell_counts.npy")[target]
    order = np.argsort(-cell_counts)[: frozen.N_PERTURBATIONS]
    order = order[np.argsort([perts[i] for i in order])]
    chosen = [perts[i] for i in order]

    path = frozen.RAW / f"{args.context}_processed_complete.h5ad"
    counts, group, var_names, diagnostics = frozen.load_target_cells(path, chosen)
    if diagnostics["max_round_trip_error"] > 1e-4:
        raise SystemExit("count recovery does not round-trip")
    control_rows = np.flatnonzero(group == -1)
    perm = rng.permutation(control_rows)
    n_ref = min(frozen.N_REFERENCE_CONTROLS, len(control_rows) // 3)
    reference_rows = np.sort(perm[:n_ref])
    pool_rows = np.sort(perm[n_ref:])
    ctrl_reference = np.asarray(counts[reference_rows].todense(), dtype=np.int64)
    ctrl_pool = np.asarray(counts[pool_rows].todense(), dtype=np.float64)
    real_blocks = [
        np.asarray(counts[np.flatnonzero(group == i)].todense(), dtype=np.int32)
        for i in range(len(chosen))
    ]
    del counts, group

    # ---- frozen V1 predictions on the source-space axis ---------------------
    m_hat = mr.MAIN_EFFECT_ESTIMATORS["M3b_basal_shrunk"](
        delta, sources, control_means=control_means, target=target
    )
    beta = {
        2: mr.tier_beta(delta, mr.source_subsets(sources, 2)[0], 2),
        1: mr.tier_beta(delta, mr.source_subsets(sources, 1)[0], 1),
        0: mr.tier_beta(delta, sources, 0),
    }
    pred_by_tier = {t: mr.predict(m_hat, beta[t], TIER_WEIGHTS[t])[order] for t in (2, 1, 0)}
    assignment = np.zeros(len(chosen), dtype=int)
    shuffled = np.random.default_rng(SEED + 1).permutation(len(chosen))
    assignment[shuffled[: ARC_PREVALENCE[2]]] = 2
    assignment[shuffled[ARC_PREVALENCE[2] : ARC_PREVALENCE[2] + ARC_PREVALENCE[1]]] = 1
    mixed = np.stack([pred_by_tier[t][i] for i, t in enumerate(assignment)])

    panel_map = build_panel_map(var_names, shared_genes)
    lib = ctrl_pool.sum(axis=1, keepdims=True)
    ctrl_log1p = np.log1p(
        1.0e4 * np.divide(ctrl_pool, lib, out=np.zeros_like(ctrl_pool), where=lib > 0)
    ).mean(axis=0)

    # ---- reference side ------------------------------------------------------
    ctrl_profile = arc_metrics.bulk_profile(ctrl_reference)
    ctrl_disp = arc_metrics.jackknife_dispersion(ctrl_reference)
    real_profiles, real_disp = frozen.bulk_block(real_blocks)
    tg = np.array([var_names.get_loc(p) if p in var_names else -1 for p in chosen])
    exclude = np.asarray(var_names.isin(set(perts)), dtype=bool)
    real_de = arc_metrics.de_table(real_blocks, ctrl_reference)
    keep = ~exclude
    true_norm = np.linalg.norm((real_profiles - ctrl_profile)[:, keep], axis=1)

    # oracle mean response of the held-out context, in CP10K, as one lfc vector
    cp10k = np.stack(
        [(b / np.maximum(b.sum(axis=1, keepdims=True), 1) * 1e4).mean(axis=0) for b in real_blocks]
    )
    pool_cp10k = (ctrl_pool / np.maximum(lib, 1) * 1e4).mean(axis=0)
    mr_lfc = np.log2((cp10k.mean(axis=0) + 1e-3) / (pool_cp10k + 1e-3))
    mr_lfc = np.where((cp10k.mean(axis=0) > 0) | (pool_cp10k > 0), mr_lfc, 0.0)

    def lfc_of(delta_rows: np.ndarray) -> np.ndarray:
        return bundle.panel_log2_fold_change(panel_map, ctrl_log1p, delta_rows)

    arms: dict[str, tuple] = {}
    for t in (2, 1, 0):
        arms[f"T{t}"] = (lfc_of(pred_by_tier[t]), np.full(len(chosen), t))
    for a in AMPLITUDES:
        arms[f"MIX_a{a}"] = (lfc_of(a * mixed), assignment)
    arms["ANCHOR_mean_response"] = (np.tile(mr_lfc, (len(chosen), 1)), assignment)

    score_rows, pert_frames = [], []

    def score(
        name: str,
        blocks,
        tier,
        *,
        real=real_de,
        ctrl=ctrl_profile,
        cdisp=ctrl_disp,
        rprof=real_profiles,
        rdisp=real_disp,
        ref_ctrl=ctrl_reference,
        amp=np.nan,
    ):
        profiles, disp = frozen.bulk_block(blocks)
        de = frozen.merge_de(
            [arc_metrics.de_table([b], ref_ctrl, tested=real.tested) for b in blocks],
            real.tested,
        )
        sc = frozen.score_generator(
            profiles,
            disp,
            de,
            real_de=real,
            real_profiles=rprof,
            real_disp=rdisp,
            ctrl_profile=ctrl,
            ctrl_disp=cdisp,
            exclude=exclude,
            target_gene=tg,
        )
        score_rows.append({"arm": name, "amplitude": amp, **sc})
        rows = per_pert_de(de, real, tg)
        frame = pd.DataFrame(rows)
        frame.insert(0, "perturbation", chosen)
        frame.insert(1, "arm", name)
        frame.insert(2, "tier", tier)
        frame["amplitude"] = amp
        frame["pred_norm"] = np.linalg.norm((profiles - ctrl)[:, keep], axis=1)
        frame["true_norm"] = true_norm
        pert_frames.append(frame)
        print(
            f"  {name:24s} FID={sc['de_wilcoxon_direction_fidelity_yield_raw']:.3f} "
            f"PDS={sc['pds_cosine']:.3f} ({time.time() - started:.0f}s)",
            flush=True,
        )

    if args.generator == "dualmoment":
        template, depths, comp_mean, comp_bulk = atlasshift_template(ctrl_pool, rng)

    def emit(lfc_row: np.ndarray, i: int) -> np.ndarray:
        if args.generator == "g1":
            return generate.transport_controls(
                ctrl_pool, lfc_row, frozen.CELLS_PER_PERT, rng=rng, smoothing=0.5
            ).astype(np.int32)
        # identical V1 effect, AtlasShift's emission: both moments scaled by 2^lfc
        mult = np.exp2(lfc_row)
        p_mean = comp_mean * mult
        p_bulk = comp_bulk * mult
        return (
            atlasshift_model()
            .dual_moment_counts(
                template,
                p_mean / p_mean.sum(),
                p_bulk / p_bulk.sum(),
                depths=depths,
                seed=SEED + i,
            )
            .astype(np.int32)
        )

    for name, (lfc, tier) in arms.items():
        blocks = [emit(lfc[i], i) for i in range(len(chosen))]
        amp = float(name.split("_a")[1]) if name.startswith("MIX_a") else 1.0
        score(name, blocks, tier, amp=amp)
        del blocks

    # ---- split-half replicate anchor (one split) -----------------------------
    half_rng = np.random.default_rng(SEED + 2)
    ref_perm = half_rng.permutation(len(ctrl_reference))
    ctrl_a = ctrl_reference[ref_perm[: len(ref_perm) // 2]]
    ctrl_b = ctrl_reference[ref_perm[len(ref_perm) // 2 :]]
    half_a, half_b = [], []
    for b in real_blocks:
        idx = half_rng.permutation(len(b))
        half_a.append(b[idx[: len(b) // 2]])
        half_b.append(b[idx[len(b) // 2 :]])
    real_b = arc_metrics.de_table(half_b, ctrl_b)
    prof_b, disp_b = frozen.bulk_block(half_b)
    score(
        "ANCHOR_split_half",
        half_a,
        assignment,
        real=real_b,
        ctrl=arc_metrics.bulk_profile(ctrl_b),
        cdisp=arc_metrics.jackknife_dispersion(ctrl_b),
        rprof=prof_b,
        rdisp=disp_b,
        ref_ctrl=ctrl_a,
    )

    pd.DataFrame(score_rows).to_csv(out / "scores_raw.csv", index=False)
    pd.concat(pert_frames).to_csv(out / "per_perturbation.csv", index=False)
    (out / "summary.json").write_text(
        json.dumps(
            {
                "context": args.context,
                "sources": [contexts[i] for i in sources],
                "n_perturbations": len(chosen),
                "n_reference_controls": int(n_ref),
                "n_pool_controls": int(len(pool_rows)),
                "tier_assignment_counts": {str(t): int((assignment == t).sum()) for t in (2, 1, 0)},
                "amplitudes": AMPLITUDES,
                "tier_weights": TIER_WEIGHTS,
                "elapsed_seconds": round(time.time() - started, 1),
            },
            indent=2,
        )
    )
    print(f"done {args.context} in {time.time() - started:.0f}s")


if __name__ == "__main__":
    main()
