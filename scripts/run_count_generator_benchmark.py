"""Count-space generator benchmark on public data (sections G, H, I).

The challenge scores raw counts, so the pseudobulk result of
``run_arc_count_space_baseline.py`` is only half a submission. This script
takes that frozen predictor, pushes it through three count generators, and
scores the generated cells with the local reimplementation of the six
``vcc2026`` metrics — on **public** data, against **real** held-out cells.

That is only possible because scPertEval's ``log1p(CP10K)`` matrices are
exactly invertible back to the integer counts they were built from
(:mod:`virtual_cell.data.counts`). The recovery is verified here before
anything is generated, and the run aborts if it does not round-trip.

The design mirrors Arc's as closely as public data allows:

* one held-out **context** (``replogle22k562``), its perturbation responses
  never read by the model, only by the scorer;
* its control cells split **disjointly** into a generator pool and a scoring
  reference, because Arc scores against held-out controls and never against
  the ones a submission was built from;
* the 300 best-measured shared perturbations, matching Arc's panel size so
  ``pds_cosine`` discriminates over a comparable field;
* exactly 400 generated cells per perturbation.

Reproduce: ``uv run python scripts/run_count_generator_benchmark.py``
"""

from __future__ import annotations

import json
import time
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from virtual_cell.analysis import loco
from virtual_cell.arc import bundle, generate
from virtual_cell.arc import metrics as arc_metrics
from virtual_cell.arc.panel import Support, build_panel_map
from virtual_cell.data import scperteval
from virtual_cell.data.counts import recover_counts, recovery_diagnostics
from virtual_cell.modelling import mean_response as mr

ROOT = Path(__file__).resolve().parents[1]
DESIGN = ROOT / "data" / "splits" / "four_context_v1"
CANONICAL = ROOT / "data" / "processed" / "four_context_v1"
RAW = ROOT / "data" / "raw" / "scperteval"
OUTDIR = ROOT / "outputs" / "arc_count_space_v1"

#: The public context held out for the count-space benchmark. Chosen before
#: any generator was run, on cell count alone: it has the most perturbations
#: with enough cells to make a single-cell ground truth worth scoring against.
TARGET_CONTEXT = "replogle22k562"

N_PERTURBATIONS = 300
CELLS_PER_PERT = 400
N_REFERENCE_CONTROLS = 2_000
SEED = 20260921

#: The operating point for the generator comparison: Tier 2 evidence (two
#: source contexts) at the frozen Tier-2 shrinkage. Chosen so the generators
#: are compared on the strongest prediction the frozen policy ever emits.
TIER = 2


def rule(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


# --------------------------------------------------------------------------
# data
# --------------------------------------------------------------------------


def load_design() -> tuple[np.ndarray, np.ndarray, list[str], list[str], list[str]]:
    delta = np.load(CANONICAL / "delta_tensor.npy").astype(np.float64)
    control_means = np.load(CANONICAL / "control_means.npy").astype(np.float64)
    contexts = (DESIGN / "contexts.txt").read_text().split()
    perts = (DESIGN / "shared_perturbations.txt").read_text().split()
    genes = (DESIGN / "shared_genes.txt").read_text().split()
    return delta, control_means, contexts, perts, genes


def load_target_cells(
    path: Path, wanted: Sequence[str]
) -> tuple[sparse.csr_matrix, np.ndarray, pd.Index, dict[str, float]]:
    """Recover integer counts for the control cells and the chosen perturbations.

    One streaming pass. ``group[i]`` indexes ``wanted`` for a perturbed cell and
    is ``-1`` for a control cell.
    """
    var_names = scperteval.read_var_names(path)
    labels = scperteval.cell_labels(path)
    row_of = {p: i for i, p in enumerate(wanted)}
    group_all = np.array(
        [row_of.get(lab, -2) if lab != scperteval.CONTROL_LABEL else -1 for lab in labels],
        dtype=np.int64,
    )

    blocks, groups = [], []
    diagnostics: dict[str, float] = {}
    for _start, chunk in scperteval.stream_row_chunks(path, chunk_size=20_000):
        offset = _start
        g_chunk = group_all[offset : offset + chunk.shape[0]]
        rows = np.flatnonzero(g_chunk >= -1)
        if not rows.size:
            continue
        sub = sparse.csr_matrix(chunk[rows])
        recovered = recover_counts(sub)
        if not diagnostics:
            diagnostics = recovery_diagnostics(sub, recovered)
        blocks.append(recovered.counts.astype(np.int32))
        groups.append(g_chunk[rows])
    return (
        sparse.vstack(blocks, format="csr"),
        np.concatenate(groups),
        var_names,
        diagnostics,
    )


# --------------------------------------------------------------------------
# structural fidelity (section H)
# --------------------------------------------------------------------------


def summarise_cells(
    lib: np.ndarray, detected: np.ndarray, spread: np.ndarray, density: float
) -> dict[str, float]:
    """The per-cell structure summary, from already-accumulated per-cell arrays.

    Shared by :func:`structure_stats` and the streaming generator loop so the
    two cannot report the same statistic differently.
    """
    return {
        "library_median": float(np.median(lib)),
        "library_q10": float(np.percentile(lib, 10)),
        "library_q90": float(np.percentile(lib, 90)),
        "library_max": float(lib.max(initial=0)),
        "genes_detected_median": float(np.median(detected)),
        "genes_detected_q10": float(np.percentile(detected, 10)),
        "genes_detected_q90": float(np.percentile(detected, 90)),
        "sparsity": float(1.0 - density),
        "heterogeneity_median": float(np.median(spread)),
    }


def structure_stats(blocks: Sequence[np.ndarray], label: str) -> dict[str, float]:
    """Per-cell structure, accumulated block by block so nothing is stacked.

    ``heterogeneity`` is each cell's distance to **its own group's** centroid in
    ``log1p(CP10K)`` space. Measuring it within a group is the point: a
    generator that emits 400 near-identical cells would still match a pooled
    spread, because the between-perturbation variation would cover for it.
    """
    lib_all, det_all, spread_all = [], [], []
    integer = non_negative = finite = True
    total, nonzero = 0, 0
    for block in blocks:
        x = np.asarray(block, dtype=np.float64)
        lib = x.sum(axis=1)
        integer = integer and bool(np.all(x == np.floor(x)))
        non_negative = non_negative and bool(np.all(x >= 0))
        finite = finite and bool(np.all(np.isfinite(x)))
        total += x.size
        nonzero += int((x > 0).sum())
        lib_all.append(lib)
        det_all.append((x > 0).sum(axis=1))
        logcpm = np.log1p(1.0e4 * x / np.where(lib > 0, lib, 1)[:, None])
        spread_all.append(np.linalg.norm(logcpm - logcpm.mean(axis=0, keepdims=True), axis=1))
    lib = np.concatenate(lib_all)
    return {
        "group": label,
        "n_cells": int(lib.size),
        "integer": integer,
        "non_negative": non_negative,
        "finite": finite,
        **summarise_cells(
            lib,
            np.concatenate(det_all),
            np.concatenate(spread_all),
            nonzero / total if total else np.nan,
        ),
    }


def realised_delta(counts: np.ndarray, control_profile: np.ndarray) -> np.ndarray:
    """The mean-``log1p(CP10K)`` response the generated cells actually carry."""
    x = np.asarray(counts, dtype=np.float64)
    lib = x.sum(axis=1, keepdims=True)
    logcpm = np.log1p(1.0e4 * np.divide(x, lib, out=np.zeros_like(x), where=lib > 0))
    return logcpm.mean(axis=0) - control_profile


# --------------------------------------------------------------------------
# scoring (section I)
# --------------------------------------------------------------------------


def bulk_block(blocks: Sequence[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    """Pseudobulk profile and jackknife dispersion for each group."""
    profiles = np.stack([arc_metrics.bulk_profile(b) for b in blocks])
    disp = np.array([arc_metrics.jackknife_dispersion(b) for b in blocks])
    return profiles, disp


def merge_de(rows: Sequence[arc_metrics.DETable], tested: np.ndarray) -> arc_metrics.DETable:
    """Stack single-perturbation DE tables into one, in the order given.

    Every metric that reads a ``DETable`` reads it row by row, and every row is
    computed from one perturbation's cells against the shared control group, so
    building the table a row at a time is identical to building it at once —
    and does not require all 300 blocks to be resident.
    """
    return arc_metrics.DETable(
        tested=tested,
        pval=np.vstack([r.pval for r in rows]),
        p_adj=np.vstack([r.p_adj for r in rows]),
        lfc=np.vstack([r.lfc for r in rows]),
    )


def score_generator(
    pred_profiles: np.ndarray,
    pred_disp: np.ndarray,
    pred_de: arc_metrics.DETable,
    *,
    real_de: arc_metrics.DETable,
    real_profiles: np.ndarray,
    real_disp: np.ndarray,
    ctrl_profile: np.ndarray,
    ctrl_disp: float,
    exclude: np.ndarray,
    target_gene: np.ndarray,
) -> dict[str, float]:
    """All six ``vcc2026`` members for one generator, in raw (unnormalised) units."""
    pds = arc_metrics.pds_cosine(
        pred_profiles - ctrl_profile[None, :],
        real_profiles - ctrl_profile[None, :],
        exclude=exclude,
    )
    mse = arc_metrics.expr_mse_unbiased_capped(
        pred_profiles,
        real_profiles,
        ctrl_profile,
        pred_dispersion=pred_disp,
        real_dispersion=real_disp,
        ctrl_dispersion=ctrl_disp,
        target_gene=target_gene,
    )
    return {
        "pds_cosine": float(np.mean(pds)),
        "expr_mse_unbiased_capped_norm": float(mse.value),
        "expr_mse_rho": float(mse.rho),
        "de_wilcoxon_direction_fidelity_yield_raw": float(
            np.nanmean(
                arc_metrics.direction_fidelity_yield(pred_de, real_de, target_gene=target_gene)
            )
        ),
        "de_wilcoxon_direction_reach_raw": float(
            np.nanmean(arc_metrics.direction_reach(pred_de, real_de, target_gene=target_gene))
        ),
        "de_wilcoxon_sig_jaccard": float(
            np.nanmean(arc_metrics.sig_jaccard(pred_de, real_de, target_gene=target_gene))
        ),
        "de_wilcoxon_lfc_nmae": float(
            np.nanmean(arc_metrics.lfc_nmae(pred_de, real_de, target_gene=target_gene))
        ),
    }


# --------------------------------------------------------------------------


def main() -> None:
    started = time.time()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    delta, control_means, contexts, perts, shared_genes = load_design()
    target_index = contexts.index(TARGET_CONTEXT)
    source_indices = tuple(i for i in range(len(contexts)) if i != target_index)

    rule("G0. SELECT THE PUBLIC COUNT-SPACE BENCHMARK")
    cell_counts = np.load(CANONICAL / "cell_counts.npy")[target_index]
    order = np.argsort(-cell_counts)[:N_PERTURBATIONS]
    order = order[np.argsort([perts[i] for i in order])]
    chosen = [perts[i] for i in order]
    source_names = [contexts[i] for i in source_indices]
    print(f"  held-out context: {TARGET_CONTEXT}  (sources {source_names})")
    print(
        f"  perturbations: {len(chosen)}  real cells/pert "
        f"min={cell_counts[order].min()} median={np.median(cell_counts[order]):.0f} "
        f"max={cell_counts[order].max()}"
    )

    rule("G1. RECOVER INTEGER COUNTS FROM THE PUBLIC log1p(CP10K) MATRIX")
    path = RAW / f"{TARGET_CONTEXT}_processed_complete.h5ad"
    counts, group, var_names, diagnostics = load_target_cells(path, chosen)
    print(f"  recovered {counts.shape[0]:,} cells x {counts.shape[1]:,} genes")
    for k, v in diagnostics.items():
        print(f"    {k:26s} {v:.6g}")
    if diagnostics["max_round_trip_error"] > 1e-4:
        raise SystemExit("count recovery does not round-trip; refusing to generate")
    pd.DataFrame([diagnostics]).to_csv(OUTDIR / "count_recovery_diagnostics.csv", index=False)

    control_rows = np.flatnonzero(group == -1)
    perm = rng.permutation(control_rows)
    reference_rows = np.sort(perm[:N_REFERENCE_CONTROLS])
    pool_rows = np.sort(perm[N_REFERENCE_CONTROLS:])
    ctrl_reference = np.asarray(counts[reference_rows].todense(), dtype=np.int64)
    ctrl_pool = np.asarray(counts[pool_rows].todense(), dtype=np.float64)
    print(
        f"  controls: {len(control_rows):,} total -> "
        f"{len(pool_rows):,} generator pool + {len(reference_rows):,} scoring reference "
        "(disjoint)"
    )

    real_blocks = [
        np.asarray(counts[np.flatnonzero(group == i)].todense(), dtype=np.int32)
        for i in range(len(chosen))
    ]
    del counts, group

    rule("G2. PREDICT THE PSEUDOBULK RESPONSE WITH THE FROZEN MODEL")
    m_hat = mr.MAIN_EFFECT_ESTIMATORS["M3b_basal_shrunk"](
        delta, source_indices, control_means=control_means, target=target_index
    )
    subset = mr.source_subsets(source_indices, mr.TIER_SOURCE_COUNT[TIER])[0]
    weights = mr.select_tier_shrinkage(
        delta,
        source_indices,
        control_means=control_means,
        main_effect=mr.MAIN_EFFECT_ESTIMATORS["M3b_basal_shrunk"],
        tiers=(2, 1),
    )
    beta = mr.tier_beta(delta, subset, TIER)
    predicted_shared = mr.predict(m_hat, beta, weights[TIER])[order]
    tier0_shared = mr.predict(m_hat, mr.tier_beta(delta, subset, 0), 0.0)[order]
    print(
        f"  m_hat=M3b_basal_shrunk  tier={TIER}  w={weights[TIER]}  "
        f"beta sources={[contexts[i] for i in subset]}"
    )

    rule("G3. MAP THE RESPONSE ONTO THE HELD-OUT CONTEXT'S OWN GENE AXIS")
    panel_map = build_panel_map(var_names, shared_genes)
    print("  " + panel_map.counts().to_string().replace("\n", "\n  "))

    lib = ctrl_pool.sum(axis=1, keepdims=True)
    ctrl_profile_log1p = np.log1p(
        1.0e4 * np.divide(ctrl_pool, lib, out=np.zeros_like(ctrl_pool), where=lib > 0)
    ).mean(axis=0)
    lfc = bundle.panel_log2_fold_change(panel_map, ctrl_profile_log1p, predicted_shared)
    lfc_tier0 = bundle.panel_log2_fold_change(panel_map, ctrl_profile_log1p, tier0_shared)

    def on_panel(delta_shared: np.ndarray) -> np.ndarray:
        """Lift a source-space delta onto the context axis, zero where unsupported."""
        take = panel_map.mask(Support.PREDICTED)
        out = np.zeros((delta_shared.shape[0], panel_map.n_panel), dtype=np.float64)
        out[:, take] = delta_shared[:, panel_map.source_index[take]]
        return out

    panel_delta = on_panel(predicted_shared)
    # Each arm is asked for a different response, so "did the generator realise
    # what it was asked for" has to be scored against that arm's own target.
    # G0 is asked for nothing, which is why its own-intent row is undefined.
    intended = {
        "G0_control_resample": np.zeros_like(panel_delta),
        "G1_transport": panel_delta,
        "G1_transport_unsmoothed": panel_delta,
        "G2_count_model": panel_delta,
        "G1_transport_tier0": on_panel(tier0_shared),
    }
    print(f"  |lfc| median {np.median(np.abs(lfc)):.4f}  max {np.abs(lfc).max():.3f}")

    rule("I0. REFERENCE SIDE: REAL CELLS, REAL CONTROLS")
    ctrl_profile = arc_metrics.bulk_profile(ctrl_reference)
    ctrl_disp = arc_metrics.jackknife_dispersion(ctrl_reference)
    real_profiles, real_disp = bulk_block(real_blocks)
    target_gene = np.array(
        [var_names.get_loc(p) if p in var_names else -1 for p in chosen], dtype=np.int64
    )
    exclude = np.asarray(var_names.isin(set(perts)), dtype=bool)
    print(f"  target genes on the axis: {(target_gene >= 0).sum()} / {len(chosen)}")
    t0 = time.time()
    real_de = arc_metrics.de_table(real_blocks, ctrl_reference)
    print(
        f"  reference DE table: {int(real_de.tested.sum())} genes tested ({time.time() - t0:.0f}s)"
    )
    real_structure = structure_stats(real_blocks, "real_perturbed")
    # Everything the reference side contributes is now summarised. Holding the
    # 300 real blocks past this point costs three gigabytes for nothing.
    del real_blocks

    rule("G4 / H / I. GENERATE, MEASURE STRUCTURE, SCORE — ONE GENERATOR AT A TIME")
    generators = {
        "G0_control_resample": lambda i: generate.resample_controls(
            ctrl_pool, CELLS_PER_PERT, rng=rng
        ),
        "G1_transport": lambda i: generate.transport_controls(
            ctrl_pool, lfc[i], CELLS_PER_PERT, rng=rng, smoothing=0.5
        ),
        "G1_transport_unsmoothed": lambda i: generate.transport_controls(
            ctrl_pool, lfc[i], CELLS_PER_PERT, rng=rng, smoothing=0.0
        ),
        "G2_count_model": lambda i: generate.count_model(
            ctrl_pool,
            bundle.predicted_profile(ctrl_profile_log1p, panel_delta[i]),
            CELLS_PER_PERT,
            rng=rng,
        ),
        "G1_transport_tier0": lambda i: generate.transport_controls(
            ctrl_pool, lfc_tier0[i], CELLS_PER_PERT, rng=rng, smoothing=0.5
        ),
    }

    # One perturbation's 400 cells exist at a time, and nothing else. Holding a
    # generator's full output is 400 x 300 x 8,563 dense counts -- four
    # gigabytes -- and every quantity below is computed per perturbation
    # anyway, so there is no reason to keep them and a machine-thrashing
    # reason not to.
    stats_rows = [
        structure_stats([ctrl_reference], "real_control_reference"),
        real_structure,
    ]
    fidelity_rows, score_rows = [], []
    g1_profiles: np.ndarray | None = None
    n_fidelity = min(50, len(chosen))

    for name, fn in generators.items():
        t_gen = t_score = 0.0
        profiles, disps, de_rows = [], [], []
        lib_all, det_all, spread_all, realised = [], [], [], []
        integer = non_negative = finite = True
        total = nonzero = 0

        for i in range(len(chosen)):
            t0 = time.time()
            block = fn(i).astype(np.int32)
            t_gen += time.time() - t0

            x = block.astype(np.float64)
            lib = x.sum(axis=1)
            integer = integer and bool(np.all(x == np.floor(x)))
            non_negative = non_negative and bool(np.all(x >= 0))
            finite = finite and bool(np.all(np.isfinite(x)))
            total += x.size
            nonzero += int((x > 0).sum())
            lib_all.append(lib)
            det_all.append((x > 0).sum(axis=1))
            logcpm = np.log1p(1.0e4 * x / np.where(lib > 0, lib, 1)[:, None])
            spread_all.append(np.linalg.norm(logcpm - logcpm.mean(axis=0, keepdims=True), axis=1))
            if i < n_fidelity:
                realised.append(logcpm.mean(axis=0) - ctrl_profile_log1p)

            t0 = time.time()
            profiles.append(arc_metrics.bulk_profile(block))
            disps.append(arc_metrics.jackknife_dispersion(block))
            de_rows.append(arc_metrics.de_table([block], ctrl_reference, tested=real_de.tested))
            t_score += time.time() - t0

        stats_rows.append(
            {
                "group": name,
                "n_cells": int(np.concatenate(lib_all).size),
                "integer": integer,
                "non_negative": non_negative,
                "finite": finite,
                **summarise_cells(
                    np.concatenate(lib_all),
                    np.concatenate(det_all),
                    np.concatenate(spread_all),
                    nonzero / total,
                ),
            }
        )

        got = np.stack(realised).ravel()
        row = {"generator": name, "realised_norm": float(np.linalg.norm(got))}
        # Two references, answering two different questions. "own" asks whether
        # the generator delivers the response it was handed -- a property of the
        # generator. "tier2" asks how much of the best available prediction each
        # arm ends up carrying -- a property of the arm, and the number that
        # makes G0 and the Tier-0 arm comparable with the rest.
        for label, target in (("own", intended[name]), ("tier2", panel_delta)):
            want = target[:n_fidelity].ravel()
            denom = float(want @ want)
            row[f"pearson_vs_{label}"] = (
                float(np.corrcoef(want, got)[0, 1]) if want.std() > 0 and got.std() > 0 else np.nan
            )
            row[f"slope_vs_{label}"] = float(want @ got / denom) if denom > 0 else np.nan
            row[f"{label}_intended_norm"] = float(np.linalg.norm(want))
        fidelity_rows.append(row)

        t0 = time.time()
        pred_profiles = np.stack(profiles)
        pred_de = merge_de(de_rows, real_de.tested)
        score_rows.append(
            {
                "generator": name,
                **score_generator(
                    pred_profiles,
                    np.array(disps),
                    pred_de,
                    real_de=real_de,
                    real_profiles=real_profiles,
                    real_disp=real_disp,
                    ctrl_profile=ctrl_profile,
                    ctrl_disp=ctrl_disp,
                    exclude=exclude,
                    target_gene=target_gene,
                ),
            }
        )
        t_score += time.time() - t0
        if name == "G1_transport":
            g1_profiles = pred_profiles
        print(f"  {name:26s} generated {t_gen:6.1f}s   scored {t_score:6.1f}s")

    rule("H. STRUCTURAL FIDELITY")
    structure = pd.DataFrame(stats_rows)
    structure.to_csv(OUTDIR / "generator_structure.csv", index=False)
    print(
        structure[
            [
                "group",
                "n_cells",
                "integer",
                "library_median",
                "genes_detected_median",
                "sparsity",
                "heterogeneity_median",
            ]
        ].to_string(index=False)
    )
    control_detected = structure.loc[0, "genes_detected_median"]
    for name in generators:
        row = structure[structure["group"] == name].iloc[0]
        change = row["genes_detected_median"] / control_detected - 1.0
        print(f"  {name:26s} genes detected vs control reference: {change:+.1%}")

    rule("H2. DOES THE GENERATOR REALISE THE INTENDED MEAN RESPONSE?")
    fidelity = pd.DataFrame(fidelity_rows)
    fidelity.to_csv(OUTDIR / "generator_mean_response_fidelity.csv", index=False)
    print(fidelity.to_string(index=False))

    rule("I. LOCAL cell-eval2 SCORING AGAINST REAL HELD-OUT CELLS")
    scores = pd.DataFrame(score_rows)
    scores.to_csv(OUTDIR / "generator_scores_raw.csv", index=False)
    print(scores.to_string(index=False))

    rule("I2. RELATIVE TO CONTROL RESAMPLING")
    floor = scores.set_index("generator").loc["G0_control_resample"]
    rel = scores.set_index("generator").copy()
    for member, higher in arc_metrics.SCORED_METRICS.items():
        sign = 1.0 if higher else -1.0
        rel[member] = sign * (rel[member] - floor[member])
    rel = rel[list(arc_metrics.SCORED_METRICS)]
    rel.to_csv(OUTDIR / "generator_scores_vs_control.csv")
    print("  improvement over G0 (positive = better than emitting controls)")
    print(rel.to_string())

    summary = {
        "target_context": TARGET_CONTEXT,
        "source_contexts": [contexts[i] for i in source_indices],
        "n_perturbations": len(chosen),
        "cells_per_perturbation": CELLS_PER_PERT,
        "tier": TIER,
        "shrinkage": {str(k): v for k, v in weights.weights.items()},
        "beta_sources": [contexts[i] for i in subset],
        "count_recovery": diagnostics,
        "panel_support": {k: int(v) for k, v in panel_map.counts().items()},
        "scores": scores.set_index("generator").to_dict(orient="index"),
        "elapsed_seconds": round(time.time() - started, 1),
    }
    (OUTDIR / "generator_summary.json").write_text(json.dumps(summary, indent=2, default=float))
    np.save(OUTDIR / "predicted_shared_delta.npy", predicted_shared.astype(np.float32))
    (OUTDIR / "benchmark_perturbations.txt").write_text("\n".join(chosen) + "\n")
    print(f"\nwritten to {OUTDIR}  ({summary['elapsed_seconds']}s)")

    # Reported for the record: per-perturbation pseudobulk agreement of the
    # prediction that drove the generators, on the held-out context's own axis.
    if g1_profiles is not None:
        per_pert = loco.per_perturbation_metrics(
            real_profiles - ctrl_profile[None, :],
            g1_profiles - ctrl_profile[None, :],
        )
        per_pert.insert(0, "perturbation", chosen)
        per_pert.to_csv(OUTDIR / "g1_per_perturbation.csv", index=False)


if __name__ == "__main__":
    main()
