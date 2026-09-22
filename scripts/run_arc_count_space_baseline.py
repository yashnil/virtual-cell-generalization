"""Arc count-space baseline — the pseudobulk half (sections A, B, D, E, F).

Freezes the research phase, re-derives the Arc target support policy from
frozen provenance, then answers the two questions the mean-response model has
free parameters for:

  D. which *feasible* context-wide perturbation baseline ``m_hat`` wins on
     public leave-one-context-out folds, and
  E. what shrinkage ``w`` each support tier deserves, chosen by nested
     validation inside the source contexts only,

and finally (F) measures the complete predictor ``m_hat + w beta_hat`` per
tier on the public held-out contexts.

No Arc perturbation outcome exists to leak; the target *context* of every fold
here is a public one, and its responses enter only at evaluation time.

Reproduce: ``uv run python scripts/run_arc_count_space_baseline.py``
"""

from __future__ import annotations

import json
import subprocess
import time
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from virtual_cell.analysis import loco
from virtual_cell.arc import metrics as arc_metrics
from virtual_cell.data import arc2026
from virtual_cell.modelling import context_main_effect as cme
from virtual_cell.modelling import mean_response as mr

ROOT = Path(__file__).resolve().parents[1]
DESIGN = ROOT / "data" / "splits" / "four_context_v1"
CANONICAL = ROOT / "data" / "processed" / "four_context_v1"
CONTROLS = ROOT / "data" / "raw" / "arc2026" / "controls"
PROVENANCE = ROOT / "data" / "provenance" / "scperteval"
INVENTORY = PROVENANCE / "public_label_inventory.json"
SPLITS = ROOT / "data" / "splits"
OUTDIR = ROOT / "outputs" / "arc_count_space_v1"

TIERS = (2, 1, 0)


def rule(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


# --------------------------------------------------------------------------
# A. freeze verification
# --------------------------------------------------------------------------


def verify_freezes() -> pd.DataFrame:
    """``shasum -c`` every frozen manifest. Any failure is fatal to the phase."""
    rows = []
    for manifest in sorted(PROVENANCE.glob("*freeze*.txt")):
        proc = subprocess.run(
            ["shasum", "-a", "256", "-c", str(manifest)],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        lines = [ln for ln in proc.stdout.splitlines() if ln.strip()]
        ok = sum(1 for ln in lines if ln.endswith(": OK"))
        bad = [ln for ln in lines if not ln.endswith(": OK")]
        rows.append(
            {
                "manifest": manifest.name,
                "files": len(lines),
                "ok": ok,
                "failed": len(bad),
                "detail": "; ".join(bad[:3]),
            }
        )
    frame = pd.DataFrame(rows)
    if frame["failed"].sum():
        raise SystemExit(f"FROZEN ARTIFACTS CHANGED:\n{frame[frame['failed'] > 0]}")
    return frame


# --------------------------------------------------------------------------
# B. target support policy
# --------------------------------------------------------------------------


def recompute_support(inventory: dict[str, dict]) -> pd.DataFrame:
    """Re-derive the frozen tiers from provenance, adding the per-target record.

    Identifier presence only. Nothing here reads a response value, so no model
    outcome can move a target between tiers — which is the standing policy.
    """
    panel = pd.Index(arc2026.load_gene_names(CONTROLS)).astype(object)
    targets = list(arc2026.load_pert_counts(CONTROLS)["target_gene"])
    panel_set = set(panel)

    gene_sets = {k: set(v["genes"]) for k, v in inventory.items()}
    pert_sets = {
        k: set(v["perturbations"]) - {"control", "non-targeting"} for k, v in inventory.items()
    }

    rows = []
    for target in targets:
        perturbed = sorted(n for n in inventory if target in pert_sets[n])
        measured = sorted(n for n in inventory if target in gene_sets[n])
        n = len(perturbed)
        # Output-gene availability: panel genes measured anywhere public.
        out_union = set().union(*(gene_sets[n_] for n_ in inventory)) & panel_set
        # Response space: panel genes on which a beta for THIS target is
        # estimable, i.e. present in every dataset that perturbed it.
        if perturbed:
            response_space = set.intersection(*(gene_sets[n_] for n_ in perturbed)) & panel_set
        else:
            response_space = set()
        rows.append(
            {
                "arc_target": target,
                "support_tier": 2 if n >= 2 else (1 if n == 1 else 0),
                "n_direct_source_contexts": n,
                "direct_source_contexts": "|".join(perturbed),
                "direct_source_cells": "|".join(
                    f"{n_}={inventory[n_]['cells']}" for n_ in perturbed
                ),
                "n_datasets_measuring_target_gene": len(measured),
                "panel_genes_output_available": len(out_union),
                "panel_genes_response_available": len(response_space),
                "target_gene_on_panel": target in panel_set,
            }
        )
    return pd.DataFrame(rows)


def compare_to_frozen(fresh: pd.DataFrame) -> dict[str, object]:
    frozen = pd.read_csv(SPLITS / "arc_target_support_v1.csv")
    merged = frozen[["arc_target", "support_tier", "n_contexts_perturbed"]].merge(
        fresh[["arc_target", "support_tier", "n_direct_source_contexts"]],
        on="arc_target",
        suffixes=("_frozen", "_fresh"),
    )
    return {
        "n_targets": int(len(merged)),
        "tier_mismatches": int(
            (merged["support_tier_frozen"] != merged["support_tier_fresh"]).sum()
        ),
        "count_mismatches": int(
            (merged["n_contexts_perturbed"] != merged["n_direct_source_contexts"]).sum()
        ),
    }


# --------------------------------------------------------------------------
# shared evaluation helpers
# --------------------------------------------------------------------------


def load_canonical() -> tuple[np.ndarray, np.ndarray, list[str], list[str], list[str]]:
    delta = np.load(CANONICAL / "delta_tensor.npy").astype(np.float64)
    control_means = np.load(CANONICAL / "control_means.npy").astype(np.float64)
    contexts = (DESIGN / "contexts.txt").read_text().split()
    perts = (DESIGN / "shared_perturbations.txt").read_text().split()
    genes = (DESIGN / "shared_genes.txt").read_text().split()
    return delta, control_means, contexts, perts, genes


def energy_explained(truth: np.ndarray, pred: np.ndarray) -> float:
    err = truth - pred
    denom = float(np.sum(truth * truth))
    return float(1.0 - np.sum(err * err) / denom) if denom > 0 else np.nan


def vector_metrics(truth: np.ndarray, pred: np.ndarray) -> dict[str, float]:
    """Metrics for a single gene-space vector (used for ``m_hat`` itself)."""
    t = np.asarray(truth, dtype=np.float64)
    p = np.asarray(pred, dtype=np.float64)
    err = t - p
    norm = np.linalg.norm(t) * np.linalg.norm(p)
    return {
        "pearson": float(np.corrcoef(t, p)[0, 1]) if p.std() > 0 else np.nan,
        "spearman": float(stats.spearmanr(t, p).statistic) if p.std() > 0 else np.nan,
        "cosine": float(t @ p / norm) if norm > 0 else np.nan,
        "mse": float(np.mean(err * err)),
        "energy_explained": energy_explained(t, p),
    }


def matrix_metrics(
    truth: np.ndarray, pred: np.ndarray, *, exclude: np.ndarray, with_spearman: bool = True
) -> dict[str, float]:
    """Per-perturbation metrics plus the two aggregate ones, for a (P, G) block."""
    per = loco.per_perturbation_metrics(truth, pred, with_spearman=with_spearman)
    pds = arc_metrics.pds_cosine(pred, truth, exclude=exclude)
    return {
        "pearson_median": float(np.nanmedian(per["pearson"])),
        # nanmedian of an all-NaN column warns; with_spearman=False makes it so.
        "spearman_median": float(np.nanmedian(per["spearman"])) if with_spearman else np.nan,
        "cosine_median": float(np.nanmedian(per["cosine"])),
        "mse_mean": float(np.mean(per["mse"])),
        "energy_explained_pooled": energy_explained(truth, pred),
        "energy_explained_median": float(np.nanmedian(per["energy_explained"])),
        "pds_cosine_mean": float(np.mean(pds)),
    }


# --------------------------------------------------------------------------
# D. the feasible context-wide baseline
# --------------------------------------------------------------------------


def main_effect_benchmark(
    delta: np.ndarray, control_means: np.ndarray, contexts: Sequence[str], exclude: np.ndarray
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Score every feasible ``m_hat`` against the hidden target main effect.

    Two views, because they answer different questions. The *main-effect* view
    asks how close ``m_hat`` is to ``m_c`` — the quantity the challenge score
    anchors zero at. The *Tier-0 prediction* view asks what ``m_hat`` alone
    scores as a complete response prediction, which is exactly what a Tier-0
    Arc target receives.
    """
    folds = loco.make_folds(list(contexts))
    main_rows, pred_rows = [], []
    for fold in folds:
        truth_main = cme.oracle_main_effect(delta, fold.target_index)
        truth_full = delta[fold.target_index]
        for name, fn in mr.MAIN_EFFECT_ESTIMATORS.items():
            m_hat = fn(
                delta,
                fold.source_indices,
                control_means=control_means,
                target=fold.target_index,
            )
            main_rows.append(
                {"fold": fold.target, "estimator": name, **vector_metrics(truth_main, m_hat)}
            )
            pred = np.broadcast_to(m_hat, truth_full.shape)
            pred_rows.append(
                {
                    "fold": fold.target,
                    "estimator": name,
                    **matrix_metrics(truth_full, pred, exclude=exclude, with_spearman=False),
                }
            )
    return pd.DataFrame(main_rows), pd.DataFrame(pred_rows)


# --------------------------------------------------------------------------
# E. tier shrinkage
# --------------------------------------------------------------------------


def shrinkage_selection(
    delta: np.ndarray,
    control_means: np.ndarray,
    contexts: Sequence[str],
    estimator_name: str,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, dict[int, float]]]:
    """Nested (inner) selection, and the outer oracle sweep for comparison only."""
    fn = mr.MAIN_EFFECT_ESTIMATORS[estimator_name]
    folds = loco.make_folds(list(contexts))
    inner_rows, outer_rows = [], []
    chosen: dict[str, dict[int, float]] = {}

    for fold in folds:
        selected = mr.select_tier_shrinkage(
            delta,
            fold.source_indices,
            control_means=control_means,
            main_effect=fn,
            tiers=(2, 1),
        )
        chosen[fold.target] = selected.weights
        for (tier, weight), value in sorted(selected.grid.items()):
            inner_rows.append(
                {"fold": fold.target, "tier": tier, "weight": weight, "inner_mse": value}
            )

        truth = delta[fold.target_index]
        m_hat = fn(
            delta, fold.source_indices, control_means=control_means, target=fold.target_index
        )
        for tier in (2, 1):
            size = mr.TIER_SOURCE_COUNT[tier]
            subsets = mr.source_subsets(fold.source_indices, size)
            for weight in mr.SHRINKAGE_GRID:
                errs = []
                for subset in subsets:
                    pred = mr.predict(m_hat, mr.centred_source_beta(delta, subset), weight)
                    diff = truth - pred
                    errs.append(float(np.mean(diff * diff)))
                outer_rows.append(
                    {
                        "fold": fold.target,
                        "tier": tier,
                        "weight": weight,
                        "outer_mse": float(np.mean(errs)),
                    }
                )
    return pd.DataFrame(inner_rows), pd.DataFrame(outer_rows), chosen


# --------------------------------------------------------------------------
# F. the complete pseudobulk predictor
# --------------------------------------------------------------------------


def complete_benchmark(
    delta: np.ndarray,
    control_means: np.ndarray,
    contexts: Sequence[str],
    exclude: np.ndarray,
    estimator_name: str,
    weights_by_fold: dict[str, dict[int, float]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """``m_hat + w beta_hat`` per tier, on every public held-out context."""
    fn = mr.MAIN_EFFECT_ESTIMATORS[estimator_name]
    folds = loco.make_folds(list(contexts))
    rows, per_pert = [], []
    for fold in folds:
        truth = delta[fold.target_index]
        m_hat = fn(
            delta, fold.source_indices, control_means=control_means, target=fold.target_index
        )
        weights = weights_by_fold[fold.target]
        for tier in TIERS:
            size = mr.TIER_SOURCE_COUNT[tier]
            weight = weights.get(tier, 0.0)
            subsets = (
                [tuple(fold.source_indices)]
                if size is None
                else mr.source_subsets(fold.source_indices, size)
            )
            block, pert_blocks = [], []
            for subset in subsets:
                beta = mr.tier_beta(delta, subset, tier)
                pred = mr.predict(m_hat, beta, weight)
                block.append(matrix_metrics(truth, pred, exclude=exclude))
                pert_blocks.append(loco.per_perturbation_metrics(truth, pred))
            summary = {k: float(np.mean([b[k] for b in block])) for k in block[0]}
            rows.append(
                {
                    "fold": fold.target,
                    "tier": tier,
                    "weight": weight,
                    "n_source_subsets": len(subsets),
                    **summary,
                }
            )
            stacked = sum(b[["pearson", "cosine", "mse", "energy_explained"]] for b in pert_blocks)
            stacked = stacked / len(pert_blocks)
            stacked["fold"] = fold.target
            stacked["tier"] = tier
            per_pert.append(stacked)
    return pd.DataFrame(rows), pd.concat(per_pert, ignore_index=True)


def mixed_panel_benchmark(
    delta: np.ndarray,
    control_means: np.ndarray,
    contexts: Sequence[str],
    exclude: np.ndarray,
    estimator_name: str,
    weights_by_fold: dict[str, dict[int, float]],
    support: pd.DataFrame,
    *,
    n_draws: int = 20,
    seed: int = 20260921,
) -> pd.DataFrame:
    """One prediction matrix mixing all three tiers at Arc's own prevalence.

    The per-tier table answers "how well does Tier 2 do", but an Arc submission
    is not a Tier-2 submission: it is 7 Tier-2 targets, 79 Tier-1 and 214
    Tier-0, scored together. That distinction is invisible for a per-
    perturbation metric — averaging the per-tier medians would be close — but
    ``pds_cosine`` is **not** per-perturbation. It ranks each prediction against
    every *other* prediction, so its value depends on the composition of the
    field: 214 predictions that are all identical (Tier 0 emits the same
    ``m_hat`` for every target) sit at distance zero from one another and drag
    the discrimination down in a way no per-tier number shows.

    So the mixture is built explicitly: each public perturbation is assigned a
    tier by a seeded draw at Arc's prevalence, one prediction matrix is
    assembled, and the metrics are computed once over it. Repeated over
    ``n_draws`` assignments, because which perturbation lands in which tier is
    arbitrary and the spread over draws is part of the answer.
    """
    fn = mr.MAIN_EFFECT_ESTIMATORS[estimator_name]
    folds = loco.make_folds(list(contexts))
    share = support["support_tier"].value_counts(normalize=True)
    probabilities = np.array([share.get(t, 0.0) for t in TIERS])
    n_pert = delta.shape[1]
    rows = []
    for fold in folds:
        truth = delta[fold.target_index]
        m_hat = fn(
            delta, fold.source_indices, control_means=control_means, target=fold.target_index
        )
        weights = weights_by_fold[fold.target]
        # One beta per tier, using that tier's own source-count structure.
        beta_by_tier = {
            tier: mr.tier_beta(
                delta,
                (
                    tuple(fold.source_indices)
                    if mr.TIER_SOURCE_COUNT[tier] is None
                    else mr.source_subsets(fold.source_indices, mr.TIER_SOURCE_COUNT[tier])[0]
                ),
                tier,
            )
            for tier in TIERS
        }
        rng = np.random.default_rng(seed)
        for draw in range(n_draws):
            assignment = rng.choice(TIERS, size=n_pert, p=probabilities)
            pred = np.empty_like(truth)
            for tier in TIERS:
                rows_of_tier = assignment == tier
                if rows_of_tier.any():
                    pred[rows_of_tier] = mr.predict(
                        m_hat, beta_by_tier[tier][rows_of_tier], weights.get(tier, 0.0)
                    )
            rows.append(
                {
                    "fold": fold.target,
                    "draw": draw,
                    "n_tier2": int((assignment == 2).sum()),
                    "n_tier1": int((assignment == 1).sum()),
                    "n_tier0": int((assignment == 0).sum()),
                    **matrix_metrics(truth, pred, exclude=exclude, with_spearman=False),
                }
            )
    return pd.DataFrame(rows)


def arc_weighted_summary(frame: pd.DataFrame, support: pd.DataFrame) -> pd.DataFrame:
    """Re-weight the per-tier public result by the Arc panel's own tier mix.

    An approximation for the per-perturbation metrics and a poor one for
    ``pds_cosine``; :func:`mixed_panel_benchmark` is the honest version.
    """
    share = support["support_tier"].value_counts(normalize=True)
    numeric = [c for c in frame.columns if frame[c].dtype.kind == "f" and c != "weight"]
    out = []
    for fold, block in frame.groupby("fold"):
        indexed = block.set_index("tier")
        row = {"fold": fold}
        for col in numeric:
            row[col] = float(sum(share.get(t, 0.0) * indexed.loc[t, col] for t in TIERS))
        out.append(row)
    return pd.DataFrame(out)


# --------------------------------------------------------------------------


def main() -> None:
    started = time.time()
    OUTDIR.mkdir(parents=True, exist_ok=True)

    rule("A. FREEZE VERIFICATION (before)")
    freezes = verify_freezes()
    print(freezes.to_string(index=False))
    freezes.to_csv(OUTDIR / "freeze_verification_before.csv", index=False)

    rule("B. ARC TARGET SUPPORT POLICY, RE-DERIVED FROM FROZEN PROVENANCE")
    with INVENTORY.open() as fh:
        inventory = json.load(fh)
    support = recompute_support(inventory)
    support.to_csv(OUTDIR / "arc_target_support_recomputed.csv", index=False)
    check = compare_to_frozen(support)
    print(f"  targets: {check['n_targets']}")
    print(f"  tier mismatches vs frozen split: {check['tier_mismatches']}")
    print(f"  context-count mismatches:        {check['count_mismatches']}")
    if check["tier_mismatches"] or check["count_mismatches"]:
        raise SystemExit("support policy is not reproducible from provenance")
    for tier in TIERS:
        block = support[support["support_tier"] == tier]
        print(
            f"  TIER {tier}: {len(block):>3} targets | "
            f"median response-space genes {block['panel_genes_response_available'].median():.0f} | "
            f"on-panel targets {int(block['target_gene_on_panel'].sum())}"
        )

    delta, control_means, contexts, perts, genes = load_canonical()
    gene_index = pd.Index(genes)
    exclude = np.asarray(gene_index.isin(set(perts)), dtype=bool)
    print(f"\n  public design: {delta.shape} contexts={contexts}")
    print(f"  perturbation-target genes on the public axis: {int(exclude.sum())} / {len(genes)}")

    rule("D. FEASIBLE CONTEXT-WIDE BASELINE m_hat")
    main_frame, tier0_frame = main_effect_benchmark(delta, control_means, contexts, exclude)
    main_frame.to_csv(OUTDIR / "main_effect_vs_oracle.csv", index=False)
    tier0_frame.to_csv(OUTDIR / "main_effect_as_prediction.csv", index=False)
    agg = main_frame.groupby("estimator")[["pearson", "cosine", "mse", "energy_explained"]].mean()
    print("\n  m_hat vs the hidden target main effect m_c (mean over 4 folds)")
    print(agg.to_string())
    agg0 = tier0_frame.groupby("estimator")[
        [
            "pearson_median",
            "cosine_median",
            "mse_mean",
            "energy_explained_pooled",
            "pds_cosine_mean",
        ]
    ].mean()
    print("\n  m_hat alone as a complete response prediction (the Tier-0 case)")
    print(agg0.to_string())

    winner = str(agg0["mse_mean"].idxmin())
    per_fold_best = tier0_frame.loc[tier0_frame.groupby("fold")["mse_mean"].idxmin()]
    print(f"\n  WINNER (lowest mean MSE as a prediction): {winner}")
    per_fold_winner = dict(zip(per_fold_best["fold"], per_fold_best["estimator"], strict=True))
    print(f"  per-fold winners: {per_fold_winner}")
    (OUTDIR / "main_effect_choice.json").write_text(
        json.dumps(
            {
                "winner": winner,
                "per_fold_winner": per_fold_winner,
                "criterion": "mean squared error of m_hat as a complete response prediction",
            },
            indent=2,
        )
    )

    rule("E. TIER-SPECIFIC SHRINKAGE (nested, sources only)")
    inner, outer, chosen = shrinkage_selection(delta, control_means, contexts, winner)
    inner.to_csv(OUTDIR / "shrinkage_inner_grid.csv", index=False)
    outer.to_csv(OUTDIR / "shrinkage_outer_grid.csv", index=False)
    print("\n  inner (selection) MSE, mean over folds")
    print(inner.groupby(["tier", "weight"])["inner_mse"].mean().unstack().to_string())
    print("\n  outer (oracle, reported only) MSE, mean over folds")
    print(outer.groupby(["tier", "weight"])["outer_mse"].mean().unstack().to_string())
    for fold, weights in chosen.items():
        print(f"  selected {fold}: tier2 w={weights[2]}  tier1 w={weights[1]}  tier0 w=0.0")
    outer_best = outer.groupby(["tier", "weight"])["outer_mse"].mean().groupby("tier").idxmin()
    print(f"  oracle-optimal (NOT used): {dict(outer_best)}")
    (OUTDIR / "shrinkage_choice.json").write_text(
        json.dumps(
            {"per_fold": {k: {str(t): w for t, w in v.items()} for k, v in chosen.items()}},
            indent=2,
        )
    )

    rule("F. COMPLETE PSEUDOBULK PREDICTOR, BY TIER")
    complete, per_pert = complete_benchmark(delta, control_means, contexts, exclude, winner, chosen)
    complete.to_csv(OUTDIR / "complete_pseudobulk_by_tier.csv", index=False)
    per_pert.to_csv(OUTDIR / "complete_pseudobulk_per_perturbation.csv", index=False)
    cols = [
        "pearson_median",
        "spearman_median",
        "cosine_median",
        "mse_mean",
        "energy_explained_pooled",
        "pds_cosine_mean",
    ]
    print(complete.set_index(["fold", "tier"])[["weight", *cols]].to_string())
    print("\n  mean over folds")
    print(complete.groupby("tier")[cols].mean().to_string())

    weighted = arc_weighted_summary(complete, support)
    weighted.to_csv(OUTDIR / "arc_tier_weighted_summary.csv", index=False)
    print("\n  Arc-prevalence-weighted per-tier averages (7 / 79 / 214 targets)")
    print(weighted.to_string(index=False))

    rule("F2. ONE MIXED PANEL AT ARC'S TIER PREVALENCE")
    mixed = mixed_panel_benchmark(delta, control_means, contexts, exclude, winner, chosen, support)
    mixed.to_csv(OUTDIR / "mixed_panel_draws.csv", index=False)
    mixed_cols = [c for c in cols if c != "spearman_median"]
    summary_mixed = mixed.groupby("fold")[mixed_cols].agg(["mean", "std"])
    print(summary_mixed.round(4).to_string())
    print("\n  mean over folds")
    print(mixed[mixed_cols].mean().round(4).to_string())
    print(
        "  pds_cosine over a mixed panel vs the Tier-2-only field: "
        f"{mixed['pds_cosine_mean'].mean():.4f} vs "
        f"{complete[complete['tier'] == 2]['pds_cosine_mean'].mean():.4f}"
    )

    rule("A. FREEZE VERIFICATION (after)")
    after = verify_freezes()
    after.to_csv(OUTDIR / "freeze_verification_after.csv", index=False)
    print(f"  {int(after['ok'].sum())} files verified across {len(after)} manifests, 0 failures")

    summary = {
        "main_effect_winner": winner,
        "shrinkage_per_fold": {k: {str(t): w for t, w in v.items()} for k, v in chosen.items()},
        "tier_counts": {str(t): int((support["support_tier"] == t).sum()) for t in TIERS},
        "support_check": check,
        "complete_by_tier": complete.groupby("tier")[cols].mean().to_dict(),
        "mixed_panel": mixed[mixed_cols].mean().to_dict(),
        "elapsed_seconds": round(time.time() - started, 1),
    }
    (OUTDIR / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\nwritten to {OUTDIR}  ({summary['elapsed_seconds']}s)")


if __name__ == "__main__":
    main()
