"""Unseen-perturbation generalization (v1).

The four-context programme answered: *given a perturbation measured elsewhere,
can its response be transferred to a context never seen?* Arc's validation
panel mostly asks something else -- 214 of its 300 targets are perturbed
nowhere in public data -- so the question becomes: *can a perturbation's
conserved effect be predicted from what is known about the gene, with no
measurement of that knockdown anywhere?*

This script builds that benchmark and runs it. Two axes are held out
independently and, in the hardest regime, together:

    P1   held-out perturbation, seen context
    P2   held-out perturbation AND held-out context   <- the Arc-like regime

Every estimator is low capacity by design: zero, nearest biological neighbour,
k-neighbour transfer, ridge. No MLP, no deep model, no gamma model.

Usage::

    uv run python scripts/run_unseen_perturbation.py
"""

from __future__ import annotations

import json
import time
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from virtual_cell.analysis import foundations
from virtual_cell.data import arc2026
from virtual_cell.modelling import context_main_effect as cme
from virtual_cell.modelling import unseen_perturbation as up
from virtual_cell.priors import catalogue
from virtual_cell.priors import features as pf

ROOT = Path(__file__).resolve().parents[1]
DESIGN = ROOT / "data" / "splits" / "four_context_v1"
CANONICAL = ROOT / "data" / "processed" / "four_context_v1"
CONTROLS = ROOT / "data" / "raw" / "arc2026" / "controls"
SPLITS = ROOT / "data" / "splits"
OUTDIR = ROOT / "outputs" / "unseen_perturbation_v1"
INVENTORY = ROOT / "data" / "provenance" / "scperteval" / "public_label_inventory.json"

N_FOLDS = 5
SEED = 20260920
KNN_K = 25
RIDGE_ALPHA = 10.0
PCA_COMPONENTS = 64


def rule(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


# --------------------------------------------------------------------------
# inputs
# --------------------------------------------------------------------------


def load_canonical() -> tuple[np.ndarray, np.ndarray, list[str], list[str], list[str]]:
    contexts = (DESIGN / "contexts.txt").read_text().split()
    perts = (DESIGN / "shared_perturbations.txt").read_text().split()
    genes = (DESIGN / "shared_genes.txt").read_text().split()
    delta = np.load(CANONICAL / "delta_tensor.npy").astype(np.float64)
    control = np.load(CANONICAL / "control_means.npy").astype(np.float64)
    return delta, control, contexts, perts, genes


# --------------------------------------------------------------------------
# D. support tiers for the 300 Arc targets
# --------------------------------------------------------------------------


def support_tiers(inventory: dict[str, dict]) -> pd.DataFrame:
    """Classify each Arc target by how many public contexts perturb it.

    Identifier presence only: a gene is 'perturbed in' a dataset when its
    symbol appears in that dataset's perturbation labels. No response data is
    read, so this classification cannot be influenced by any outcome.
    """
    targets = list(arc2026.load_pert_counts(CONTROLS)["target_gene"])
    rows = []
    for target in targets:
        observed, measured, cells = [], [], []
        for name, rec in sorted(inventory.items()):
            perturbed = target in set(rec["perturbations"]) - {"control", "non-targeting"}
            if perturbed:
                observed.append(name)
            if target in set(rec["genes"]):
                measured.append(name)
            cells.append(f"{name}={rec['cells']}")
        n = len(observed)
        rows.append(
            {
                "arc_target": target,
                "support_tier": 2 if n >= 2 else (1 if n == 1 else 0),
                "n_contexts_perturbed": n,
                "datasets_perturbed": "|".join(observed),
                "n_datasets_measured": len(measured),
                "datasets_measured": "|".join(measured),
                "dataset_cells": "|".join(cells),
            }
        )
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# G. prior feature blocks
# --------------------------------------------------------------------------


def build_feature_blocks(
    vocabulary: Sequence[str],
    control_means: np.ndarray,
    genes: Sequence[str],
    contexts: Sequence[str],
) -> dict[str, pf.FeatureBlock]:
    blocks: dict[str, pf.FeatureBlock] = {}

    print("  basal ...", flush=True)
    blocks["basal"] = pf.basal_features(vocabulary, control_means, genes, context_names=contexts)

    for name in ("hallmark", "reactome"):
        print(f"  {name} ...", flush=True)
        src = catalogue.source(name)
        blocks[name] = pf.pathway_features(vocabulary, ROOT / src.paths[0], name=name)

    print("  string ...", flush=True)
    src = catalogue.source("string")
    blocks["string"] = pf.string_features(vocabulary, ROOT / src.paths[0], ROOT / src.paths[1])

    print("  depmap ...", flush=True)
    src = catalogue.source("depmap")
    blocks["depmap"] = pf.depmap_features(vocabulary, ROOT / src.paths[0])

    return blocks


def prepare_matrix(block: pf.FeatureBlock) -> np.ndarray:
    """Standardise a block and make its missingness explicit.

    Uncovered rows are filled with the covered rows' mean -- which after
    standardisation is zero, i.e. 'no information' -- and an indicator column
    records that the fill happened, so a model can learn that absence is itself
    a signal rather than being told a false value.
    """
    x = np.array(block.matrix, dtype=np.float64, copy=True)
    covered = block.covered
    if covered.any():
        mean = np.nanmean(x[covered], axis=0)
        std = np.nanstd(x[covered], axis=0)
    else:
        mean = np.zeros(x.shape[1])
        std = np.ones(x.shape[1])
    std = np.where(std > 0, std, 1.0)
    x = np.where(np.isfinite(x), x, mean[None, :])
    x[~covered] = mean[None, :]
    x = (x - mean[None, :]) / std[None, :]
    return np.column_stack([x, covered.astype(np.float64)])


# --------------------------------------------------------------------------
# E. feasible context main effect
# --------------------------------------------------------------------------


def evaluate_main_effect(
    delta: np.ndarray, control_means: np.ndarray, contexts: Sequence[str]
) -> pd.DataFrame:
    n_ctx = delta.shape[0]
    rows = []
    for target in range(n_ctx):
        sources = [c for c in range(n_ctx) if c != target]
        truth = cme.oracle_main_effect(delta, target)
        for name, fn in cme.FEASIBLE_ESTIMATORS.items():
            pred = fn(delta, sources, control_means=control_means, target=target)
            metrics = up.evaluate_predictions(pred[None, :], truth[None, :])
            rows.append(
                {
                    "target_context": contexts[target],
                    "estimator": name,
                    "pearson": float(metrics["pearson"][0]),
                    "cosine": float(metrics["cosine"][0]),
                    "unexplained_fraction": float(metrics["unexplained_fraction"][0]),
                    "norm_ratio": float(np.linalg.norm(pred) / np.linalg.norm(truth)),
                }
            )
        # The zero prediction, for reference: this is what "emit the control"
        # does, and the VCC scale makes it strictly worse than the baseline.
        metrics = up.evaluate_predictions(np.zeros((1, truth.size)), truth[None, :])
        rows.append(
            {
                "target_context": contexts[target],
                "estimator": "Z_zero_response",
                "pearson": np.nan,
                "cosine": 0.0,
                "unexplained_fraction": float(metrics["unexplained_fraction"][0]),
                "norm_ratio": 0.0,
            }
        )
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# F / H / I. the two-axis benchmark
# --------------------------------------------------------------------------


def run_benchmark(
    delta: np.ndarray,
    control_means: np.ndarray,
    contexts: Sequence[str],
    matrices: dict[str, np.ndarray],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    n_ctx, n_pert, _ = delta.shape
    splits = up.make_two_axis_splits(n_ctx, n_pert, n_folds=N_FOLDS, seed=SEED)

    summary_rows: list[dict] = []
    per_pert_rows: list[dict] = []

    for family, features in matrices.items():
        n_comp = min(PCA_COMPONENTS, features.shape[1])
        for split in splits:
            test = list(split.test_perturbations)
            train = list(split.train_perturbations)
            _, mu_train = up.conserved_beta(delta, split.source_contexts, train)

            for est_name, estimator in up.ESTIMATORS.items():
                kwargs: dict = {}
                if est_name == "U2_knn":
                    kwargs["k"] = KNN_K
                if est_name == "U3_ridge":
                    kwargs["alpha"] = RIDGE_ALPHA
                beta_hat = up.fit_predict_split(
                    delta, features, split, estimator, n_components=n_comp, **kwargs
                )

                if split.regime is up.Regime.P1:
                    # Truth is the conserved effect, in the training frame.
                    truth = delta[:, test].mean(axis=0) - mu_train
                    metrics = up.evaluate_predictions(beta_hat, truth)
                    eval_context = "all"
                    response_metrics = None
                else:
                    target = split.target_context
                    sources = list(split.source_contexts)
                    # Scientific view: the beta part, with the oracle main effect removed.
                    truth = delta[target][test] - cme.oracle_main_effect(delta, target)
                    metrics = up.evaluate_predictions(beta_hat, truth)
                    # Deployable view: everything a submission would have to emit,
                    # using only a FEASIBLE main-effect estimate.
                    m_hat = cme.basal_weighted_mean(
                        delta, sources, control_means=control_means, target=target
                    )
                    response_metrics = up.evaluate_predictions(
                        beta_hat + m_hat[None, :], delta[target][test]
                    )
                    eval_context = contexts[target]

                row = {
                    "regime": split.regime.value,
                    "prior_family": family,
                    "estimator": est_name,
                    "eval_context": eval_context,
                    "n_test": len(test),
                    "beta_pearson": float(np.nanmedian(metrics["pearson"])),
                    "beta_cosine": float(np.nanmedian(metrics["cosine"])),
                    "beta_unexplained": float(np.nanmedian(metrics["unexplained_fraction"])),
                }
                if response_metrics is not None:
                    row["response_pearson"] = float(np.nanmedian(response_metrics["pearson"]))
                    row["response_unexplained"] = float(
                        np.nanmedian(response_metrics["unexplained_fraction"])
                    )
                summary_rows.append(row)

                for i, p in enumerate(test):
                    per_pert_rows.append(
                        {
                            "regime": split.regime.value,
                            "prior_family": family,
                            "estimator": est_name,
                            "eval_context": eval_context,
                            "perturbation_index": int(p),
                            "beta_pearson": float(metrics["pearson"][i]),
                            "beta_unexplained": float(metrics["unexplained_fraction"][i]),
                        }
                    )

    return pd.DataFrame(summary_rows), pd.DataFrame(per_pert_rows)


# --------------------------------------------------------------------------
# J. support distance
# --------------------------------------------------------------------------


def reference_ceilings(delta: np.ndarray, contexts: Sequence[str]) -> pd.DataFrame:
    """What a perfect and a merely-good beta would score in each regime.

    Without these the benchmark numbers float free. Two reference points
    bracket what is achievable:

    ``oracle_beta``
        the conserved effect computed from every context, target included.
        Its residual in P2 is pure ``gamma``, so this is the ceiling any
        beta-predicting model can reach -- and the frozen conclusion is that
        ``gamma`` is not zero-shot predictable.
    ``seen_perturbation_transfer``
        the conserved effect from the SOURCE contexts only. This is what the
        four-context programme achieved for a perturbation measured elsewhere,
        so it is the Tier-2 ceiling and the thing Tier-0 must be compared with.
    """
    n_ctx, n_pert, _ = delta.shape
    splits = up.make_two_axis_splits(n_ctx, n_pert, n_folds=N_FOLDS, seed=SEED)
    rows = []
    for split in splits:
        if split.regime is not up.Regime.P2:
            continue
        target = split.target_context
        test = list(split.test_perturbations)
        truth = delta[target][test] - cme.oracle_main_effect(delta, target)

        oracle = delta[:, test].mean(axis=0) - delta.mean(axis=(0, 1))
        sources = list(split.source_contexts)
        transfer = delta[np.ix_(sources, test)].mean(axis=0) - delta[sources].mean(axis=(0, 1))
        # Frozen conclusion 1: the point-prediction baseline is SCALE-CALIBRATED
        # conserved transfer, with the scalar fitted inner leave-one-source-out
        # on the sources alone. Comparing against raw transfer would understate
        # the Tier-2 ceiling, because raw transfer overshoots in magnitude.
        scale = foundations.fit_response_shrinkage(delta, sources)

        for name, pred in (
            ("oracle_beta", oracle),
            ("seen_perturbation_transfer_raw", transfer),
            ("seen_perturbation_transfer_scaled", scale * transfer),
        ):
            m = up.evaluate_predictions(pred, truth)
            rows.append(
                {
                    "eval_context": contexts[target],
                    "reference": name,
                    "beta_pearson": float(np.nanmedian(m["pearson"])),
                    "beta_unexplained": float(np.nanmedian(m["unexplained_fraction"])),
                    "scale": float(scale),
                }
            )
    return pd.DataFrame(rows)


def support_scores(features: np.ndarray, train: Sequence[int], query: Sequence[int]) -> np.ndarray:
    """Max cosine similarity from each query row to any training row."""
    a = np.asarray(features, dtype=np.float64)
    an = np.linalg.norm(a, axis=1, keepdims=True)
    unit = a / np.where(an > 0, an, 1.0)
    return (unit[list(query)] @ unit[list(train)].T).max(axis=1)


def main() -> None:
    started = time.time()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    delta, control_means, contexts, perts, genes = load_canonical()
    print(f"canonical tensor {delta.shape}  contexts={contexts}")

    with INVENTORY.open() as fh:
        inventory = json.load(fh)

    rule("D. SUPPORT TIERS FOR THE 300 ARC TARGETS")
    tiers = support_tiers(inventory)
    SPLITS.mkdir(parents=True, exist_ok=True)
    tiers.to_csv(SPLITS / "arc_target_support_v1.csv", index=False)
    counts = tiers["support_tier"].value_counts().sort_index()
    for tier in (0, 1, 2):
        n = int(counts.get(tier, 0))
        print(f"  TIER {tier}: {n:>3} / 300  ({100 * n / 300:.1f}%)")
    print(f"  written to {SPLITS / 'arc_target_support_v1.csv'}")

    rule("G. PRIOR SOURCES")
    arc_targets = list(tiers["arc_target"])
    vocabulary = sorted(set(perts) | set(arc_targets))
    print(
        f"  feature vocabulary: {len(vocabulary)} genes "
        f"({len(perts)} trainable perturbations + {len(arc_targets)} Arc targets)"
    )
    blocks = build_feature_blocks(vocabulary, control_means, genes, contexts)

    vocab_index = pd.Index(vocabulary)
    train_rows = vocab_index.get_indexer(pd.Index(perts))
    arc_rows = vocab_index.get_indexer(pd.Index(arc_targets))

    cov_rows = []
    for name, block in blocks.items():
        src = catalogue.source(name)
        cov_rows.append(
            {
                "prior_family": name,
                "version": src.version,
                "licence": src.licence,
                "leakage": src.leakage.value,
                "n_features": block.n_features,
                "coverage_all": block.coverage,
                "coverage_training_perturbations": float(block.covered[train_rows].mean()),
                "coverage_arc_targets": float(block.covered[arc_rows].mean()),
                "missing_fraction": float(np.isnan(block.matrix).mean()),
            }
        )
    coverage = pd.DataFrame(cov_rows)
    print(coverage.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    coverage.to_csv(OUTDIR / "prior_coverage.csv", index=False)

    matrices = {name: prepare_matrix(block) for name, block in blocks.items()}
    combined = np.column_stack([matrices[n] for n in sorted(matrices)])
    matrices["combined"] = combined
    print(f"\n  combined feature width: {combined.shape[1]}")

    rule("E. FEASIBLE CONTEXT MAIN EFFECT (public leave-one-context-out)")
    main_effect = evaluate_main_effect(delta, control_means, contexts)
    print(main_effect.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    main_effect.to_csv(OUTDIR / "context_main_effect.csv", index=False)

    rule("F / H / I. TWO-AXIS UNSEEN-PERTURBATION BENCHMARK")
    train_matrices = {k: v[train_rows] for k, v in matrices.items()}
    summary, per_pert = run_benchmark(delta, control_means, contexts, train_matrices)
    summary.to_csv(OUTDIR / "benchmark_summary.csv", index=False)
    per_pert.to_csv(OUTDIR / "benchmark_per_perturbation.csv", index=False)

    for regime in ("P1", "P2"):
        print(f"\n  --- {regime} ---")
        block = summary[summary["regime"] == regime]
        cols = ["beta_pearson", "beta_unexplained"]
        if regime == "P2":
            cols += ["response_pearson", "response_unexplained"]
        pivot = block.groupby(["prior_family", "estimator"])[cols].median().reset_index()
        print(pivot.to_string(index=False, float_format=lambda v: f"{v:+.4f}"))

    rule("CEILINGS: WHAT A PERFECT BETA WOULD SCORE IN P2")
    ceilings = reference_ceilings(delta, contexts)
    print(
        ceilings.groupby("reference")[["beta_pearson", "beta_unexplained"]]
        .median()
        .to_string(float_format=lambda v: f"{v:+.4f}")
    )
    ceilings.to_csv(OUTDIR / "ceilings.csv", index=False)
    print("\n  The residual of oracle_beta in P2 is pure gamma. Nothing that")
    print("  predicts only beta can go below it, however good the prior.")

    rule("J. SUPPORT DISTANCE AND TIER-0 COVERAGE")
    dist_rows = []
    for family, matrix in matrices.items():
        arc_support = support_scores(matrix, train_rows, arc_rows)
        held_support = support_scores(
            matrix, train_rows[: len(train_rows) // 2], train_rows[len(train_rows) // 2 :]
        )
        for tier in (0, 1, 2):
            mask = (tiers["support_tier"] == tier).to_numpy()
            if not mask.any():
                continue
            dist_rows.append(
                {
                    "prior_family": family,
                    "group": f"arc_tier{tier}",
                    "n": int(mask.sum()),
                    "median_max_similarity": float(np.median(arc_support[mask])),
                    "q10": float(np.quantile(arc_support[mask], 0.10)),
                }
            )
        dist_rows.append(
            {
                "prior_family": family,
                "group": "held_out_training_gene",
                "n": len(held_support),
                "median_max_similarity": float(np.median(held_support)),
                "q10": float(np.quantile(held_support, 0.10)),
            }
        )
    support = pd.DataFrame(dist_rows)
    print(support.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    support.to_csv(OUTDIR / "support_distance.csv", index=False)

    print("\n  Arc targets with NO representation at all in each prior:")
    zero_rows = []
    for family, matrix in matrices.items():
        arc_support = support_scores(matrix, train_rows, arc_rows)
        for tier in (0, 1, 2):
            mask = (tiers["support_tier"] == tier).to_numpy()
            zero_rows.append(
                {
                    "prior_family": family,
                    "tier": tier,
                    "n": int(mask.sum()),
                    "n_zero_support": int((arc_support[mask] <= 0).sum()),
                    "pct_zero_support": float(100 * (arc_support[mask] <= 0).mean()),
                }
            )
    zero = pd.DataFrame(zero_rows)
    print(zero.to_string(index=False, float_format=lambda v: f"{v:.1f}"))
    zero.to_csv(OUTDIR / "zero_support.csv", index=False)

    rule("J7. IS THERE A USABLE UNCERTAINTY SCORE FOR AN UNSEEN PERTURBATION?")
    agreement = {}
    splits_p2 = [
        s
        for s in up.make_two_axis_splits(delta.shape[0], delta.shape[1], n_folds=N_FOLDS, seed=SEED)
        if s.regime is up.Regime.P2
    ]
    for family, matrix in train_matrices.items():
        per_index = {}
        for split in splits_p2:
            if split.target_context != 0:
                continue
            beta_train, _ = up.conserved_beta(
                delta, split.source_contexts, split.train_perturbations
            )
            vals = up.neighbour_agreement(
                matrix[list(split.train_perturbations)],
                beta_train,
                matrix[list(split.test_perturbations)],
                k=KNN_K,
            )
            for idx, v in zip(split.test_perturbations, vals, strict=True):
                per_index[idx] = v
        agreement[family] = per_index

    acc_rows = []
    for family in matrices:
        scores = support_scores(matrices[family], train_rows, train_rows)
        for regime in ("P1", "P2"):
            block = per_pert[
                (per_pert["regime"] == regime)
                & (per_pert["prior_family"] == family)
                & (per_pert["estimator"] == "U2_knn")
            ]
            if block.empty:
                continue
            agg = block.groupby("perturbation_index")["beta_pearson"].median()
            s_vec = scores[agg.index.to_numpy()]
            ok = np.isfinite(agg.to_numpy()) & np.isfinite(s_vec)
            if ok.sum() < 20:
                continue
            rho = stats.spearmanr(s_vec[ok], agg.to_numpy()[ok]).statistic
            row = {
                "prior_family": family,
                "regime": regime,
                "n": int(ok.sum()),
                "spearman_support_vs_accuracy": float(rho),
                "spearman_neighbour_agreement_vs_accuracy": np.nan,
            }
            if family in agreement:
                a_vec = np.array([agreement[family].get(i, np.nan) for i in agg.index.to_numpy()])
                ok2 = np.isfinite(agg.to_numpy()) & np.isfinite(a_vec)
                if ok2.sum() >= 20:
                    row["spearman_neighbour_agreement_vs_accuracy"] = float(
                        stats.spearmanr(a_vec[ok2], agg.to_numpy()[ok2]).statistic
                    )
            acc_rows.append(row)
    uncertainty = pd.DataFrame(acc_rows)
    print(uncertainty.to_string(index=False, float_format=lambda v: f"{v:+.4f}"))
    uncertainty.to_csv(OUTDIR / "support_vs_accuracy.csv", index=False)

    summary_json = {
        "tier_counts": {str(k): int(v) for k, v in counts.items()},
        "n_folds": N_FOLDS,
        "seed": SEED,
        "knn_k": KNN_K,
        "ridge_alpha": RIDGE_ALPHA,
        "pca_components": PCA_COMPONENTS,
        "runtime_minutes": (time.time() - started) / 60.0,
    }
    with (OUTDIR / "summary.json").open("w") as fh:
        json.dump(summary_json, fh, indent=2, sort_keys=True)
    print(f"\nWrote {OUTDIR}   [{(time.time() - started) / 60:.1f} min]")


if __name__ == "__main__":
    main()
