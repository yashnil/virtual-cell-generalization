"""The frozen external-validation protocol, applied to a new target context.

``scripts/run_arch1_external.py`` is frozen, so this module does not import it.
It is a transcription of the same protocol, and the transcription is checked
rather than asserted: ``scripts/run_external_validation.py`` runs it on
``arch1`` first and refuses to continue unless every number matches the frozen
``arch1_external.csv`` to within floating-point tolerance. Only then is a new
context scored.

Nothing here is fitted to a new dataset. The estimator, its hyperparameters,
the prior representation, the preprocessing and the confidence statistic are
all constants, fixed before ``arch1`` was scored and unchanged since.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from virtual_cell.analysis import foundations, loco
from virtual_cell.data import scperteval
from virtual_cell.modelling import context_main_effect as cme
from virtual_cell.modelling import unseen_perturbation as up
from virtual_cell.priors import catalogue
from virtual_cell.priors import features as pf

__all__ = [
    "PRIOR_FAMILY",
    "KNN_K",
    "RIDGE_ALPHA",
    "PCA_COMPONENTS",
    "ExternalResult",
    "prepare_matrix",
    "run_external_benchmark",
    "target_reliability",
]

#: Frozen on the internal benchmark, before any external context was scored.
PRIOR_FAMILY = "string"
KNN_K = 25
RIDGE_ALPHA = 10.0
PCA_COMPONENTS = 64


def prepare_matrix(block: pf.FeatureBlock) -> np.ndarray:
    """Standardise a feature block, with missingness carried as its own column."""
    x = np.array(block.matrix, dtype=np.float64, copy=True)
    covered = block.covered
    mean = np.nanmean(x[covered], axis=0) if covered.any() else np.zeros(x.shape[1])
    std = np.nanstd(x[covered], axis=0) if covered.any() else np.ones(x.shape[1])
    std = np.where(std > 0, std, 1.0)
    x = np.where(np.isfinite(x), x, mean[None, :])
    x[~covered] = mean[None, :]
    return np.column_stack([(x - mean[None, :]) / std[None, :], covered.astype(np.float64)])


def target_reliability(
    path: Path,
    *,
    genes: Sequence[str],
    perturbations: Sequence[str],
    seed: int = 0,
    max_perturbations: int | None = None,
) -> pd.DataFrame:
    """Split-half reliability of a target context's own responses.

    This is the control that decides whether an external benchmark can answer
    anything at all. A negative result on a target whose own responses do not
    reproduce is uninterpretable: the ceiling on any correlation is
    ``sqrt(rho)``, so a dataset with ``rho = 0.12`` caps every method at about
    0.35 however good it is.

    Cells are split in half per perturbation **and** for the control, so both
    halves carry independent control noise -- splitting only the perturbation
    side leaves a shared control component and inflates the estimate.
    """
    rng = np.random.default_rng(seed)
    chosen = list(perturbations)
    if max_perturbations is not None and len(chosen) > max_perturbations:
        chosen = sorted(rng.choice(chosen, size=max_perturbations, replace=False).tolist())

    matrix, group = scperteval.load_perturbation_cells(
        path, genes=genes, perturbations=[*chosen, scperteval.CONTROL_LABEL]
    )
    matrix = matrix.tocsr()
    control_id = len(chosen)

    control_rows = rng.permutation(np.flatnonzero(group == control_id))
    half = len(control_rows) // 2
    ctrl_a = np.asarray(matrix[control_rows[:half]].mean(axis=0)).ravel()
    ctrl_b = np.asarray(matrix[control_rows[half:]].mean(axis=0)).ravel()

    a = np.full((len(chosen), len(genes)), np.nan)
    b = a.copy()
    n_cells = np.zeros(len(chosen), dtype=int)
    for i in range(len(chosen)):
        rows = np.flatnonzero(group == i)
        n_cells[i] = len(rows)
        if len(rows) < 4:
            continue
        shuffled = rng.permutation(rows)
        cut = len(shuffled) // 2
        a[i] = np.asarray(matrix[shuffled[:cut]].mean(axis=0)).ravel() - ctrl_a
        b[i] = np.asarray(matrix[shuffled[cut:]].mean(axis=0)).ravel() - ctrl_b

    ok = np.isfinite(a).all(axis=1)
    rho = np.full(len(chosen), np.nan)
    rho[ok] = loco.split_half_reliability(a[ok], b[ok])
    corrected = loco.spearman_brown(rho)
    magnitude = np.full(len(chosen), np.nan)
    magnitude[ok] = np.linalg.norm((a[ok] + b[ok]) / 2.0, axis=1)

    return pd.DataFrame(
        {
            "perturbation": chosen,
            "n_cells": n_cells,
            "split_half_rho": rho,
            "spearman_brown": corrected,
            "reliability_ceiling": np.sqrt(np.clip(corrected, 0.0, None)),
            "delta_norm": magnitude,
        }
    )


@dataclass
class ExternalResult:
    """Everything one external context produces under the frozen protocol."""

    dataset: str
    n_unseen: int
    n_seen: int
    n_shared_genes: int
    prior_coverage_train: float
    prior_coverage_unseen: float
    main_effect: pd.DataFrame
    oracle_scalar: float
    oracle_scalar_residual: float
    main_effect_estimator: str
    unseen: pd.DataFrame
    seen: pd.DataFrame
    confidence: pd.DataFrame
    shrinkage_scale: float
    unseen_labels: list[str] = field(default_factory=list)


def _eligible(target_perturbations: Sequence[str], seen_anywhere: set[str]) -> list[str]:
    """Perturbations unseen under the same definition ``arch1`` used.

    'Unseen' means absent from **every** source context's perturbation set, so
    the exclusion set is the union of the four, not their intersection. A gene
    measured as an output somewhere is still unseen if it was never knocked
    down.
    """
    return [p for p in target_perturbations if p not in seen_anywhere]


def run_external_benchmark(
    dataset: str,
    *,
    path: Path,
    delta: np.ndarray,
    control_means: np.ndarray,
    source_genes: Sequence[str],
    source_perturbations: Sequence[str],
    seen_anywhere: set[str],
    root: Path,
) -> ExternalResult:
    """Score one external target context under the frozen protocol."""
    target_genes = set(scperteval.read_var_names(path))
    labels = set(map(str, scperteval.cell_labels(path)))
    target_perts = sorted(labels - {"control", "non-targeting"})

    genes = [g for g in source_genes if g in target_genes]
    unseen = _eligible(target_perts, seen_anywhere)
    seen = [p for p in target_perts if p in set(source_perturbations)]

    pseudo = scperteval.pseudobulk(path, genes=genes, perturbations=unseen, context=dataset)
    observed = pseudo.delta

    keep = pd.Index(source_genes).get_indexer(pd.Index(genes))
    delta_k = delta[:, :, keep]
    control_k = control_means[:, keep]

    # The target context's own control cells are the only thing read from it,
    # exactly as Arc permits its controls to be read.
    stacked = np.vstack([control_k, pseudo.control_mean[None, :]])
    target_row = stacked.shape[0] - 1
    sources = list(range(delta_k.shape[0]))
    padded = np.concatenate([delta_k, np.zeros((1, *delta_k.shape[1:]))], axis=0)

    oracle = observed.mean(axis=0)
    rows = []
    for name, fn in cme.FEASIBLE_ESTIMATORS.items():
        pred = fn(padded, sources, control_means=stacked, target=target_row)
        m = up.evaluate_predictions(pred[None, :], oracle[None, :])
        rows.append(
            {
                "estimator": name,
                "pearson": float(m["pearson"][0]),
                "unexplained_fraction": float(m["unexplained_fraction"][0]),
                "norm_ratio": float(np.linalg.norm(pred) / np.linalg.norm(oracle)),
            }
        )
    main_effect = pd.DataFrame(rows)

    pooled = cme.source_pooled_mean(padded, sources)
    s_star = float(pooled @ oracle) / float(pooled @ pooled)
    residual = float(((oracle - s_star * pooled) ** 2).sum() / (oracle**2).sum())

    best_me = str(main_effect.loc[main_effect["unexplained_fraction"].idxmin(), "estimator"])
    m_hat = cme.FEASIBLE_ESTIMATORS[best_me](
        padded, sources, control_means=stacked, target=target_row
    )

    vocabulary = sorted(set(source_perturbations) | set(unseen))
    src = catalogue.source(PRIOR_FAMILY)
    block = pf.string_features(vocabulary, root / src.paths[0], root / src.paths[1])
    matrix = prepare_matrix(block)
    vocab = pd.Index(vocabulary)
    train_rows = vocab.get_indexer(pd.Index(source_perturbations))
    test_rows = vocab.get_indexer(pd.Index(unseen))

    beta_train = delta_k.mean(axis=0) - delta_k.mean(axis=(0, 1))
    scale = foundations.fit_response_shrinkage(delta_k, sources)
    train_f, test_f = matrix[train_rows], matrix[test_rows]
    pca = pf.fit_pca(train_f, min(PCA_COMPONENTS, train_f.shape[1]))
    train_p, test_p = pf.apply_pca(train_f, pca), pf.apply_pca(test_f, pca)

    truth_beta = observed - oracle[None, :]
    results = []
    for est_name, estimator in up.ESTIMATORS.items():
        kwargs: dict = {}
        if est_name == "U2_knn":
            kwargs["k"] = KNN_K
        if est_name == "U3_ridge":
            kwargs["alpha"] = RIDGE_ALPHA
        beta_hat = estimator(train_p, beta_train, test_p, **kwargs)
        b = up.evaluate_predictions(beta_hat, truth_beta)
        r = up.evaluate_predictions(beta_hat + m_hat[None, :], observed)
        results.append(
            {
                "estimator": est_name,
                "beta_pearson": float(np.nanmedian(b["pearson"])),
                "beta_unexplained": float(np.nanmedian(b["unexplained_fraction"])),
                "response_pearson": float(np.nanmedian(r["pearson"])),
                "response_unexplained": float(np.nanmedian(r["unexplained_fraction"])),
            }
        )
    unseen_table = pd.DataFrame(results)

    seen_rows = []
    if seen:
        seen_pb = scperteval.pseudobulk(path, genes=genes, perturbations=seen, context=dataset)
        seen_obs = seen_pb.delta
        pos = pd.Index(source_perturbations).get_indexer(pd.Index(seen))
        transfer = delta_k[:, pos].mean(axis=0) - delta_k.mean(axis=(0, 1))
        seen_truth = seen_obs - seen_obs.mean(axis=0)[None, :]
        raw = up.evaluate_predictions(transfer, seen_truth)
        scaled = up.evaluate_predictions(scale * transfer, seen_truth)
        seen_rows.append(
            {
                "n_seen": len(seen),
                "scale": float(scale),
                "raw_pearson": float(np.nanmedian(raw["pearson"])),
                "raw_unexplained": float(np.nanmedian(raw["unexplained_fraction"])),
                "scaled_pearson": float(np.nanmedian(scaled["pearson"])),
                "scaled_unexplained": float(np.nanmedian(scaled["unexplained_fraction"])),
            }
        )

    agreement = up.neighbour_agreement(train_p, beta_train, test_p, k=KNN_K)
    knn_beta = up.predict_knn(train_p, beta_train, test_p, k=KNN_K)
    per_pert = up.evaluate_predictions(knn_beta, truth_beta)
    norm = np.linalg.norm(train_p, axis=1, keepdims=True)
    unit_train = train_p / np.where(norm > 0, norm, 1.0)
    tnorm = np.linalg.norm(test_p, axis=1, keepdims=True)
    unit_test = test_p / np.where(tnorm > 0, tnorm, 1.0)
    support = (unit_test @ unit_train.T).max(axis=1)

    confidence = pd.DataFrame(
        {
            "perturbation": unseen,
            "neighbour_agreement": agreement,
            "support_distance": support,
            "beta_pearson": per_pert["pearson"],
            "beta_unexplained": per_pert["unexplained_fraction"],
            "string_covered": block.covered[test_rows],
        }
    )

    return ExternalResult(
        dataset=dataset,
        n_unseen=len(unseen),
        n_seen=len(seen),
        n_shared_genes=len(genes),
        prior_coverage_train=float(block.covered[train_rows].mean()),
        prior_coverage_unseen=float(block.covered[test_rows].mean()),
        main_effect=main_effect,
        oracle_scalar=s_star,
        oracle_scalar_residual=residual,
        main_effect_estimator=best_me,
        unseen=unseen_table,
        seen=pd.DataFrame(seen_rows),
        confidence=confidence,
        shrinkage_scale=float(scale),
        unseen_labels=list(unseen),
    )
