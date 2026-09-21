"""Scoring the frozen unseen-perturbation predictor on a delta matrix.

The frozen external protocol (:mod:`virtual_cell.modelling.external_benchmark`)
reads cells from an ``.h5ad`` and pseudobulks them. Feng et al. ship per-line
**log fold changes** instead, so the cell-level entry point does not apply. This
module provides the same scoring over an already-formed delta matrix.

It does not redefine a single constant or estimator: ``PRIOR_FAMILY``,
``KNN_K``, ``RIDGE_ALPHA``, ``PCA_COMPONENTS``, ``prepare_matrix`` and the
estimators are imported from the frozen modules. The equivalence is checked
rather than claimed -- ``scripts/run_feng_multicontext.py`` first feeds arch1's
pseudobulk through this path and requires the frozen arch1 numbers back before
any Feng line is scored.

One representation caveat, stated because it cannot be removed: our source
deltas are differences of log-normalised pseudobulk means, whereas Feng's are
log fold changes from their own model. Both are log-space contrasts against a
non-targeting control, but they are not the identical estimator, and a
correlation computed across the two inherits that difference.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from virtual_cell.analysis import foundations
from virtual_cell.modelling import unseen_perturbation as up
from virtual_cell.modelling.external_benchmark import (
    KNN_K,
    PCA_COMPONENTS,
    PRIOR_FAMILY,
    RIDGE_ALPHA,
    prepare_matrix,
)
from virtual_cell.priors import catalogue
from virtual_cell.priors import features as pf

__all__ = [
    "SourceBlock",
    "build_source_block",
    "string_matrix",
    "score_unseen",
    "score_measured_transfer",
    "cross_line_agreement",
]


class SourceBlock:
    """The source tensor restricted to a gene axis, plus what the protocol derives."""

    def __init__(
        self,
        delta: np.ndarray,
        source_genes: Sequence[str],
        source_perturbations: Sequence[str],
        genes: Sequence[str],
    ) -> None:
        keep = pd.Index(source_genes).get_indexer(pd.Index(genes))
        if (keep < 0).any():
            raise ValueError("every gene must be present on the source axis")
        self.genes = list(genes)
        self.perturbations = list(source_perturbations)
        self.delta = np.asarray(delta, dtype=np.float64)[:, :, keep]
        self.sources = list(range(self.delta.shape[0]))
        self.beta = self.delta.mean(axis=0) - self.delta.mean(axis=(0, 1))
        self.scale = foundations.fit_response_shrinkage(self.delta, self.sources)


def build_source_block(
    delta: np.ndarray,
    source_genes: Sequence[str],
    source_perturbations: Sequence[str],
    genes: Sequence[str],
) -> SourceBlock:
    return SourceBlock(delta, source_genes, source_perturbations, genes)


def string_matrix(vocabulary: Sequence[str], root: Path) -> tuple[np.ndarray, np.ndarray]:
    """The frozen STRING representation over a gene vocabulary."""
    src = catalogue.source(PRIOR_FAMILY)
    block = pf.string_features(vocabulary, root / src.paths[0], root / src.paths[1])
    return prepare_matrix(block), block.covered


def score_unseen(
    observed: np.ndarray,
    test_perturbations: Sequence[str],
    block: SourceBlock,
    matrix: np.ndarray,
    vocabulary: Sequence[str],
    *,
    main_effect: np.ndarray | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run the frozen estimators on one target context's unseen perturbations.

    ``observed`` is ``(n_test, n_genes)`` on ``block.genes``. The target's own
    mean over perturbations is removed to form the beta-level truth, exactly as
    the frozen protocol does. ``main_effect`` enables the response-level view;
    pass ``None`` when no feasible estimate exists for this context.
    """
    observed = np.asarray(observed, dtype=np.float64)
    vocab = pd.Index(vocabulary)
    train_rows = vocab.get_indexer(pd.Index(block.perturbations))
    test_rows = vocab.get_indexer(pd.Index(test_perturbations))

    train_f, test_f = matrix[train_rows], matrix[test_rows]
    pca = pf.fit_pca(train_f, min(PCA_COMPONENTS, train_f.shape[1]))
    train_p, test_p = pf.apply_pca(train_f, pca), pf.apply_pca(test_f, pca)

    oracle = observed.mean(axis=0)
    truth = observed - oracle[None, :]

    rows = []
    per_pert: dict[str, np.ndarray] = {}
    for name, estimator in up.ESTIMATORS.items():
        kwargs: dict = {}
        if name == "U2_knn":
            kwargs["k"] = KNN_K
        if name == "U3_ridge":
            kwargs["alpha"] = RIDGE_ALPHA
        beta_hat = estimator(train_p, block.beta, test_p, **kwargs)
        metrics = up.evaluate_predictions(beta_hat, truth)
        row = {
            "estimator": name,
            "n_test": len(test_perturbations),
            "beta_pearson": float(np.nanmedian(metrics["pearson"])),
            "beta_spearman": _median_spearman(beta_hat, truth),
            "beta_cosine": float(np.nanmedian(metrics["cosine"])),
            "beta_unexplained": float(np.nanmedian(metrics["unexplained_fraction"])),
        }
        if main_effect is not None:
            resp = up.evaluate_predictions(beta_hat + main_effect[None, :], observed)
            row["response_pearson"] = float(np.nanmedian(resp["pearson"]))
            row["response_unexplained"] = float(np.nanmedian(resp["unexplained_fraction"]))
        rows.append(row)
        if name == "U2_knn":
            per_pert["knn_pearson"] = metrics["pearson"]
            per_pert["knn_unexplained"] = metrics["unexplained_fraction"]

    agreement = up.neighbour_agreement(train_p, block.beta, test_p, k=KNN_K)
    support = _max_cosine(test_p, train_p)
    detail = pd.DataFrame(
        {
            "perturbation": list(test_perturbations),
            "neighbour_agreement": agreement,
            "support_distance": support,
            **per_pert,
        }
    )
    return pd.DataFrame(rows), detail


def score_measured_transfer(
    observed: np.ndarray, seen_perturbations: Sequence[str], block: SourceBlock
) -> dict[str, float]:
    """The positive control: transfer a perturbation measured in the sources.

    Without this, a failure of the unseen arm cannot be attributed to the prior
    rather than to the target measurement -- which is exactly what went wrong
    on ``kaden25rpe1``.
    """
    observed = np.asarray(observed, dtype=np.float64)
    pos = pd.Index(block.perturbations).get_indexer(pd.Index(seen_perturbations))
    if (pos < 0).any():
        raise ValueError("every seen perturbation must be a source perturbation")
    transfer = block.delta[:, pos].mean(axis=0) - block.delta.mean(axis=(0, 1))
    truth = observed - observed.mean(axis=0)[None, :]
    raw = up.evaluate_predictions(transfer, truth)
    scaled = up.evaluate_predictions(block.scale * transfer, truth)
    return {
        "n_seen": len(seen_perturbations),
        "scale": float(block.scale),
        "raw_pearson": float(np.nanmedian(raw["pearson"])),
        "raw_unexplained": float(np.nanmedian(raw["unexplained_fraction"])),
        "scaled_pearson": float(np.nanmedian(scaled["pearson"])),
        "scaled_spearman": _median_spearman(block.scale * transfer, truth),
        "scaled_cosine": float(np.nanmedian(scaled["cosine"])),
        "scaled_unexplained": float(np.nanmedian(scaled["unexplained_fraction"])),
    }


def cross_line_agreement(stack: np.ndarray) -> dict[str, np.ndarray]:
    """How much a perturbation's response actually changes between lines.

    ``stack`` is ``(n_lines, n_perturbations, n_genes)`` with NaN where a line
    does not carry a perturbation. Returns, per perturbation, the mean pairwise
    correlation between lines and the fraction of energy that is line-specific.
    """
    n_lines, n_perts, _ = stack.shape
    agreement = np.full(n_perts, np.nan)
    line_specific = np.full(n_perts, np.nan)
    for p in range(n_perts):
        block = stack[:, p, :]
        ok = np.isfinite(block).all(axis=1)
        if ok.sum() < 2:
            continue
        present = block[ok]
        centred = present - present.mean(axis=1, keepdims=True)
        norm = np.linalg.norm(centred, axis=1)
        good = norm > 0
        if good.sum() < 2:
            continue
        unit = centred[good] / norm[good][:, None]
        gram = unit @ unit.T
        off = gram[~np.eye(gram.shape[0], dtype=bool)]
        agreement[p] = float(off.mean())
        consensus = present[good].mean(axis=0)
        total = float((present[good] ** 2).sum())
        if total > 0:
            line_specific[p] = float(((present[good] - consensus[None, :]) ** 2).sum() / total)
    del n_lines
    return {"cross_line_agreement": agreement, "line_specific_fraction": line_specific}


def _max_cosine(query: np.ndarray, reference: np.ndarray) -> np.ndarray:
    qn = np.linalg.norm(query, axis=1, keepdims=True)
    rn = np.linalg.norm(reference, axis=1, keepdims=True)
    q = query / np.where(qn > 0, qn, 1.0)
    r = reference / np.where(rn > 0, rn, 1.0)
    return (q @ r.T).max(axis=1)


def _median_spearman(predicted: np.ndarray, observed: np.ndarray) -> float:
    from scipy import stats

    values = [
        stats.spearmanr(p, o).statistic
        for p, o in zip(np.asarray(predicted), np.asarray(observed), strict=True)
    ]
    return float(np.nanmedian(values)) if values else float("nan")
