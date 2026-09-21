"""Predicting a conserved perturbation effect for a gene never perturbed.

The four-context programme predicted ``beta_p`` by averaging measurements of
that same perturbation in other contexts. For 214 of Arc's 300 targets no such
measurement exists anywhere, so ``beta_p`` has to come from what is known about
the gene itself. This module holds the estimators and the two-axis protocol
that tests whether that is possible at all.

Two axes are held out independently:

* **context** -- the target context's responses are hidden entirely, as in the
  existing leave-one-context-out work;
* **perturbation** -- the test perturbation is removed from *every* context, so
  no measurement of it exists anywhere in training.

Arc's validation panel is predominantly both at once. A protocol that holds out
only one axis will report a number that does not transfer, so
:func:`nested_evaluation` refuses to run without both.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

__all__ = [
    "Regime",
    "TwoAxisSplit",
    "perturbation_folds",
    "make_two_axis_splits",
    "conserved_beta",
    "context_main_effect",
    "predict_zero",
    "predict_nearest_neighbour",
    "predict_knn",
    "predict_ridge",
    "ESTIMATORS",
    "fit_predict_split",
    "neighbour_agreement",
    "evaluate_predictions",
]

from enum import StrEnum


class Regime(StrEnum):
    """Which axes a split holds out."""

    #: Held-out perturbation, seen context.
    P1 = "P1"
    #: Held-out perturbation and held-out context. The Arc-like regime.
    P2 = "P2"


@dataclass(frozen=True)
class TwoAxisSplit:
    """One evaluation fold, with both axes named explicitly."""

    regime: Regime
    #: Index of the held-out context, or ``None`` in :attr:`Regime.P1`.
    target_context: int | None
    #: Contexts whose responses training may read.
    source_contexts: tuple[int, ...]
    #: Perturbation indices removed from every context during training.
    test_perturbations: tuple[int, ...]
    #: Perturbation indices training may read.
    train_perturbations: tuple[int, ...]
    #: Contexts the predictions are scored in.
    eval_contexts: tuple[int, ...]

    def __post_init__(self) -> None:
        if set(self.test_perturbations) & set(self.train_perturbations):
            raise ValueError("a perturbation cannot be both train and test")
        if self.regime is Regime.P2:
            if self.target_context is None:
                raise ValueError("P2 requires a held-out target context")
            if self.target_context in self.source_contexts:
                raise ValueError("the target context must not be a source context")


def perturbation_folds(
    n_perturbations: int, n_folds: int = 5, *, seed: int = 0
) -> list[np.ndarray]:
    """Partition perturbation indices into ``n_folds`` disjoint groups.

    The partition is global: a perturbation in a fold is withheld from every
    context at once, which is what makes it unseen rather than merely
    unseen-here.
    """
    if n_folds < 2:
        raise ValueError("n_folds must be at least 2")
    order = np.random.default_rng(seed).permutation(n_perturbations)
    return [np.sort(part) for part in np.array_split(order, n_folds)]


def make_two_axis_splits(
    n_contexts: int,
    n_perturbations: int,
    *,
    n_folds: int = 5,
    seed: int = 0,
    regimes: Sequence[Regime] = (Regime.P1, Regime.P2),
) -> list[TwoAxisSplit]:
    """Enumerate every fold of the requested regimes."""
    folds = perturbation_folds(n_perturbations, n_folds, seed=seed)
    contexts = tuple(range(n_contexts))
    splits: list[TwoAxisSplit] = []
    for test in folds:
        train = tuple(int(p) for p in range(n_perturbations) if p not in set(test.tolist()))
        test_t = tuple(int(p) for p in test)
        if Regime.P1 in regimes:
            splits.append(
                TwoAxisSplit(
                    regime=Regime.P1,
                    target_context=None,
                    source_contexts=contexts,
                    test_perturbations=test_t,
                    train_perturbations=train,
                    eval_contexts=contexts,
                )
            )
        if Regime.P2 in regimes:
            for target in contexts:
                sources = tuple(c for c in contexts if c != target)
                splits.append(
                    TwoAxisSplit(
                        regime=Regime.P2,
                        target_context=target,
                        source_contexts=sources,
                        test_perturbations=test_t,
                        train_perturbations=train,
                        eval_contexts=(target,),
                    )
                )
    return splits


# --------------------------------------------------------------------------
# the quantities being predicted
# --------------------------------------------------------------------------


def conserved_beta(
    delta: np.ndarray, contexts: Sequence[int], perturbations: Sequence[int]
) -> tuple[np.ndarray, np.ndarray]:
    """``(beta, mu)`` estimated from a sub-block of the response tensor.

    ``delta`` is ``(n_contexts, n_perturbations, n_genes)``. Both ``mu`` and
    ``beta`` are computed from the named contexts and perturbations only, so a
    held-out perturbation contributes to neither -- including through ``mu``,
    which is the easy leak to miss.
    """
    block = np.asarray(delta)[np.ix_(list(contexts), list(perturbations))]
    mu = block.mean(axis=(0, 1))
    beta = block.mean(axis=0) - mu
    return beta, mu


def context_main_effect(
    delta: np.ndarray, context: int, perturbations: Sequence[int]
) -> np.ndarray:
    """``m_c = mean_p delta[c, p]`` over the named perturbations.

    This is the ORACLE context main effect. For Arc's A, B and C it is hidden,
    so it may be used as an evaluation target and never as a model input.
    """
    return np.asarray(delta)[context][list(perturbations)].mean(axis=0)


# --------------------------------------------------------------------------
# estimators
# --------------------------------------------------------------------------


def predict_zero(
    train_features: np.ndarray, train_beta: np.ndarray, test_features: np.ndarray, **_: object
) -> np.ndarray:
    """``U0``: no perturbation-specific effect."""
    del train_features
    return np.zeros((test_features.shape[0], train_beta.shape[1]))


def _cosine_similarity(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    an = np.linalg.norm(a, axis=1, keepdims=True)
    bn = np.linalg.norm(b, axis=1, keepdims=True)
    a = a / np.where(an > 0, an, 1.0)
    b = b / np.where(bn > 0, bn, 1.0)
    return a @ b.T


def predict_nearest_neighbour(
    train_features: np.ndarray, train_beta: np.ndarray, test_features: np.ndarray, **_: object
) -> np.ndarray:
    """``U1``: copy the conserved effect of the most similar training gene."""
    sim = _cosine_similarity(test_features, train_features)
    return train_beta[np.argmax(sim, axis=1)]


def predict_knn(
    train_features: np.ndarray,
    train_beta: np.ndarray,
    test_features: np.ndarray,
    *,
    k: int = 10,
    **_: object,
) -> np.ndarray:
    """``U2``: similarity-weighted average over the ``k`` nearest training genes.

    Negative similarities are clipped away rather than used as negative
    weights: a gene being *unlike* the target says nothing about which
    direction its response should be subtracted in.
    """
    sim = _cosine_similarity(test_features, train_features)
    k = int(min(k, train_features.shape[0]))
    out = np.zeros((test_features.shape[0], train_beta.shape[1]))
    for i in range(test_features.shape[0]):
        top = np.argpartition(-sim[i], k - 1)[:k]
        w = np.clip(sim[i][top], 0.0, None)
        total = w.sum()
        if total <= 0:
            continue
        out[i] = (w / total) @ train_beta[top]
    return out


def predict_ridge(
    train_features: np.ndarray,
    train_beta: np.ndarray,
    test_features: np.ndarray,
    *,
    alpha: float = 1.0,
    **_: object,
) -> np.ndarray:
    """``U3``: multi-output ridge from the feature representation to ``beta``.

    Solved in whichever of the primal or dual form is smaller, which matters
    because the feature blocks here are routinely wider than they are tall.
    """
    x = np.asarray(train_features, dtype=np.float64)
    y = np.asarray(train_beta, dtype=np.float64)
    xm, ym = x.mean(axis=0), y.mean(axis=0)
    xc, yc = x - xm, y - ym
    n, d = xc.shape
    if d <= n:
        gram = xc.T @ xc + alpha * np.eye(d)
        coef = np.linalg.solve(gram, xc.T @ yc)
    else:
        gram = xc @ xc.T + alpha * np.eye(n)
        coef = xc.T @ np.linalg.solve(gram, yc)
    return (np.asarray(test_features, dtype=np.float64) - xm) @ coef + ym


#: The low-capacity estimators, in the order the protocol reports them.
ESTIMATORS = {
    "U0_zero": predict_zero,
    "U1_nearest": predict_nearest_neighbour,
    "U2_knn": predict_knn,
    "U3_ridge": predict_ridge,
}


def neighbour_agreement(
    train_features: np.ndarray,
    train_beta: np.ndarray,
    test_features: np.ndarray,
    *,
    k: int = 25,
) -> np.ndarray:
    """How much the ``k`` nearest training genes agree with each other.

    The frozen confidence result for a *seen* perturbation was raw source
    agreement: how well the contexts that measured it agree. A Tier-0
    perturbation has no sources, so that statistic is undefined. This is its
    natural analogue -- not agreement between measurements of the same gene,
    but agreement between the genes standing in for it. If the neighbours a
    prediction is averaging disagree, the average is not worth trusting.

    Returned as the mean pairwise cosine between the neighbours' conserved
    effects, which is bounded and needs no reference measurement.
    """
    sim = _cosine_similarity(test_features, train_features)
    k = int(min(k, train_features.shape[0]))
    beta = np.asarray(train_beta, dtype=np.float64)
    norm = np.linalg.norm(beta, axis=1, keepdims=True)
    unit = beta / np.where(norm > 0, norm, 1.0)

    out = np.full(test_features.shape[0], np.nan)
    for i in range(test_features.shape[0]):
        top = np.argpartition(-sim[i], k - 1)[:k]
        block = unit[top]
        gram = block @ block.T
        off = gram[~np.eye(k, dtype=bool)]
        if off.size:
            out[i] = float(off.mean())
    return out


# --------------------------------------------------------------------------
# the one entry point that touches the response tensor
# --------------------------------------------------------------------------


def fit_predict_split(
    delta: np.ndarray,
    features: np.ndarray,
    split: TwoAxisSplit,
    estimator,
    *,
    n_components: int | None = None,
    **kwargs: object,
) -> np.ndarray:
    """Fit on a split's training block and predict ``beta`` for its test genes.

    This is the *only* function in the protocol that reads ``delta``, and it
    reads exactly ``delta[source_contexts][train_perturbations]``. That
    narrowness is the point: it makes the leakage claim checkable by
    corrupting everything outside that block and requiring the output to be
    bit-identical, which is what the test-suite does.

    Any projection is fitted on training rows only, for the same reason.
    """
    beta_train, _ = conserved_beta(delta, split.source_contexts, split.train_perturbations)
    train_features = np.asarray(features, dtype=np.float64)[list(split.train_perturbations)]
    test_features = np.asarray(features, dtype=np.float64)[list(split.test_perturbations)]

    if n_components is not None:
        from virtual_cell.priors.features import apply_pca, fit_pca

        pca = fit_pca(train_features, n_components)
        train_features = apply_pca(train_features, pca)
        test_features = apply_pca(test_features, pca)

    return estimator(train_features, beta_train, test_features, **kwargs)


# --------------------------------------------------------------------------
# evaluation
# --------------------------------------------------------------------------


def _rowwise_pearson(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    ac = a - a.mean(axis=1, keepdims=True)
    bc = b - b.mean(axis=1, keepdims=True)
    num = (ac * bc).sum(axis=1)
    den = np.linalg.norm(ac, axis=1) * np.linalg.norm(bc, axis=1)
    return np.where(den > 0, num / np.where(den > 0, den, 1.0), np.nan)


def evaluate_predictions(predicted: np.ndarray, observed: np.ndarray) -> dict[str, np.ndarray]:
    """Per-row agreement, plus the variance actually explained.

    ``r`` alone can look healthy while the prediction has the wrong scale, so
    the fraction of observed energy left unexplained is reported beside it.
    """
    predicted = np.asarray(predicted, dtype=np.float64)
    observed = np.asarray(observed, dtype=np.float64)
    residual = ((observed - predicted) ** 2).sum(axis=1)
    signal = (observed**2).sum(axis=1)
    safe = np.where(signal > 0, signal, 1.0)
    return {
        "pearson": _rowwise_pearson(predicted, observed),
        "cosine": np.diagonal(_cosine_similarity(predicted, observed)),
        "unexplained_fraction": np.where(signal > 0, residual / safe, np.nan),
    }
