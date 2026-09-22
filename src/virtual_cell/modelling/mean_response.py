r"""The tiered mean-response model: ``delta_hat[c,p] = m_hat[c] + w[p] beta_hat[p]``.

This is the model the Arc track builds on, and it is deliberately the smallest
one the frozen evidence supports. Three phases of external validation
established what it may and may not contain:

* **direct perturbation evidence transfers** across contexts (Feng: measured
  transfer positive in 19 of 19 lines, pooled ``r = +0.320``), so a
  perturbation measured somewhere public gets a ``beta_hat``;
* **prior-derived prediction of unseen perturbations does not** (arch1 and all
  19 Feng lines at ``r ~ 0``), so a perturbation measured nowhere gets
  ``beta_hat = 0`` exactly — no STRING, no DepMap, no pathway, no embedding;
* **the score's zero is the context main effect** ``m_c = mu + alpha_c``, which
  Arc hides, so ``m_hat`` must be estimated from source responses and target
  controls alone.

The three tiers are the frozen Arc support classification
(``data/splits/arc_target_support_v1.csv``): Tier 2 is direct evidence in two
or more public contexts, Tier 1 in exactly one, Tier 0 in none.

Why ``beta_hat`` is centred
---------------------------
``beta_hat`` must carry the *perturbation-specific* part and nothing else,
because the context-wide part is ``m_hat``'s job and counting it twice would
inflate every prediction by a copy of the template. So the source responses are
centred over perturbations before they are transferred:

.. math::
    \hat\beta_p = \operatorname{mean}_{s \in S_p}
        \bigl(\delta[s,p] - \operatorname{mean}_q \delta[s,q]\bigr)

With the balanced-design decomposition ``delta[c,p] = mu + alpha_c + beta_p +
gamma[c,p]`` this estimates ``beta_p + mean_{s} gamma[s,p]``: the conserved
effect plus whatever interaction the source contexts happen to share. It is
therefore biased toward the sources and away from the target, which is exactly
what the shrinkage ``w`` exists to absorb.

Why ``w`` is a shrinkage and not a fit
--------------------------------------
``gamma`` is unpredictable zero-shot — that is the programme's central negative
result — so the variance-optimal point prediction of a quantity we can only
estimate up to an unpredictable perturbation is a *shrunk* version of the
estimate, not a corrected one. ``w`` has one degree of freedom per tier, chosen
from a five-point grid by nested validation on source contexts only. Tier 0 is
fixed at 0 and is not tunable.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from itertools import combinations

import numpy as np

from virtual_cell.modelling import context_main_effect as cme

__all__ = [
    "SHRINKAGE_GRID",
    "TIER_SOURCE_COUNT",
    "MAIN_EFFECT_ESTIMATORS",
    "TierWeights",
    "centred_source_beta",
    "tier_beta",
    "source_subsets",
    "predict",
    "fit_main_effect_scale",
    "source_pooled_shrunk_mean",
    "zero_main_effect",
    "select_tier_shrinkage",
]

EPS = 1e-12

#: The conservative shrinkage grid. Deliberately coarse: with four public
#: contexts a finer grid would be fitting noise in the third decimal.
SHRINKAGE_GRID: tuple[float, ...] = (0.0, 0.25, 0.5, 0.75, 1.0)

#: How many public contexts supply ``beta_hat`` in each tier. ``None`` means the
#: tier has no direct evidence at all.
TIER_SOURCE_COUNT: dict[int, int | None] = {2: 2, 1: 1, 0: None}


# --------------------------------------------------------------------------
# m_hat: the feasible context-wide perturbation baseline
# --------------------------------------------------------------------------


def zero_main_effect(
    delta: np.ndarray,
    sources: Sequence[int],
    *,
    control_means: np.ndarray | None = None,
    target: int | None = None,
) -> np.ndarray:
    """``M0``: no context-wide perturbation offset at all.

    The floor. It is what a submission emits if it believes perturbing a cell
    has no average effect, and it is the reference the other three must beat.
    """
    del sources, control_means, target
    return np.zeros(np.asarray(delta).shape[-1], dtype=np.float64)


def fit_main_effect_scale(
    delta: np.ndarray,
    sources: Sequence[int],
    base: Callable[..., np.ndarray],
    *,
    control_means: np.ndarray | None = None,
) -> float:
    """One scalar, fitted leave-one-source-out, rescaling a main-effect estimate.

    For each source in turn, its main effect is predicted from the others by
    ``base`` and the least-squares scale that would have been right is recorded;
    the median is returned. Nothing outside ``sources`` is read, so the scalar
    is available at Arc inference time.
    """
    scales: list[float] = []
    for held in sources:
        rest = [c for c in sources if c != held]
        if not rest:
            continue
        pred = base(delta, rest, control_means=control_means, target=held)
        truth = cme.oracle_main_effect(delta, held)
        denom = float(pred @ pred)
        if denom > EPS:
            scales.append(float(pred @ truth) / denom)
    return float(np.median(scales)) if scales else 1.0


def source_pooled_shrunk_mean(
    delta: np.ndarray,
    sources: Sequence[int],
    *,
    control_means: np.ndarray | None = None,
    target: int | None = None,
) -> np.ndarray:
    """``M3a``: :func:`~virtual_cell.modelling.context_main_effect.source_pooled_mean`
    rescaled by one leave-one-source-out scalar."""
    del target
    scale = fit_main_effect_scale(delta, sources, cme.source_pooled_mean)
    return scale * cme.source_pooled_mean(delta, sources, control_means=control_means)


#: The candidates, in increasing capacity. ``M3b`` is the frozen ``E2`` of
#: :mod:`virtual_cell.modelling.context_main_effect`, reused rather than
#: reimplemented. None of them exceeds one fitted scalar, per the foundations
#: finding that basal expression carries no usable directional information
#: about the response template.
MAIN_EFFECT_ESTIMATORS: dict[str, Callable[..., np.ndarray]] = {
    "M0_zero": zero_main_effect,
    "M1_source_pooled": cme.source_pooled_mean,
    "M2_basal_weighted": cme.basal_weighted_mean,
    "M3a_source_pooled_shrunk": source_pooled_shrunk_mean,
    "M3b_basal_shrunk": cme.basal_shrunk_mean,
}


# --------------------------------------------------------------------------
# beta_hat: perturbation-specific effect from direct measurements only
# --------------------------------------------------------------------------


def centred_source_beta(delta: np.ndarray, sources: Sequence[int]) -> np.ndarray:
    """Perturbation-specific effect shared by ``sources``, centred over perturbations.

    Each source's own main effect is removed *before* averaging, so a source
    with an unusually large template cannot leak that template into
    ``beta_hat``.
    """
    d = np.asarray(delta, dtype=np.float64)[list(sources)]
    centred = d - d.mean(axis=1, keepdims=True)
    return centred.mean(axis=0)


def source_subsets(sources: Sequence[int], size: int) -> list[tuple[int, ...]]:
    """Every subset of ``sources`` of the given size, in frozen order.

    Tier 2 on Arc means *two* public contexts, not three, so the public
    benchmark must estimate ``beta_hat`` from two of its three available
    sources to be a like-for-like analogue. Averaging over all such subsets
    removes the arbitrariness of picking one.
    """
    if size < 1 or size > len(sources):
        raise ValueError(f"cannot take subsets of size {size} from {len(sources)} sources")
    return [tuple(c) for c in combinations(sorted(sources), size)]


def tier_beta(delta: np.ndarray, sources: Sequence[int], tier: int) -> np.ndarray:
    """``beta_hat`` for one tier, given the source contexts it may read.

    Tier 0 returns exact zeros regardless of what ``sources`` contains: the
    frozen policy is that a perturbation with no direct evidence gets no
    perturbation-specific prediction, and that is enforced here rather than
    left to the caller.
    """
    d = np.asarray(delta, dtype=np.float64)
    if tier == 0:
        return np.zeros((d.shape[1], d.shape[2]), dtype=np.float64)
    if tier not in TIER_SOURCE_COUNT:
        raise ValueError(f"unknown tier {tier!r}; known: {sorted(TIER_SOURCE_COUNT)}")
    return centred_source_beta(d, sources)


# --------------------------------------------------------------------------
# the model
# --------------------------------------------------------------------------


def predict(main_effect: np.ndarray, beta: np.ndarray, weight: float) -> np.ndarray:
    """``delta_hat[p] = m_hat + w beta_hat[p]``."""
    m = np.asarray(main_effect, dtype=np.float64)
    b = np.asarray(beta, dtype=np.float64)
    if b.ndim != 2 or b.shape[1] != m.shape[-1]:
        raise ValueError(f"beta {b.shape} is not (P, {m.shape[-1]})")
    return m[None, :] + float(weight) * b


@dataclass(frozen=True)
class TierWeights:
    """The selected shrinkage per tier, and the grid search that produced it."""

    weights: dict[int, float]
    #: ``(tier, weight) -> mean squared error`` over the nested inner folds.
    grid: dict[tuple[int, float], float]
    #: Inner folds contributing to the selection, for reporting.
    n_inner_folds: int

    def __getitem__(self, tier: int) -> float:
        return self.weights[tier]


def select_tier_shrinkage(
    delta: np.ndarray,
    sources: Sequence[int],
    *,
    control_means: np.ndarray,
    main_effect: Callable[..., np.ndarray],
    tiers: Sequence[int] = (2, 1),
    grid: Sequence[float] = SHRINKAGE_GRID,
) -> TierWeights:
    """Choose ``w`` per tier by nested validation **inside the source contexts**.

    One inner fold per source: that source is held out, the remaining sources
    play the role of the public evidence, and the tier's source-count structure
    is preserved (Tier 2 draws ``beta_hat`` from two inner sources, Tier 1 from
    one). The outer target context is never read, at any point, by anything in
    this function — which is what makes the selected weights usable on Arc.

    The objective is mean squared error. Correlation is not usable here:
    ``w`` rescales ``beta_hat`` only, so it cannot change the Pearson or
    Spearman correlation of the perturbation-specific part at all, and an
    MSE-free criterion would select arbitrarily.
    """
    d = np.asarray(delta, dtype=np.float64)
    src = list(sources)
    if len(src) < 2:
        raise ValueError("nested shrinkage selection needs at least two source contexts")

    totals: dict[tuple[int, float], float] = {(t, float(w)): 0.0 for t in tiers for w in grid}
    n_folds = 0
    for held in src:
        inner_sources = [c for c in src if c != held]
        truth = d[held]
        m_hat = main_effect(d, inner_sources, control_means=control_means, target=held)
        n_folds += 1
        for tier in tiers:
            size = TIER_SOURCE_COUNT[tier]
            if size is None:
                continue
            if size > len(inner_sources):
                # The inner loop cannot reproduce this tier's source count;
                # fall back to every inner source, and say so via n_inner_folds.
                subsets = [tuple(inner_sources)]
            else:
                subsets = source_subsets(inner_sources, size)
            for weight in grid:
                err = 0.0
                for subset in subsets:
                    pred = predict(m_hat, centred_source_beta(d, subset), weight)
                    diff = truth - pred
                    err += float(np.mean(diff * diff))
                totals[(tier, float(weight))] += err / len(subsets)

    weights: dict[int, float] = {}
    for tier in tiers:
        if TIER_SOURCE_COUNT[tier] is None:
            weights[tier] = 0.0
            continue
        weights[tier] = min(grid, key=lambda w, t=tier: totals[(t, float(w))])
    weights[0] = 0.0
    return TierWeights(
        weights=weights,
        grid={k: v / max(n_folds, 1) for k, v in totals.items()},
        n_inner_folds=n_folds,
    )
