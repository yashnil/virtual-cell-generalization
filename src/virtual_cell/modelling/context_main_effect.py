"""Estimating a context's mean perturbation response without seeing it.

The VCC 2026 score anchors zero at the context's mean perturbation response,
``m_c = mean_p delta[c, p]``. In the decomposition
``delta[c,p] = mu + alpha_c + beta_p + gamma[c,p]`` the centring constraints
kill the last two terms under an average over ``p``, so ``m_c = mu + alpha_c``:
the context main effect, and nothing else.

That makes ``m_c`` worth estimating on its own, but it is **hidden** for Arc's
A, B and C -- it is a summary of exactly the perturbation outcomes the
challenge withholds. The distinction this module enforces:

``oracle_main_effect``
    ``m_c`` computed from the context's own responses. An evaluation target.
    Using it as a model input, at any remove, invalidates the result.

``FEASIBLE_ESTIMATORS``
    Functions of the target context's **control cells** and the **source
    contexts' responses** only. These are what could actually be applied to
    Arc.

A widespread shortcut is to assert that having the target controls makes the
score-0 baseline automatic. It does not: controls give ``alpha`` in *basal*
space, and ``m_c`` lives in *response* space. Whether one predicts the other is
an empirical question, answered on public leave-one-context-out folds.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

__all__ = [
    "oracle_main_effect",
    "source_pooled_mean",
    "basal_weighted_mean",
    "basal_shrunk_mean",
    "FEASIBLE_ESTIMATORS",
]


def oracle_main_effect(delta: np.ndarray, context: int) -> np.ndarray:
    """``m_c`` from the context's own responses. EVALUATION TARGET ONLY."""
    return np.asarray(delta)[context].mean(axis=0)


def _source_main_effects(delta: np.ndarray, sources: Sequence[int]) -> np.ndarray:
    return np.stack([np.asarray(delta)[c].mean(axis=0) for c in sources])


def source_pooled_mean(
    delta: np.ndarray,
    sources: Sequence[int],
    *,
    control_means: np.ndarray | None = None,
    target: int | None = None,
) -> np.ndarray:
    """``E0``: the unweighted mean of the source contexts' main effects.

    Reads no target information whatsoever, so it is the floor every other
    feasible estimator must beat to justify its extra machinery.
    """
    del control_means, target
    return _source_main_effects(delta, sources).mean(axis=0)


def _basal_weights(control_means: np.ndarray, target: int, sources: Sequence[int]) -> np.ndarray:
    basal = np.asarray(control_means, dtype=np.float64)
    centred = basal - basal.mean(axis=1, keepdims=True)
    norm = np.linalg.norm(centred, axis=1)
    unit = centred / np.where(norm > 0, norm, 1.0)[:, None]
    sim = unit[list(sources)] @ unit[target]
    w = np.clip(sim, 0.0, None)
    return w / w.sum() if w.sum() > 0 else np.full(len(sources), 1.0 / len(sources))


def basal_weighted_mean(
    delta: np.ndarray,
    sources: Sequence[int],
    *,
    control_means: np.ndarray,
    target: int,
) -> np.ndarray:
    """``E1``: source main effects weighted by basal similarity to the target.

    Uses the target's control cells only, which is a permitted Arc input.
    """
    weights = _basal_weights(control_means, target, sources)
    return weights @ _source_main_effects(delta, sources)


def basal_shrunk_mean(
    delta: np.ndarray,
    sources: Sequence[int],
    *,
    control_means: np.ndarray,
    target: int,
) -> np.ndarray:
    """``E2``: ``E1`` rescaled by a single scalar fitted on the sources.

    The scalar is fitted by leave-one-source-out: for each source, predict its
    main effect from the others and take the least-squares scale that would
    have been right, then use the median. One free parameter, fitted without
    touching the target -- deliberately at the bottom of the capacity ladder,
    because the foundations phase found no evidence for a richer basal-to-alpha
    map.
    """
    scales: list[float] = []
    for held in sources:
        rest = [c for c in sources if c != held]
        if not rest:
            continue
        pred = basal_weighted_mean(delta, rest, control_means=control_means, target=held)
        truth = oracle_main_effect(delta, held)
        denom = float(pred @ pred)
        if denom > 0:
            scales.append(float(pred @ truth) / denom)
    scale = float(np.median(scales)) if scales else 1.0
    return scale * basal_weighted_mean(delta, sources, control_means=control_means, target=target)


#: Estimators that read only the target's controls and the sources' responses.
FEASIBLE_ESTIMATORS = {
    "E0_source_pooled": source_pooled_mean,
    "E1_basal_weighted": basal_weighted_mean,
    "E2_basal_shrunk": basal_shrunk_mean,
}
