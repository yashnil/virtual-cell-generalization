r"""Leave-one-context-out zero-shot transfer: splits, baselines, metrics, reliability.

Diagnostic machinery for the zero-shot recoverability phase. **No learned model
lives here** — every baseline is a closed-form function of source-context
responses and basal (control) expression.

Leakage rules
-------------
For held-out context ``c*`` the following are ALLOWED as prediction inputs:

* basal / control cells of ``c*`` (its control mean profile),
* perturbation identity,
* perturbation responses of the three source contexts,
* external priors that do not encode target perturbation outcomes.

The following are FORBIDDEN as prediction inputs:

* any perturbation-response measurement from ``c*``,
* feature selection, normalisation or hyperparameter choices fitted to ``c*``
  perturbation responses,
* the four-context ``beta`` or ``gamma``, because both are computed from all
  four contexts and therefore contain the held-out response.

Target perturbation responses enter **only** at evaluation time. Functions that
consume them are named ``*_evaluation_only`` or documented as such.

Key algebra
-----------
Write the frozen four-context decomposition as
``delta[c,p] = mu + alpha_c + beta_p + gamma[c,p]`` with the balanced-design
side conditions ``sum_c alpha_c = 0`` and ``sum_c gamma[c,p] = 0``.

With three sources ``S`` and held-out ``c*``, averaging over ``S`` gives
``mean_S alpha = -alpha_{c*}/3`` and ``mean_S gamma[.,p] = -gamma[c*,p]/3``, so
the source-mean prediction is

.. math::
    A[p] = \operatorname{mean}_{c \in S} delta[c,p]
         = mu - \frac{alpha_{c*}}{3} + beta_p - \frac{gamma[c*,p]}{3}

Two consequences matter and are asserted in the tests:

1. **The residual of source-mean transfer is exactly**
   ``delta[c*,p] - A[p] = (4/3) (alpha_{c*} + gamma[c*,p])``.
   Centring it over perturbations removes ``alpha_{c*}`` and leaves exactly
   ``(4/3) gamma[c*,p]``. So "predict the centred transfer residual" and
   "predict gamma" are the *same problem* up to a constant factor.
2. **The source mean carries a negative gamma component** (``-gamma[c*,p]/3``).
   Source-mean transfer is therefore not gamma-neutral: it is mildly
   anti-correlated with the held-out interaction by construction.

Reliability correction
----------------------
Additive independent-noise model: an observed half response is
``h = L + e`` with latent ``L`` and noise ``e`` independent of ``L`` and
between halves, all treated as vectors over genes.

* ``corr(h1, h2) = Var(L) / (Var(L) + Var(e)) =: rho_half`` — the observed
  split-half correlation *is* the reliability of a single half.
* For a perfect latent predictor ``P = L``,
  ``corr(P, h1) = Var(L) / sqrt(Var(L) (Var(L) + Var(e))) = sqrt(rho_half)``.
  **The ceiling is the square root of reliability, not reliability itself.**
* Hence for any predictor, ``corr(P, h) = rho_PL * sqrt(rho_half)``, so the
  disattenuated estimate of the latent correlation is
  ``rho_PL = corr(P, h) / sqrt(rho_half)``.
* Averaging the two independent halves halves the noise variance, so the
  reliability of the full (both-half) estimate follows Spearman-Brown:
  ``rho_full = 2 rho_half / (1 + rho_half)``.

Assumptions, all of which fail quietly if ignored: noise independent of the
latent signal, independent between halves, and ``rho_half > 0``. The estimator
is validated against synthetic data with a known latent signal in
``tests/test_loco.py``.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

EPS = 1e-12


# --------------------------------------------------------------------------
# folds
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class LocoFold:
    """One leave-one-context-out evaluation."""

    target: str
    sources: tuple[str, ...]
    target_index: int
    source_indices: tuple[int, ...]

    @property
    def name(self) -> str:
        return f"{'+'.join(self.sources)}__to__{self.target}"


def make_folds(contexts: Sequence[str]) -> list[LocoFold]:
    """One fold per context, each holding that context out. Order is frozen."""
    if len(contexts) < 2:
        raise ValueError("Need at least two contexts for leave-one-context-out.")
    folds = []
    for i, target in enumerate(contexts):
        sources = tuple(c for j, c in enumerate(contexts) if j != i)
        folds.append(
            LocoFold(
                target=target,
                sources=sources,
                target_index=i,
                source_indices=tuple(j for j in range(len(contexts)) if j != i),
            )
        )
    return folds


# --------------------------------------------------------------------------
# basal (control) similarity — the only context feature allowed
# --------------------------------------------------------------------------


def basal_similarity(control_means: np.ndarray) -> np.ndarray:
    """Pearson correlation between context basal profiles. Controls only."""
    return np.corrcoef(np.asarray(control_means, dtype=np.float64))


def basal_affine_weights(
    control_means: np.ndarray, target: int, sources: Sequence[int]
) -> np.ndarray:
    r"""Weights expressing the target's basal profile as an affine mix of sources.

    Solves ``min_w ||b_target - B w||^2`` subject to ``sum_i w_i = 1`` in closed
    form. Parameter-free: nothing is tuned, and only control profiles are used,
    so no target perturbation response can enter.
    """
    B = np.asarray(control_means, dtype=np.float64)[list(sources)].T  # (G, S)
    b = np.asarray(control_means, dtype=np.float64)[target]
    n = B.shape[1]
    G = B.T @ B
    ones = np.ones(n)
    # KKT system for the equality-constrained least squares
    kkt = np.zeros((n + 1, n + 1))
    kkt[:n, :n] = 2 * G
    kkt[:n, n] = ones
    kkt[n, :n] = ones
    rhs = np.zeros(n + 1)
    rhs[:n] = 2 * (B.T @ b)
    rhs[n] = 1.0
    return np.linalg.solve(kkt, rhs)[:n]


def basal_simplex_weights(
    control_means: np.ndarray, target: int, sources: Sequence[int]
) -> np.ndarray:
    """As :func:`basal_affine_weights` but constrained to be non-negative.

    Uses non-negative least squares on the same basal-only objective, then
    renormalises to sum to one. Still parameter-free.
    """
    from scipy.optimize import nnls

    B = np.asarray(control_means, dtype=np.float64)[list(sources)].T
    b = np.asarray(control_means, dtype=np.float64)[target]
    w, _ = nnls(B, b)
    total = w.sum()
    if total <= EPS:
        return np.full(len(sources), 1.0 / len(sources))
    return w / total


# --------------------------------------------------------------------------
# zero-shot baselines (source responses + basal only)
# --------------------------------------------------------------------------


def source_mean(delta: np.ndarray, sources: Sequence[int]) -> np.ndarray:
    """Baseline A — conserved transfer: the mean source response per perturbation."""
    return np.asarray(delta, dtype=np.float64)[list(sources)].mean(axis=0)


def single_source(delta: np.ndarray, source: int) -> np.ndarray:
    """Baseline B — predict the target with one source context's response."""
    return np.asarray(delta, dtype=np.float64)[source]


def nearest_context(
    control_means: np.ndarray, target: int, sources: Sequence[int]
) -> tuple[int, np.ndarray]:
    """Baseline C support: the basal-nearest source and the similarity row."""
    sim = basal_similarity(control_means)[target]
    best = max(sources, key=lambda s: sim[s])
    return best, sim


def weighted_source(delta: np.ndarray, sources: Sequence[int], weights: np.ndarray) -> np.ndarray:
    """Baseline D — basal-similarity-weighted combination of source responses."""
    w = np.asarray(weights, dtype=np.float64)
    if w.shape != (len(sources),):
        raise ValueError(f"weights shape {w.shape} != ({len(sources)},)")
    return np.tensordot(w, np.asarray(delta, dtype=np.float64)[list(sources)], axes=(0, 0))


# --------------------------------------------------------------------------
# metrics
# --------------------------------------------------------------------------


def _safe_pearson(a: np.ndarray, b: np.ndarray) -> float:
    if a.std() < EPS or b.std() < EPS:
        return np.nan
    r = float(np.corrcoef(a, b)[0, 1])
    return r if np.isfinite(r) else np.nan


def per_perturbation_metrics(
    observed: np.ndarray, predicted: np.ndarray, *, with_spearman: bool = True
) -> pd.DataFrame:
    """Gene-wise metrics for every perturbation.

    ``observed`` and ``predicted`` are ``(P, G)``. Metrics are computed across
    genes within each perturbation.
    """
    observed = np.asarray(observed, dtype=np.float64)
    predicted = np.asarray(predicted, dtype=np.float64)
    if observed.shape != predicted.shape:
        raise ValueError(f"shape mismatch {observed.shape} vs {predicted.shape}")
    rows = []
    for i in range(observed.shape[0]):
        y, x = observed[i], predicted[i]
        err = y - x
        denom = float(y @ y)
        rows.append(
            {
                "pearson": _safe_pearson(y, x),
                "spearman": float(stats.spearmanr(y, x).statistic) if with_spearman else np.nan,
                "cosine": float(y @ x / (np.linalg.norm(y) * np.linalg.norm(x) + EPS)),
                "mse": float(err @ err / len(y)),
                "energy_explained": float(1.0 - (err @ err) / denom) if denom > EPS else np.nan,
                "obs_norm": float(np.linalg.norm(y)),
                "pred_norm": float(np.linalg.norm(x)),
            }
        )
    return pd.DataFrame(rows)


def template_removed(matrix: np.ndarray, template: np.ndarray) -> np.ndarray:
    """Subtract a shared template row from every perturbation."""
    return np.asarray(matrix, dtype=np.float64) - np.asarray(template, dtype=np.float64)[None, :]


def target_template_evaluation_only(observed: np.ndarray) -> np.ndarray:
    """Mean observed target response over perturbations.

    **Evaluation-only.** This is computed from held-out responses and must never
    be fed to a predictor; it exists so that template-removed metrics can be
    reported the way the cross-context literature does.
    """
    return np.asarray(observed, dtype=np.float64).mean(axis=0)


# --------------------------------------------------------------------------
# reliability
# --------------------------------------------------------------------------


def split_half_reliability(half_a: np.ndarray, half_b: np.ndarray) -> np.ndarray:
    """Per-perturbation ``corr(h1, h2)`` across genes — reliability of one half."""
    a = np.asarray(half_a, dtype=np.float64)
    b = np.asarray(half_b, dtype=np.float64)
    if a.shape != b.shape:
        raise ValueError(f"shape mismatch {a.shape} vs {b.shape}")
    return np.array([_safe_pearson(a[i], b[i]) for i in range(a.shape[0])])


def spearman_brown(rho_half: np.ndarray | float) -> np.ndarray | float:
    """Reliability of the pooled estimate from the reliability of one half."""
    r = np.asarray(rho_half, dtype=np.float64)
    return 2.0 * r / (1.0 + r)


def disattenuate(
    r_observed: np.ndarray | float,
    rho: np.ndarray | float,
    *,
    min_rho: float = 0.05,
    clip: tuple[float, float] | None = (-1.5, 1.5),
) -> np.ndarray:
    r"""Correct an observed correlation for measurement noise in the target.

    ``rho_latent = r_observed / sqrt(rho)`` where ``rho`` is the reliability of
    the *measurement the correlation was computed against* (see the module
    docstring for the derivation).

    Returns NaN where ``rho <= min_rho``: at low reliability the denominator is
    both tiny and badly estimated, so the correction amplifies noise without
    bound and the value is not interpretable. Results are clipped to ``clip``
    because sampling error can push a corrected value past 1.
    """
    r = np.asarray(r_observed, dtype=np.float64)
    rho_arr = np.asarray(rho, dtype=np.float64)
    out = np.full(np.broadcast(r, rho_arr).shape, np.nan, dtype=np.float64)
    usable = np.isfinite(r) & np.isfinite(rho_arr) & (rho_arr > min_rho)
    out[usable] = (r * np.ones_like(out))[usable] / np.sqrt((rho_arr * np.ones_like(out))[usable])
    if clip is not None:
        out = np.clip(out, *clip)
    return out


def ceiling_from_reliability(rho: np.ndarray | float) -> np.ndarray:
    """Maximum attainable correlation against a measurement of reliability ``rho``."""
    r = np.asarray(rho, dtype=np.float64)
    return np.sqrt(np.clip(r, 0.0, None))


def reliability_normalised_evaluation(
    predicted: np.ndarray,
    half_a: np.ndarray,
    half_b: np.ndarray,
    *,
    min_rho: float = 0.05,
) -> pd.DataFrame:
    """Evaluate a prediction against two independent target halves.

    Reports the raw correlation against each half, the split-half reliability of
    the target, the implied ceiling ``sqrt(rho_half)``, and the disattenuated
    correlation. Nothing here is fitted; the halves are the only extra input.
    """
    r_a = np.array(
        [
            _safe_pearson(np.asarray(half_a)[i], np.asarray(predicted)[i])
            for i in range(np.asarray(predicted).shape[0])
        ]
    )
    r_b = np.array(
        [
            _safe_pearson(np.asarray(half_b)[i], np.asarray(predicted)[i])
            for i in range(np.asarray(predicted).shape[0])
        ]
    )
    rho_half = split_half_reliability(half_a, half_b)
    r_mean = np.nanmean(np.vstack([r_a, r_b]), axis=0)
    return pd.DataFrame(
        {
            "r_half_a": r_a,
            "r_half_b": r_b,
            "r_half_mean": r_mean,
            "rho_half": rho_half,
            "rho_full_spearman_brown": spearman_brown(rho_half),
            "ceiling": ceiling_from_reliability(rho_half),
            "r_disattenuated": disattenuate(r_mean, rho_half, min_rho=min_rho),
        }
    )


# --------------------------------------------------------------------------
# bootstrap
# --------------------------------------------------------------------------


def bootstrap_median(
    values: np.ndarray, *, n_boot: int = 2000, seed: int = 0, alpha: float = 0.05
) -> dict[str, float]:
    """Median with a percentile bootstrap CI over perturbations."""
    v = np.asarray(values, dtype=np.float64)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return {"median": np.nan, "lo": np.nan, "hi": np.nan, "n": 0}
    rng = np.random.default_rng(seed)
    boot = np.median(rng.choice(v, size=(n_boot, v.size), replace=True), axis=1)
    return {
        "median": float(np.median(v)),
        "lo": float(np.percentile(boot, 100 * alpha / 2)),
        "hi": float(np.percentile(boot, 100 * (1 - alpha / 2))),
        "n": int(v.size),
    }
