r"""Reliability-aware zero-shot transferability / confidence estimation.

The point predictor is **frozen**: ``B = scale-calibrated source-only conserved
transfer``, exactly as established and used in the pathway phases. Nothing here
changes ``B``. This module only predicts *how trustworthy* ``B`` will be for a
given (context, perturbation), using information available before the target
perturbation is ever measured.

Reliability-aware targets
-------------------------
With two independent target half-responses ``h1 = L + e1``, ``h2 = L + e2``
(``e_i`` zero-mean, independent of ``L`` and of each other) and a prediction
``B`` that is independent of the target noise — which source-only predictions
are by construction:

.. math::
    E\langle h_1, h_2\rangle = \|L\|^2 \qquad
    E\langle h_1 - B, h_2 - B\rangle = \|L - B\|^2

because every noise cross-term has zero expectation. The first is
**signal energy**, the second **reproducible residual energy** — an unbiased
estimate of the true squared prediction error, *not* inflated by measurement
noise the way ``\|h - B\|^2`` is.

.. math::
    D = \frac{\langle h_1 - B, h_2 - B\rangle}{\langle h_1, h_2\rangle}

``D ~ 0`` transfer explains most reproducible response energy; ``D ~ 1`` it is no
better than predicting zero; ``D > 1`` it is **worse** than predicting zero.
``D`` is never clipped — the ``>1`` regime is real and was observed throughout
this project.

The denominator is itself a noisy estimate of ``\|L\|^2`` and becomes unstable as
it approaches zero. ``D`` is therefore a **secondary** target, evaluated only on
perturbations passing a stability rule that is derived from synthetic simulation
*before* any outer fold is touched (:func:`derive_stability_rule`). The
**primary** target is the reproducible residual energy itself, which has no
denominator and no instability.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

EPS = 1e-12

#: Predeclared ridge penalties for the confidence model.
CONFIDENCE_ALPHAS: tuple[float, ...] = (0.1, 1.0, 10.0, 100.0, 1000.0)

#: Predeclared coverage levels for risk-coverage analysis.
COVERAGE_LEVELS: tuple[float, ...] = (1.00, 0.75, 0.50, 0.25, 0.10)

#: Predeclared experiment budgets for the prioritisation analysis.
BUDGETS: tuple[float, ...] = (0.05, 0.10, 0.20, 0.30, 0.50)

#: Frozen feature set. Every entry is computable from source contexts plus
#: target *basal* expression, before the target perturbation is measured.
FEATURE_NAMES: tuple[str, ...] = (
    "source_agreement",
    "source_min_agreement",
    "source_sd_ratio",
    "source_sign_consistency",
    "source_reliability",
    "source_cells_mean",
    "source_cells_min",
    "source_magnitude",
    "source_reliable_energy",
    "source_magnitude_spread",
    "target_gene_basal",
    "basal_similarity_mean",
)

#: Predeclared ablation groups (G).
ABLATIONS: dict[str, tuple[str, ...]] = {
    "agreement_only": ("source_agreement",),
    "agreement_plus_quality": (
        "source_agreement",
        "source_reliability",
        "source_cells_mean",
        "source_cells_min",
    ),
    "agreement_plus_magnitude": (
        "source_agreement",
        "source_magnitude",
        "source_reliable_energy",
        "source_magnitude_spread",
    ),
    "agreement_plus_basal": (
        "source_agreement",
        "target_gene_basal",
        "basal_similarity_mean",
    ),
    "full": FEATURE_NAMES,
}


# --------------------------------------------------------------------------
# reliability-aware quantities
# --------------------------------------------------------------------------


def signal_energy(h1: np.ndarray, h2: np.ndarray) -> np.ndarray:
    """``<h1, h2>`` — unbiased estimate of ``||L||^2`` per perturbation."""
    a = np.asarray(h1, dtype=np.float64)
    b = np.asarray(h2, dtype=np.float64)
    if a.shape != b.shape:
        raise ValueError(f"shape mismatch {a.shape} vs {b.shape}")
    return np.einsum("...g,...g->...", a, b)


def residual_energy(h1: np.ndarray, h2: np.ndarray, prediction: np.ndarray) -> np.ndarray:
    """``<h1 - B, h2 - B>`` — unbiased estimate of ``||L - B||^2``."""
    a = np.asarray(h1, dtype=np.float64)
    b = np.asarray(h2, dtype=np.float64)
    p = np.asarray(prediction, dtype=np.float64)
    if a.shape != b.shape or a.shape != p.shape:
        raise ValueError("shape mismatch between halves and prediction")
    return np.einsum("...g,...g->...", a - p, b - p)


def normalised_d(
    h1: np.ndarray, h2: np.ndarray, prediction: np.ndarray, *, min_signal: float
) -> tuple[np.ndarray, np.ndarray]:
    """``(D, stable)`` with ``D`` masked to NaN where the denominator is unusable.

    ``min_signal`` must come from :func:`derive_stability_rule`, fixed before any
    outer fold is evaluated. ``D`` is **not** clipped.
    """
    sig = signal_energy(h1, h2)
    res = residual_energy(h1, h2, prediction)
    stable = sig > min_signal
    out = np.full(sig.shape, np.nan)
    out[stable] = res[stable] / sig[stable]
    return out, stable


def derive_stability_rule(
    *,
    n_genes: int,
    noise_sd: float,
    signal_levels: Sequence[float] = tuple(np.linspace(0.05, 2.0, 40)),
    n_rep: int = 400,
    tolerance: float = 0.25,
    coverage: float = 0.90,
    seed: int = 0,
) -> dict[str, float]:
    """Smallest signal energy at which ``D`` is usable, from synthetic simulation.

    For a grid of latent signal strengths, simulates independent halves with a
    known true ``D`` and measures how often the estimate lands within
    ``tolerance`` of the truth. Returns the smallest **estimated** signal energy
    at which at least ``coverage`` of replicates are within tolerance.

    This is run once, on synthetic data, **before** any outer fold is evaluated,
    and its output is then frozen.
    """
    rng = np.random.default_rng(seed)
    records = []
    for level in signal_levels:
        latent = rng.normal(scale=level, size=(n_rep, n_genes))
        pred = 0.5 * latent  # a known, fixed-quality predictor
        true_d = 0.25  # ||L - 0.5L||^2 / ||L||^2
        h1 = latent + rng.normal(scale=noise_sd, size=latent.shape)
        h2 = latent + rng.normal(scale=noise_sd, size=latent.shape)
        sig = signal_energy(h1, h2)
        res = residual_energy(h1, h2, pred)
        with np.errstate(divide="ignore", invalid="ignore"):
            d_hat = np.where(np.abs(sig) > EPS, res / sig, np.nan)
        ok = np.abs(d_hat - true_d) <= tolerance
        records.append(
            {
                "level": float(level),
                "median_signal": float(np.median(sig)),
                "within_tolerance": float(np.nanmean(ok)),
            }
        )
    passing = [r for r in records if r["within_tolerance"] >= coverage]
    threshold = min((r["median_signal"] for r in passing), default=float("inf"))
    return {
        "min_signal_energy": float(threshold),
        "tolerance": float(tolerance),
        "coverage": float(coverage),
        "noise_sd": float(noise_sd),
        "n_genes": int(n_genes),
        "curve": records,
    }


# --------------------------------------------------------------------------
# frozen confidence features (source-only + target basal)
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class SourceStats:
    """Per-(context, perturbation) source-side quantities, computed once."""

    pair_corr: np.ndarray = field(repr=False)  # (C, C, P)
    reliability: np.ndarray = field(repr=False)  # (C, P)
    cells: np.ndarray = field(repr=False)  # (C, P)
    magnitude: np.ndarray = field(repr=False)  # (C, P)
    gene_basal: np.ndarray = field(repr=False)  # (C, P)
    basal_sim: np.ndarray = field(repr=False)  # (C, C)


def source_agreement(stats: SourceStats, sources: Sequence[int]) -> np.ndarray:
    """**Frozen definition**: mean pairwise gene-level Pearson among source responses.

    Unchanged from every previous phase, and not redefined after modelling.
    """
    src = list(sources)
    pairs = [(a, b) for i, a in enumerate(src) for b in src[i + 1 :]]
    if not pairs:
        raise ValueError("Source agreement needs at least two source contexts.")
    return np.nanmean(np.vstack([stats.pair_corr[a, b] for a, b in pairs]), axis=0)


def build_confidence_features(
    Y: np.ndarray,
    halves: np.ndarray,
    target: int,
    sources: Sequence[int],
    *,
    stats: SourceStats,
) -> np.ndarray:
    """``(P, F)`` design matrix. Never reads ``Y[target]`` or the target's halves."""
    src = list(sources)
    Ys = np.asarray(Y, dtype=np.float64)[src]
    A = Ys.mean(axis=0)
    pairs = [(a, b) for i, a in enumerate(src) for b in src[i + 1 :]]
    pair_stack = np.vstack([stats.pair_corr[a, b] for a, b in pairs])

    a_norm = np.linalg.norm(A, axis=1)
    sd_norm = np.linalg.norm(Ys.std(axis=0), axis=1)
    sign_consistency = np.abs(A).sum(axis=1) / (np.abs(Ys).mean(axis=0).sum(axis=1) + EPS)
    src_h1 = np.asarray(halves, dtype=np.float64)[:, 0][:, src].mean(axis=1).mean(axis=0)
    src_h2 = np.asarray(halves, dtype=np.float64)[:, 1][:, src].mean(axis=1).mean(axis=0)

    columns = {
        "source_agreement": np.nanmean(pair_stack, axis=0),
        "source_min_agreement": np.nanmin(pair_stack, axis=0),
        "source_sd_ratio": sd_norm / (a_norm + EPS),
        "source_sign_consistency": sign_consistency,
        "source_reliability": np.nanmean(stats.reliability[src], axis=0),
        "source_cells_mean": stats.cells[src].mean(axis=0),
        "source_cells_min": stats.cells[src].min(axis=0),
        "source_magnitude": a_norm,
        "source_reliable_energy": signal_energy(src_h1, src_h2),
        "source_magnitude_spread": stats.magnitude[src].std(axis=0),
        "target_gene_basal": np.nan_to_num(stats.gene_basal[target]),
        "basal_similarity_mean": np.full(
            Y.shape[1], float(np.mean([stats.basal_sim[target, s] for s in src]))
        ),
    }
    X = np.stack([columns[name] for name in FEATURE_NAMES], axis=1)
    return np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)


# --------------------------------------------------------------------------
# low-capacity confidence models
# --------------------------------------------------------------------------


class Standardiser:
    """Mean/SD standardisation fitted on training rows only."""

    def __init__(self) -> None:
        self.mean_: np.ndarray | None = None
        self.scale_: np.ndarray | None = None

    def fit(self, X: np.ndarray) -> Standardiser:
        Xa = np.asarray(X, dtype=np.float64)
        self.mean_ = Xa.mean(axis=0)
        sd = Xa.std(axis=0)
        self.scale_ = np.where(sd > EPS, sd, 1.0)
        return self

    def transform(self, X: np.ndarray) -> np.ndarray:
        if self.mean_ is None:
            raise RuntimeError("Standardiser used before fit().")
        return (np.asarray(X, dtype=np.float64) - self.mean_) / self.scale_


def fit_ridge(X: np.ndarray, y: np.ndarray, alpha: float) -> np.ndarray:
    Xa = np.asarray(X, dtype=np.float64)
    ya = np.asarray(y, dtype=np.float64).ravel()
    n, f = Xa.shape
    Xb = np.hstack([np.ones((n, 1)), Xa])
    pen = np.eye(f + 1) * alpha
    pen[0, 0] = 0.0
    return np.linalg.solve(Xb.T @ Xb + pen, Xb.T @ ya)


def predict_ridge(X: np.ndarray, coef: np.ndarray) -> np.ndarray:
    return coef[0] + np.asarray(X, dtype=np.float64) @ coef[1:]


def fit_monotone(x: np.ndarray, y: np.ndarray):
    """Isotonic (monotone, non-parametric, one-dimensional) calibration."""
    from sklearn.isotonic import IsotonicRegression

    ok = np.isfinite(x) & np.isfinite(y)
    iso = IsotonicRegression(increasing="auto", out_of_bounds="clip")
    iso.fit(np.asarray(x, dtype=np.float64)[ok], np.asarray(y, dtype=np.float64)[ok])
    return iso


# --------------------------------------------------------------------------
# selective prediction
# --------------------------------------------------------------------------


def risk_coverage(
    confidence: np.ndarray, risk: np.ndarray, coverages: Sequence[float] = COVERAGE_LEVELS
) -> list[dict[str, float]]:
    """Mean/median risk among the most-confident fraction, at each coverage.

    ``confidence`` is higher-is-better. Perturbations with non-finite confidence
    or risk are dropped before ranking.
    """
    c = np.asarray(confidence, dtype=np.float64)
    r = np.asarray(risk, dtype=np.float64)
    ok = np.isfinite(c) & np.isfinite(r)
    c, r = c[ok], r[ok]
    order = np.argsort(-c, kind="stable")
    r_sorted = r[order]
    out = []
    for cov in coverages:
        k = max(1, int(round(cov * r_sorted.size)))
        sel = r_sorted[:k]
        out.append(
            {
                "coverage": float(cov),
                "n": int(k),
                "mean_risk": float(np.mean(sel)),
                "median_risk": float(np.median(sel)),
            }
        )
    return out


def area_under_risk_coverage(confidence: np.ndarray, risk: np.ndarray) -> float:
    """Mean of the running-average risk over all coverage levels (lower is better)."""
    c = np.asarray(confidence, dtype=np.float64)
    r = np.asarray(risk, dtype=np.float64)
    ok = np.isfinite(c) & np.isfinite(r)
    c, r = c[ok], r[ok]
    if c.size == 0:
        return float("nan")
    r_sorted = r[np.argsort(-c, kind="stable")]
    running = np.cumsum(r_sorted) / np.arange(1, r_sorted.size + 1)
    return float(np.mean(running))


def prioritisation_capture(
    confidence: np.ndarray, risk: np.ndarray, budgets: Sequence[float] = BUDGETS
) -> list[dict[str, float]]:
    """Fraction of total reproducible error captured by the least-confident budget.

    Risk is clipped at zero *for the capture accounting only*, because a
    reproducible-error estimate can be slightly negative from noise and negative
    contributions would make "fraction of total error" meaningless. The
    underlying risk values are never modified.
    """
    c = np.asarray(confidence, dtype=np.float64)
    r = np.asarray(risk, dtype=np.float64)
    ok = np.isfinite(c) & np.isfinite(r)
    c, r = c[ok], np.clip(r[ok], 0.0, None)
    total = float(r.sum())
    order = np.argsort(c, kind="stable")  # ascending confidence = worst first
    r_worst = r[order]
    out = []
    for b in budgets:
        k = max(1, int(round(b * r.size)))
        captured = float(r_worst[:k].sum())
        out.append(
            {
                "budget": float(b),
                "n": int(k),
                "captured_fraction": captured / total if total > EPS else np.nan,
                "random_expectation": float(b),
            }
        )
    return out


def calibration_curve(
    confidence: np.ndarray, quality: np.ndarray, n_bins: int = 5
) -> list[dict[str, float]]:
    """Observed quality within equal-count bins of predicted confidence."""
    c = np.asarray(confidence, dtype=np.float64)
    q = np.asarray(quality, dtype=np.float64)
    ok = np.isfinite(c) & np.isfinite(q)
    c, q = c[ok], q[ok]
    if c.size < n_bins:
        return []
    edges = np.quantile(c, np.linspace(0, 1, n_bins + 1))
    out = []
    for i in range(n_bins):
        lo, hi = edges[i], edges[i + 1]
        m = (c >= lo) & (c <= hi) if i == n_bins - 1 else (c >= lo) & (c < hi)
        if m.sum() == 0:
            continue
        out.append(
            {
                "bin": i,
                "n": int(m.sum()),
                "mean_confidence": float(np.mean(c[m])),
                "mean_quality": float(np.mean(q[m])),
                "median_quality": float(np.median(q[m])),
            }
        )
    return out
