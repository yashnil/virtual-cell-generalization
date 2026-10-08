"""N4: source agreement vs simple statistics and an exchangeable random-effects null.

Preregistered in ``reports/n1_n4_protocol.md`` §3. Quality targets use per-repeat
unbiased energies (independent halves) averaged over repeats.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from scipy import stats

EPS = 1e-12


def rowwise_pearson(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    ac = a - a.mean(axis=-1, keepdims=True)
    bc = b - b.mean(axis=-1, keepdims=True)
    den = np.linalg.norm(ac, axis=-1) * np.linalg.norm(bc, axis=-1)
    return np.where(den > EPS, np.einsum("...g,...g->...", ac, bc) / np.maximum(den, EPS), np.nan)


def agreement(S: np.ndarray) -> np.ndarray:
    """Frozen definition: mean pairwise gene-level Pearson among sources ``(S, P, G)``."""
    n = S.shape[0]
    pairs = [rowwise_pearson(S[a], S[b]) for a in range(n) for b in range(a + 1, n)]
    return np.nanmean(np.vstack(pairs), axis=0)


def quality_d(
    h1: Sequence[np.ndarray], h2: Sequence[np.ndarray], B: np.ndarray, *, min_signal: float
) -> tuple[np.ndarray, np.ndarray]:
    """``(−D, signal)`` from per-repeat independent halves, energies averaged over repeats."""
    sig = np.mean([np.einsum("pg,pg->p", a, b) for a, b in zip(h1, h2, strict=True)], axis=0)
    res = np.mean(
        [np.einsum("pg,pg->p", a - B, b - B) for a, b in zip(h1, h2, strict=True)], axis=0
    )
    q = np.full(sig.shape, np.nan)
    ok = sig > min_signal
    q[ok] = -res[ok] / sig[ok]
    return q, sig


def noise_variance(h1: Sequence[np.ndarray], h2: Sequence[np.ndarray]) -> np.ndarray:
    """Per-perturbation per-gene noise variance of a full-depth response.

    Two independent half-depth estimates differ by noise of variance 4 σ²_full.
    """
    G = h1[0].shape[1]
    return np.mean([np.sum((a - b) ** 2, axis=1) for a, b in zip(h1, h2, strict=True)], axis=0) / (
        4 * G
    )


def fit_exchangeable(S_c: np.ndarray, sigma2_src: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Plug-in ``m_p`` and ``τ²_p`` from centred sources ``(S, P, G)`` (protocol §3.4 T3)."""
    n_s, _, G = S_c.shape
    m = S_c.mean(axis=0)
    between = np.sum((S_c - m) ** 2, axis=(0, 2)) / ((n_s - 1) * G)
    tau2 = np.maximum(0.0, between - sigma2_src.mean(axis=0))
    return m, tau2


def simulate_null(
    m: np.ndarray,
    tau2: np.ndarray,
    sigma2_src: np.ndarray,
    sigma2_tgt: np.ndarray,
    *,
    scale: float,
    min_signal: float,
    rng: np.random.Generator,
) -> float:
    """One exchangeable-world replicate: Spearman(agreement', −D')."""
    n_s = sigma2_src.shape[0]
    P, G = m.shape
    sd_tau = np.sqrt(tau2)[:, None]
    S = np.empty((n_s, P, G), dtype=np.float32)
    for s in range(n_s):
        S[s] = m + sd_tau * rng.standard_normal((P, G), dtype=np.float32)
        S[s] += np.sqrt(sigma2_src[s])[:, None] * rng.standard_normal((P, G), dtype=np.float32)
    latent_t = m + sd_tau * rng.standard_normal((P, G), dtype=np.float32)
    sd_half = np.sqrt(2 * sigma2_tgt)[:, None]
    h1 = latent_t + sd_half * rng.standard_normal((P, G), dtype=np.float32)
    h2 = latent_t + sd_half * rng.standard_normal((P, G), dtype=np.float32)
    a = agreement(S.astype(np.float64))
    B = scale * S.mean(axis=0, dtype=np.float64)
    q, _ = quality_d([h1.astype(np.float64)], [h2.astype(np.float64)], B, min_signal=min_signal)
    ok = np.isfinite(a) & np.isfinite(q)
    return float(stats.spearmanr(a[ok], q[ok]).statistic)


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    ok = np.isfinite(x) & np.isfinite(y)
    return float(stats.spearmanr(x[ok], y[ok]).statistic)


def _rank(x: np.ndarray) -> np.ndarray:
    return stats.rankdata(x)


def paired_bootstrap_spearman_diff(
    score: np.ndarray, other: np.ndarray, q: np.ndarray, *, n_boot: int, rng: np.random.Generator
) -> tuple[float, float, float]:
    """``ρ(score, q) − |ρ(other, q)|`` and its 95 % paired-bootstrap interval."""
    ok = np.isfinite(score) & np.isfinite(other) & np.isfinite(q)
    s, o, y = score[ok], other[ok], q[ok]
    point = spearman(s, y) - abs(spearman(o, y))
    n = len(y)
    vals = np.empty(n_boot)
    for b in range(n_boot):
        i = rng.integers(0, n, n)
        rs, ro, ry = _rank(s[i]), _rank(o[i]), _rank(y[i])
        vals[b] = np.corrcoef(rs, ry)[0, 1] - abs(np.corrcoef(ro, ry)[0, 1])
    return point, float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))
