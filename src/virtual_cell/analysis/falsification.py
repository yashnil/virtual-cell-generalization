r"""Representation falsification: is pathway gamma biology, or just aggregation?

Every representation studied here — biological pathways, size/geometry-matched
random gene sets, and dense random projections — is a single linear map
``W`` of shape ``(k, G)`` applied to gene-level responses as ``R @ W.T``. Putting
them all in one form is the point: it makes the comparison *exactly* matched on
output dimensionality and on the arithmetic applied, so any remaining difference
is attributable to **which genes are grouped together**, not to how many
dimensions survive.

Null constructions, and exactly what each preserves
---------------------------------------------------

``permuted``
    Randomly permutes the *gene axis* of the real membership matrix. Because
    ``(M P)(M P)^T = M M^T`` for a permutation ``P``, this preserves **exactly**:
    the number of sets, every set size, every pairwise set-set overlap, and the
    multiset of per-gene membership counts (how many sets each gene belongs to).
    It destroys only which specific genes are grouped together. This is the
    strictest available null and is the primary one.

``resampled``
    Draws each set independently and uniformly from the shared gene universe
    with the observed set sizes. Preserves the number of sets and the size
    distribution. Does **not** preserve set-set overlap structure or the
    per-gene membership-count distribution — genes become near-uniformly used.

``gaussian``
    A dense i.i.d. normal projection to the same output dimensionality.
    Preserves only the dimensionality. Answers "are k-dimensional vectors simply
    easier?" rather than anything about gene grouping.

None of these constructions may be adjusted after seeing a result.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

EPS = 1e-12


# --------------------------------------------------------------------------
# vectorised correlation (the battery runs hundreds of replicates)
# --------------------------------------------------------------------------


def rowwise_pearson(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Pearson correlation along the last axis, broadcast over leading axes.

    Returns NaN wherever either row is constant, matching the scalar helpers used
    elsewhere rather than silently emitting a divide-by-zero.
    """
    x = np.asarray(a, dtype=np.float64)
    y = np.asarray(b, dtype=np.float64)
    if x.shape != y.shape:
        raise ValueError(f"shape mismatch {x.shape} vs {y.shape}")
    xc = x - x.mean(axis=-1, keepdims=True)
    yc = y - y.mean(axis=-1, keepdims=True)
    num = np.einsum("...g,...g->...", xc, yc)
    den = np.sqrt(np.einsum("...g,...g->...", xc, xc) * np.einsum("...g,...g->...", yc, yc))
    out = np.full(num.shape, np.nan)
    ok = den > EPS
    out[ok] = num[ok] / den[ok]
    return out


# --------------------------------------------------------------------------
# representations as linear maps
# --------------------------------------------------------------------------


def membership_to_weights(membership: np.ndarray) -> np.ndarray:
    """Turn a binary membership matrix into a mean-aggregation map ``W``."""
    M = np.asarray(membership, dtype=np.float64)
    counts = M.sum(axis=1, keepdims=True)
    if np.any(counts <= 0):
        raise ValueError("A gene set has no member genes.")
    return M / counts


def permuted_membership(membership: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Gene-label permutation null. Preserves sizes, overlaps and gene degrees."""
    M = np.asarray(membership, dtype=np.float64)
    return M[:, rng.permutation(M.shape[1])]


def resampled_membership(
    sizes: Sequence[int], n_genes: int, rng: np.random.Generator
) -> np.ndarray:
    """Size-matched independent resampling null. Preserves sizes only."""
    M = np.zeros((len(sizes), n_genes))
    for i, s in enumerate(sizes):
        if s > n_genes:
            raise ValueError(f"set size {s} exceeds the gene universe {n_genes}")
        M[i, rng.choice(n_genes, size=int(s), replace=False)] = 1.0
    return M


def gaussian_projection(k: int, n_genes: int, rng: np.random.Generator) -> np.ndarray:
    """Dense random projection to ``k`` dimensions. Preserves dimensionality only."""
    W = rng.normal(size=(k, n_genes))
    return W / (np.linalg.norm(W, axis=1, keepdims=True) + EPS)


def project(responses: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Apply a linear representation to the last axis of ``responses``."""
    return np.asarray(responses, dtype=np.float64) @ np.asarray(weights, dtype=np.float64).T


# --------------------------------------------------------------------------
# decomposition and recoverability, vectorised
# --------------------------------------------------------------------------


def gamma_of(tensor: np.ndarray) -> np.ndarray:
    """Interaction component of a balanced ``(C, P, K)`` tensor."""
    D = np.asarray(tensor, dtype=np.float64)
    mu = D.mean(axis=(0, 1))
    alpha = D.mean(axis=1) - mu
    beta = D.mean(axis=0) - mu
    return D - mu[None, None, :] - alpha[:, None, :] - beta[None, :, :]


def centred(matrix: np.ndarray) -> np.ndarray:
    """Remove the per-representation mean over perturbations."""
    m = np.asarray(matrix, dtype=np.float64)
    return m - m.mean(axis=0, keepdims=True)


def gamma_recoverability(
    projected: np.ndarray,
    target: int,
    sources: Sequence[int],
    predictor_index: int | None = None,
    weights: np.ndarray | None = None,
) -> np.ndarray:
    """Per-perturbation ``r(gamma_true, gamma_hat)`` in a projected space.

    ``gamma_true`` is the interaction of the **full** four-context projected
    tensor (evaluation only). ``gamma_hat`` is the difference between the centred
    prediction and the centred source mean, both built from source contexts only.
    Supply either ``predictor_index`` (a single source) or ``weights`` (a
    combination over ``sources``).
    """
    Dp = np.asarray(projected, dtype=np.float64)
    gamma_true = gamma_of(Dp)[target]
    A = Dp[list(sources)].mean(axis=0)
    if predictor_index is not None:
        P = Dp[predictor_index]
    elif weights is not None:
        P = np.tensordot(np.asarray(weights, dtype=np.float64), Dp[list(sources)], axes=(0, 0))
    else:
        raise ValueError("Provide predictor_index or weights.")
    return rowwise_pearson(gamma_true, centred(P) - centred(A))


def gamma_reliability(projected_halves: np.ndarray, target: int) -> np.ndarray:
    """Split-half reliability of projected gamma, averaged over repeats.

    ``projected_halves`` has shape ``(repeats, 2, C, P, K)``. Every context is
    split — gamma depends on all of them, so splitting only the target would
    leave an identical source component in both halves and inflate the estimate.
    """
    H = np.asarray(projected_halves, dtype=np.float64)
    per_repeat = [
        rowwise_pearson(gamma_of(H[r, 0])[target], gamma_of(H[r, 1])[target])
        for r in range(H.shape[0])
    ]
    return np.nanmean(np.vstack(per_repeat), axis=0)


def spearman_brown(rho_half: np.ndarray) -> np.ndarray:
    r = np.asarray(rho_half, dtype=np.float64)
    return 2.0 * r / (1.0 + r)


def empirical_percentile(observed: float, null_values: Sequence[float]) -> dict[str, float]:
    """Where the observed value sits in a null distribution.

    Reports the one-sided p-value with the standard ``(count + 1) / (n + 1)``
    correction, which never returns exactly zero and is the honest statement for
    a finite number of replicates.
    """
    nulls = np.asarray([v for v in null_values if np.isfinite(v)], dtype=np.float64)
    if nulls.size == 0 or not np.isfinite(observed):
        return {
            "null_mean": np.nan,
            "null_sd": np.nan,
            "percentile": np.nan,
            "p_value": np.nan,
            "z": np.nan,
            "n_null": 0,
        }
    sd = float(nulls.std(ddof=1)) if nulls.size > 1 else np.nan
    n_ge = int(np.sum(nulls >= observed))
    return {
        "null_mean": float(nulls.mean()),
        "null_sd": sd,
        "percentile": float((nulls < observed).mean() * 100.0),
        "p_value": float((n_ge + 1) / (nulls.size + 1)),
        "z": float((observed - nulls.mean()) / sd) if sd and sd > EPS else np.nan,
        "n_null": int(nulls.size),
    }
