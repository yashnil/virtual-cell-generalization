"""C2: low-complexity cell realisation of a fixed mean prediction (ours, competition track).

Every generator here takes the same expected composition the C1 atlas predicts and only
changes how integer cells are drawn from it. Predeclared in
``reports/competition_v2/c2_predeclaration.md``.

* ``g1_counts``: real control donors, binomial thinning for down-regulated genes and
  depth-proportional Poisson additions for up-regulated ones. The identity at ``r = 1``.
* ``g2_counts``: multinomial at the target composition with the donors' depths.
* ``g3_counts``: gamma-Poisson (negative binomial) at the target composition with per-gene
  dispersion fitted on control cells.
* ``ideal_de_table``: no cells; an analytic Welch test at the predicted means under the
  reference's own per-gene noise (the mean-response ceiling of section D).

The frozen G0 (C1 dual-moment emitter) stays in :mod:`generator`. The arc-track
generators in :mod:`virtual_cell.arc.generate` are unrelated and unchanged.
"""

from __future__ import annotations

import numpy as np
from scipy import stats

from virtual_cell.arc import metrics as M

#: Mask below which a composition entry is treated as absent when forming ratios.
_TINY = 1e-300


def effect_ratio(target: np.ndarray, control: np.ndarray) -> np.ndarray:
    """Per-gene multiplier ``r = target / control``; 1 where the control is zero."""
    target = np.asarray(target, dtype=np.float64)
    control = np.asarray(control, dtype=np.float64)
    if target.shape != control.shape:
        raise ValueError("target and control compositions differ in shape")
    return np.divide(target, control, out=np.ones_like(target), where=control > _TINY)


def g1_counts(
    donors: np.ndarray, ratio: np.ndarray, control: np.ndarray, *, rng: np.random.Generator
) -> np.ndarray:
    """Transport real donor cells by a per-gene multiplier without redrawing them.

    ``donors`` are raw integer control cells. Genes with ``ratio < 1`` are binomially
    thinned (``Bin(x, r)``); genes with ``ratio > 1`` receive
    ``Poisson((r - 1) * depth_i * control_g)`` extra counts, so an up-regulated gene can
    appear in a cell that had none. Both steps are unbiased for ``r * E[x]`` when
    ``control`` is the donors' mean composition. At ``ratio == 1`` the cell is unchanged.
    """
    x = np.asarray(donors)
    r = np.asarray(ratio, dtype=np.float64)
    c = np.asarray(control, dtype=np.float64)
    if x.ndim != 2 or r.shape != (x.shape[1],) or c.shape != r.shape:
        raise ValueError("donor, ratio and control axes differ")
    if (r < 0).any() or not np.isfinite(r).all() or (c < 0).any():
        raise ValueError("ratio and control must be finite and nonnegative")
    if (x < 0).any():
        raise ValueError("donors must be nonnegative counts")
    x = x.astype(np.int64)
    out = x.copy()
    down = np.flatnonzero(r < 1)
    if down.size:
        out[:, down] = rng.binomial(x[:, down], r[down][None, :])
    up = np.flatnonzero(r > 1)
    if up.size:
        depth = x.sum(axis=1).astype(np.float64)
        lam = depth[:, None] * ((r[up] - 1) * c[up])[None, :]
        out[:, up] += rng.poisson(lam)
    return out


def g2_counts(depths: np.ndarray, composition: np.ndarray, *, rng: np.random.Generator):
    """``Multinomial(d_i, p)`` per cell."""
    p = np.asarray(composition, dtype=np.float64)
    if (p < 0).any() or p.sum() <= 0:
        raise ValueError("composition must be nonnegative with positive mass")
    p = p / p.sum()
    d = np.asarray(depths, dtype=np.int64)
    return rng.multinomial(d, p).astype(np.int64)


def fit_dispersion(controls: np.ndarray, composition: np.ndarray) -> np.ndarray:
    """Per-gene NB dispersion ``phi`` by moments around ``mu_ig = d_i * c_g``; floored at 0.

    ``phi_g = max(0, (sum (x - mu)^2 - sum mu) / sum mu^2)``, the pooled
    method-of-moments estimate for ``var = mu + phi * mu^2`` with cell-specific means.
    """
    x = np.asarray(controls, dtype=np.float64)
    c = np.asarray(composition, dtype=np.float64)
    depth = x.sum(axis=1)
    num = np.zeros(x.shape[1])
    den = np.zeros(x.shape[1])
    lin = np.zeros(x.shape[1])
    for left in range(0, len(x), 512):
        blk = x[left : left + 512]
        mu = depth[left : left + 512, None] * c[None, :]
        num += ((blk - mu) ** 2).sum(axis=0)
        lin += mu.sum(axis=0)
        den += (mu**2).sum(axis=0)
    phi = np.divide(num - lin, den, out=np.zeros_like(num), where=den > 0)
    return np.clip(phi, 0.0, None)


def g3_counts(
    depths: np.ndarray, composition: np.ndarray, phi: np.ndarray, *, rng: np.random.Generator
):
    """Gamma-Poisson draws: ``Poisson(d_i p_g G_ig)``, ``G ~ Gamma(1/phi_g, phi_g)``.

    ``phi_g == 0`` degenerates to Poisson. Depth is preserved only in expectation.
    """
    p = np.asarray(composition, dtype=np.float64)
    p = p / p.sum()
    phi = np.asarray(phi, dtype=np.float64)
    mu = np.asarray(depths, dtype=np.float64)[:, None] * p[None, :]
    over = phi > 0
    rate = mu.copy()
    if over.any():
        shape = 1.0 / phi[over]
        rate[:, over] = mu[:, over] * rng.gamma(
            shape[None, :], 1.0 / shape[None, :], size=(len(mu), int(over.sum()))
        )
    return rng.poisson(rate).astype(np.int64)


def ideal_de_table(
    p_cpm: np.ndarray, reference: np.ndarray, tested: np.ndarray, *, n_cells: int = 400
) -> M.DETable:
    """Analytic DE for idealised cells with exact predicted means and real-control noise.

    ``p_cpm`` is ``(n_targets, n_genes)`` mean-CPM compositions on the full fold axis.
    For each tested gene, the predicted per-cell CPM mean is ``1e6 * p``. Its variance is
    the reference per-cell CPM variance scaled by ``(m_pred / m_ctrl)^2`` (a constant
    coefficient of variation). ``z`` is Welch's statistic at ``n_cells`` against the
    reference, the two-sided normal p-value is taken, and BH is applied per perturbation.
    The log fold change follows :func:`virtual_cell.arc.metrics.de_table`.
    """
    tested = np.asarray(tested, dtype=bool)
    ref = M._cpm(reference)[:, tested]
    m_c = ref.mean(axis=0)
    v_c = ref.var(axis=0, ddof=1)
    n_ref = ref.shape[0]
    m_p = 1e6 * np.asarray(p_cpm, dtype=np.float64)[:, tested]
    scale = np.divide(m_p, m_c[None, :], out=np.ones_like(m_p), where=m_c[None, :] > 0)
    v_p = v_c[None, :] * scale**2
    se = np.sqrt(v_p / n_cells + v_c[None, :] / n_ref)
    with np.errstate(invalid="ignore", divide="ignore"):
        z = np.where(se > 0, (m_p - m_c[None, :]) / se, 0.0)
    pval = np.clip(2 * stats.norm.sf(np.abs(z)), 0.0, 1.0)
    p_adj = np.vstack([M._benjamini_hochberg(row) for row in pval])
    lfc = np.log2((m_p + M.LFC_EPSILON) / (m_c[None, :] + M.LFC_EPSILON))
    return M.DETable(tested=tested, pval=pval, p_adj=p_adj, lfc=lfc)


def cell_diagnostics(counts: np.ndarray, *, rng: np.random.Generator, n_pairs_cells: int = 100):
    """Per-group structure: gene-wise mean / variance of raw counts and of log-CPM, and the
    median pairwise cosine distance between cells (log1p CP10K, a fixed-size subsample)."""
    x = np.asarray(counts, dtype=np.float64)
    lib = x.sum(axis=1)
    logcpm = np.log1p(1e4 * x / np.where(lib > 0, lib, 1)[:, None])
    take = rng.choice(len(x), size=min(n_pairs_cells, len(x)), replace=False)
    sub = logcpm[take]
    norms = np.linalg.norm(sub, axis=1)
    sim = (sub / np.maximum(norms, 1e-30)[:, None]) @ (sub / np.maximum(norms, 1e-30)[:, None]).T
    iu = np.triu_indices(len(sub), k=1)
    return {
        "gene_mean": x.mean(axis=0),
        "gene_var": x.var(axis=0, ddof=1),
        "logcpm_var": logcpm.var(axis=0, ddof=1),
        "zero_fraction": float((x == 0).mean()),
        "cosine_distance_median": float(np.median(1 - sim[iu])),
    }
