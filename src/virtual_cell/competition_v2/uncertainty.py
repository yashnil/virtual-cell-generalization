"""C6: gene-level measurement uncertainty and measurement-aware fusion (ours).

Predeclared in ``reports/competition_v2/c6_predeclaration.md``. C1's effects are unchanged;
this module adds:

* per-group cell moments (mean, variance, n) of per-cell CPM, with seeded split halves,
  for in-memory or streamed cells;
* the delta-method SE of C1's own log2fc effect (``delta_sigma2``), and the
  CD4 summary-statistic SE (``cd4_sigma2``);
* panel mapping identical to :func:`atlas.source_effect`'s lookup;
* weighted per-cell fusion (U1 / U2), the pooled-per-gene DerSimonian–Laird ``tau2``, and
  empirical-Bayes multipliers (U3).

With unit weights and no multipliers, :func:`fuse` equals C1 equal fusion up to
floating-point order (tested against :func:`fusion.fused_effects`).
"""

from __future__ import annotations

import numpy as np
from scipy import sparse

from virtual_cell.competition_v2 import fusion
from virtual_cell.competition_v2.sources import MINIMUM_CELLS, PRIOR_COUNTS

LN2 = np.log(2.0)


# ----------------------------------------------------------------------------- moments


def half_assignment(labels: np.ndarray, groups, *, seed: int) -> np.ndarray:
    """Per-cell half (0/1) for every group, in ``groups`` order; -1 for other cells."""
    rng = np.random.default_rng(seed)
    half = np.full(len(labels), -1, dtype=np.int8)
    for g in groups:
        rows = rng.permutation(np.flatnonzero(labels == g))
        h = len(rows) // 2
        half[rows[:h]], half[rows[h:]] = 0, 1
    return half


class Moments:
    """Streaming per-group sums of per-cell CPM and CPM², overall and per half."""

    def __init__(self, groups, n_genes: int):
        self.groups = list(groups)
        G = len(self.groups)
        self.n = np.zeros(G)
        self.s1 = np.zeros((G, n_genes))
        self.s2 = np.zeros((G, n_genes))
        self.hn = np.zeros((G, 2))
        self.h1 = np.zeros((G, 2, n_genes))
        self.h2 = np.zeros((G, 2, n_genes))

    def update(self, x: sparse.spmatrix, ids: np.ndarray, half: np.ndarray) -> None:
        x = sparse.csr_matrix(x, dtype=np.float64)
        lib = np.asarray(x.sum(axis=1)).ravel()
        if (lib <= 0).any():
            raise ValueError("zero-depth cell")
        cpm = sparse.diags(1e6 / lib) @ x
        sq = cpm.multiply(cpm)
        ok = ids >= 0
        G = len(self.groups)

        def member(rows_mask, codes):
            r = np.flatnonzero(rows_mask)
            return sparse.csr_matrix((np.ones(len(r)), (codes[r], r)), shape=(G, len(ids)))

        m = member(ok, ids)
        self.n += np.asarray(m.sum(axis=1)).ravel()
        self.s1 += (m @ cpm).toarray()
        self.s2 += (m @ sq).toarray()
        for h in (0, 1):
            mh = member(ok & (half == h), ids)
            self.hn[:, h] += np.asarray(mh.sum(axis=1)).ravel()
            self.h1[:, h] += (mh @ cpm).toarray()
            self.h2[:, h] += (mh @ sq).toarray()

    @staticmethod
    def _mv(n, s1, s2):
        nn = np.maximum(n, 1)[..., None]
        mean = s1 / nn
        var = np.maximum(s2 - nn * mean**2, 0) / np.maximum(nn - 1, 1)
        return mean, var

    def result(self) -> dict:
        mean, var = self._mv(self.n, self.s1, self.s2)
        hm, hv = self._mv(self.hn, self.h1, self.h2)
        return {
            "groups": np.asarray(self.groups),
            "n": self.n,
            "mean": mean,
            "var": var,
            "half_n": self.hn,
            "half_mean": hm,
            "half_var": hv,
        }


# ----------------------------------------------------------------------------- sigma2


def delta_sigma2(m_shrunk, c, f, v_t, n_t, v_c, n_c) -> np.ndarray:
    """Delta-method variance of ``log2((m'+1)/(c+1))`` with ``m' = f m_t + (1-f) c``."""
    m_shrunk = np.asarray(m_shrunk, dtype=np.float64)
    c = np.asarray(c, dtype=np.float64)
    f = np.asarray(f, dtype=np.float64)
    a = f / ((m_shrunk + 1) * LN2)
    b = ((1 - f) / (m_shrunk + 1) - 1 / (c + 1)) / LN2
    with np.errstate(invalid="ignore", divide="ignore"):
        return a**2 * v_t / np.maximum(n_t, 1) + b**2 * v_c / np.maximum(n_c, 1)


def shrink_fraction(target_count_sums: np.ndarray, prior: float = PRIOR_COUNTS) -> np.ndarray:
    total = np.asarray(target_count_sums, dtype=np.float64).sum(axis=1)
    return total / (total + prior)


def cd4_sigma2(stats: dict, minimum_cells: int = MINIMUM_CELLS) -> np.ndarray:
    """``Σ_c lfcSE²_c / k²`` over the usable conditions, per (source target, source gene)."""
    usable = stats["available"] & stats["quality_pass"] & (stats["n_cells"] >= minimum_cells)
    se = np.asarray(stats["lfcSE"], dtype=np.float64)
    lfc = np.asarray(stats["log2fc"], dtype=np.float64)
    ok = usable[:, :, None] & np.isfinite(se) & np.isfinite(lfc)
    k = ok.sum(axis=0)
    num = np.where(ok, se**2, 0).sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(k > 0, num / np.maximum(k, 1) ** 2, np.nan)


def to_panel(values, src_targets, src_genes, targets, genes, *, row_ok=None, col_ok=None):
    """Map a (source target x source gene) array onto the panel, NaN where absent."""
    ti = {str(t): i for i, t in enumerate(src_targets)}
    gi = {str(g): i for i, g in enumerate(src_genes)}
    rows = np.array([ti.get(str(t), -1) for t in targets])
    cols = np.array([gi.get(str(g), -1) for g in genes])
    if row_ok is not None:
        rows = np.where((rows >= 0) & row_ok[np.maximum(rows, 0)], rows, -1)
    if col_ok is not None:
        cols = np.where((cols >= 0) & col_ok[np.maximum(cols, 0)], cols, -1)
    out = np.full((len(targets), len(genes)), np.nan)
    r, c = np.flatnonzero(rows >= 0), np.flatnonzero(cols >= 0)
    if len(r) and len(c):
        out[np.ix_(r, c)] = np.asarray(values, dtype=np.float64)[np.ix_(rows[r], cols[c])]
    return out


# ----------------------------------------------------------------------------- fusion


def fuse(comp, space: str, control_probability, weights: dict, mult: dict | None = None):
    """Per-cell weighted mean over measuring sources; ``weights``/``mult`` are (T x G)."""
    mult = mult or {}
    parts = comp[space]
    shape = next(iter(parts.values()))[0].shape if parts else comp["CD4"][0].shape
    num = np.zeros(shape)
    den = np.zeros(shape)
    for name, (eff, msk) in parts.items():
        e = eff.astype(np.float64) * mult.get(name, 1.0)
        w = np.where(msk, weights[name], 0.0)
        num += w * e
        den += w
    if "CD4" in comp:
        eff, avail = comp["CD4"]
        e = eff.astype(np.float64) * mult.get("CD4", 1.0)
        moved = e if space == "log2fc" else fusion.cd4_to_bulk_delta(e, control_probability)
        w = np.where(avail, weights["CD4"], 0.0)
        num += w * np.where(avail, moved, 0.0)
        den += w
    return np.divide(num, den, out=np.zeros(shape), where=den > 0)


def pooled_tau2(effects: dict, sigma2: dict, masks: dict) -> np.ndarray:
    """DerSimonian–Laird tau² per gene, pooled over all (p) with >= 2 measuring sources."""
    names = list(effects)
    shape = next(iter(effects.values())).shape
    W = np.zeros((len(names), *shape))
    E = np.zeros_like(W)
    for i, n in enumerate(names):
        ok = masks[n] & np.isfinite(sigma2[n]) & (sigma2[n] > 0)
        W[i] = np.where(ok, 1 / np.where(ok, sigma2[n], 1), 0)
        E[i] = np.where(ok, effects[n], 0)
    k = (W > 0).sum(axis=0)
    S1 = W.sum(axis=0)
    S2 = (W**2).sum(axis=0)
    mean = np.divide((W * E).sum(axis=0), S1, out=np.zeros(shape), where=S1 > 0)
    Q = (W * (E - mean) ** 2).sum(axis=0)
    multi = k >= 2
    num = np.where(multi, Q - (k - 1), 0).sum(axis=0)
    den = np.where(multi, S1 - np.divide(S2, S1, out=np.zeros(shape), where=S1 > 0), 0).sum(axis=0)
    return np.maximum(np.divide(num, den, out=np.zeros(shape[1]), where=den > 0), 0)


def eb_multiplier(effect, sigma2, usable_rows) -> tuple[np.ndarray, np.ndarray]:
    """Per-source EB shrinkage ``τ²/(τ²+σ²)`` with ``τ²_g = max(0, mean e² − mean σ²)``.

    The moments run over ``usable_rows`` (the source's usable retained targets). Cells
    with non-finite σ² get a multiplier of 1 (C1 treatment). Returns
    ``(multiplier, tau2_g)``.
    """
    e = np.asarray(effect, dtype=np.float64)[usable_rows]
    s = np.asarray(sigma2, dtype=np.float64)[usable_rows]
    ok = np.isfinite(s)
    cnt = ok.sum(axis=0)
    me2 = np.divide(np.where(ok, e**2, 0).sum(axis=0), cnt, out=np.zeros(e.shape[1]), where=cnt > 0)
    ms2 = np.divide(np.where(ok, s, 0).sum(axis=0), cnt, out=np.zeros(e.shape[1]), where=cnt > 0)
    tau2 = np.maximum(me2 - ms2, 0)
    full = np.asarray(sigma2, dtype=np.float64)
    with np.errstate(invalid="ignore", divide="ignore"):
        m = tau2[None, :] / (tau2[None, :] + full)
    m = np.where(np.isfinite(full), np.nan_to_num(m, nan=0.0), 1.0)
    return m, tau2
