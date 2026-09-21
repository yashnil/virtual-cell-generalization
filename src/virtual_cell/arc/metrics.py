"""A local reimplementation of the six scored ``vcc2026`` metrics.

Definitions follow the published metric reference (``cell-eval2``
``docs/vcc2026_metrics/vcc2026-metrics-brief.md``, 2026/08/19, competition
``rule_version`` 3) and the packaged ``configs/vcc2026.yaml``. Parameter
defaults below are that config's values, written out rather than inherited so
that reading this file tells you what is being computed.

This exists so a candidate predictor can be measured before it is submitted.
It is not the official scorer: the official one is authoritative, and any
disagreement is this module's bug. The invariants the reference states
analytically (a zero-effect prediction scoring exactly 0.5 on ``pds_cosine``, a
zero-fold-change prediction scoring exactly 1 on ``de_wilcoxon_lfc_nmae``) are
asserted in the test-suite, which is the check that this transcription is
faithful.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats

__all__ = [
    "BULK_TARGET_SUM",
    "CELL_TARGET_SUM",
    "MIN_CPM",
    "P_ADJ_THRESHOLD",
    "REACH_PURITY_FLOOR",
    "MIN_GATE_SIZE",
    "LFC_EPSILON",
    "SCORED_METRICS",
    "DETable",
    "bulk_profile",
    "jackknife_dispersion",
    "pds_cosine",
    "expr_mse_unbiased_capped",
    "de_table",
    "direction_fidelity_yield",
    "direction_reach",
    "sig_jaccard",
    "lfc_nmae",
    "scale_score",
]

#: ``bulk_target_sum`` — the per-group normalisation of sections 1 and 2.
BULK_TARGET_SUM = 5.0e4
#: ``target_sum`` — the per-cell normalisation of the differential-expression table.
CELL_TARGET_SUM = 1.0e6
#: ``filter_gene_min_cpm_cell`` — the low-expression gate, read on control cells only.
MIN_CPM = 5.0
#: ``p_adj_threshold`` — Benjamini-Hochberg, applied per perturbation.
P_ADJ_THRESHOLD = 0.05
#: ``REACH_PURITY_FLOOR``
REACH_PURITY_FLOOR = 0.9
#: ``min_gate_size`` for ``de_wilcoxon_lfc_nmae``.
MIN_GATE_SIZE = 10
#: ``epsilon`` in the log2 fold change.
LFC_EPSILON = 1.0e-9

#: The six members of the ``vcc2026`` profile, and whether higher is better.
SCORED_METRICS: dict[str, bool] = {
    "pds_cosine": True,
    "expr_mse_unbiased_capped_norm": False,
    "de_wilcoxon_direction_fidelity_yield_raw": True,
    "de_wilcoxon_direction_reach_raw": True,
    "de_wilcoxon_sig_jaccard": True,
    "de_wilcoxon_lfc_nmae": False,
}


def bulk_profile(counts: np.ndarray, *, target_sum: float = BULK_TARGET_SUM) -> np.ndarray:
    """Pseudobulk a group of cells: sum, normalise to ``target_sum``, ``log1p``."""
    totals = np.asarray(counts, dtype=np.float64).sum(axis=0)
    denom = totals.sum()
    if denom <= 0:
        return np.zeros_like(totals)
    return np.log1p(target_sum * totals / denom)


def jackknife_dispersion(counts: np.ndarray, *, target_sum: float = BULK_TARGET_SUM) -> float:
    """Delete-one jackknife estimate of a group profile's sampling variance.

    ``C_p = (n-1)/n * sum_g sum_i (v_ig - vbar_g)^2`` with ``v_ig`` the profile
    recomputed with cell ``i`` removed. Zero for fewer than two cells.
    """
    y = np.asarray(counts, dtype=np.float64)
    n = y.shape[0]
    if n < 2:
        return 0.0
    totals = y.sum(axis=0)
    lib = y.sum(axis=1)
    s = lib.sum()
    denom = s - lib
    if np.any(denom <= 0):
        return 0.0
    v = np.log1p(target_sum * (totals[None, :] - y) / denom[:, None])
    centred = v - v.mean(axis=0, keepdims=True)
    return float((n - 1) / n * np.sum(centred**2))


def _cosine_distance(pred: np.ndarray, real: np.ndarray) -> np.ndarray:
    """``1 - cos`` between every predicted row and every measured row.

    A row of zero norm sits at distance 1 from everything, per the reference.
    """
    pn = np.linalg.norm(pred, axis=1)
    rn = np.linalg.norm(real, axis=1)
    d = np.ones((pred.shape[0], real.shape[0]), dtype=np.float64)
    ok_p = pn > 0
    ok_r = rn > 0
    if ok_p.any() and ok_r.any():
        sim = (pred[ok_p] / pn[ok_p, None]) @ (real[ok_r] / rn[ok_r, None]).T
        block = 1.0 - sim
        rows = np.flatnonzero(ok_p)[:, None]
        cols = np.flatnonzero(ok_r)[None, :]
        d[rows, cols] = block
    return d


def pds_cosine(
    pred_delta: np.ndarray,
    real_delta: np.ndarray,
    *,
    exclude: np.ndarray | None = None,
) -> np.ndarray:
    """Per-perturbation discrimination score.

    ``pred_delta`` and ``real_delta`` are ``(n_perturbations, n_genes)``
    effects against the *reference* control profile, in matching row order.
    ``exclude`` is a boolean gene mask of the 300 panel target genes, removed
    from every distance so the feature space is identical for every pair.

    Returns ``1 - k_p/(n-1)`` per perturbation, with midrank ties.
    """
    pred = np.asarray(pred_delta, dtype=np.float64)
    real = np.asarray(real_delta, dtype=np.float64)
    if pred.shape != real.shape:
        raise ValueError("pred_delta and real_delta must have the same shape")
    n = pred.shape[0]
    if n < 2:
        raise ValueError("pds_cosine needs at least two perturbations")
    if exclude is not None:
        keep = ~np.asarray(exclude, dtype=bool)
        pred, real = pred[:, keep], real[:, keep]

    d = _cosine_distance(pred, real)
    own = np.diagonal(d)[:, None]
    closer = (d < own).sum(axis=1)
    tied = (d == own).sum(axis=1) - 1
    k = closer + 0.5 * tied
    return 1.0 - k / (n - 1)


@dataclass(frozen=True)
class ExprMse:
    """The parts of ``expr_mse_unbiased_capped``, kept separable for auditing."""

    numerator: np.ndarray
    denominator: np.ndarray
    rho: float

    @property
    def value(self) -> float:
        return float(self.numerator.sum() / self.denominator.sum())


def expr_mse_unbiased_capped(
    pred_profiles: np.ndarray,
    real_profiles: np.ndarray,
    ctrl_profile: np.ndarray,
    *,
    pred_dispersion: np.ndarray,
    real_dispersion: np.ndarray,
    ctrl_dispersion: float,
    target_gene: np.ndarray | None = None,
) -> ExprMse:
    """Sampling-noise-corrected expression error, as a fraction of the measured effect.

    ``target_gene[p]`` is the column index of perturbation ``p``'s own target
    gene, excluded from its two squared distances; ``-1`` excludes nothing.
    """
    pred = np.asarray(pred_profiles, dtype=np.float64)
    real = np.asarray(real_profiles, dtype=np.float64)
    ctrl = np.asarray(ctrl_profile, dtype=np.float64)
    n_pert, n_genes = pred.shape
    c_hat = np.asarray(pred_dispersion, dtype=np.float64)
    c_real = np.asarray(real_dispersion, dtype=np.float64)

    keep = np.ones((n_pert, n_genes), dtype=bool)
    if target_gene is not None:
        tg = np.asarray(target_gene, dtype=np.int64)
        has = tg >= 0
        keep[np.flatnonzero(has), tg[has]] = False
    g_p = keep.sum(axis=1).astype(np.float64)

    err = (((pred - real) ** 2) * keep).sum(axis=1)
    eff = (((real - ctrl[None, :]) ** 2) * keep).sum(axis=1)
    credit = np.minimum(c_hat, c_real)

    # rho bounds the credited correction by the submission's own spread.
    w = 1.0 / g_p
    wsum = (w[:, None] * keep).sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        wmean = np.where(wsum > 0, ((w[:, None] * pred) * keep).sum(axis=0) / wsum, 0.0)
    b = float((w[:, None] * ((pred - wmean[None, :]) ** 2) * keep).sum())
    spread = float((credit / g_p).sum())
    rho = 1.0 if spread <= 0 else min(1.0, b / spread)

    numerator = (err - rho * credit - c_real) / g_p
    denominator = (eff - c_real - ctrl_dispersion) / g_p
    return ExprMse(numerator=numerator, denominator=denominator, rho=rho)


def _benjamini_hochberg(p: np.ndarray) -> np.ndarray:
    """BH-adjusted p-values over a one-dimensional array."""
    p = np.asarray(p, dtype=np.float64)
    m = p.size
    if m == 0:
        return p
    order = np.argsort(p, kind="stable")
    ranked = p[order] * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(adj)
    out[order] = np.clip(adj, 0.0, 1.0)
    return out


@dataclass(frozen=True)
class DETable:
    """A differential-expression table, computed identically for both sides.

    Rows are perturbations, columns the genes surviving the control-side
    low-expression gate. ``tested`` is that gate as a mask on the full axis.
    """

    tested: np.ndarray
    pval: np.ndarray
    p_adj: np.ndarray
    lfc: np.ndarray

    @property
    def significant(self) -> np.ndarray:
        return self.p_adj < P_ADJ_THRESHOLD

    @property
    def adjudicable(self) -> np.ndarray:
        """Genes the *reference* assigns a defined, non-zero fold change."""
        return np.isfinite(self.lfc) & (self.lfc != 0.0)


def _cpm(counts: np.ndarray) -> np.ndarray:
    x = np.asarray(counts, dtype=np.float64)
    lib = x.sum(axis=1, keepdims=True)
    return np.divide(x * CELL_TARGET_SUM, lib, out=np.zeros_like(x), where=lib > 0)


def de_table(
    pert_counts: list[np.ndarray],
    ctrl_counts: np.ndarray,
    *,
    tested: np.ndarray | None = None,
) -> DETable:
    """Wilcoxon rank-sum table for each perturbation against the control group.

    ``tested`` is the control-side low-expression gate. It is an argument
    because the reference's control defines it for *both* sides; pass the
    reference's mask when scoring a prediction.
    """
    ctrl = _cpm(ctrl_counts)
    if tested is None:
        tested = ctrl.mean(axis=0) > MIN_CPM
    tested = np.asarray(tested, dtype=bool)
    ctrl_t = ctrl[:, tested]
    ctrl_mean = ctrl_t.mean(axis=0)

    n_pert = len(pert_counts)
    pval = np.ones((n_pert, int(tested.sum())), dtype=np.float64)
    p_adj = np.ones_like(pval)
    lfc = np.zeros_like(pval)
    for i, block in enumerate(pert_counts):
        x = _cpm(block)[:, tested]
        with np.errstate(invalid="ignore"):
            res = stats.mannwhitneyu(
                x, ctrl_t, axis=0, alternative="two-sided", method="asymptotic"
            )
        pv = np.nan_to_num(np.asarray(res.pvalue, dtype=np.float64), nan=1.0)
        pval[i] = pv
        p_adj[i] = _benjamini_hochberg(pv)
        lfc[i] = np.log2((x.mean(axis=0) + LFC_EPSILON) / (ctrl_mean + LFC_EPSILON))
    return DETable(tested=tested, pval=pval, p_adj=p_adj, lfc=lfc)


def _drop_target(mask: np.ndarray, target: int) -> np.ndarray:
    if target < 0:
        return mask
    out = mask.copy()
    out[target] = False
    return out


def direction_fidelity_yield(
    pred: DETable, real: DETable, *, target_gene: np.ndarray | None = None
) -> np.ndarray:
    """``F_p = k / max(n_pred, n_real)``; NaN where the reference found nothing."""
    n_pert = real.lfc.shape[0]
    tg = np.full(n_pert, -1) if target_gene is None else np.asarray(target_gene)
    out = np.full(n_pert, np.nan)
    adjud = real.adjudicable
    for p in range(n_pert):
        real_sig = _drop_target(real.significant[p], tg[p])
        pred_sig = _drop_target(pred.significant[p] & adjud[p], tg[p])
        n_real, n_pred = int(real_sig.sum()), int(pred_sig.sum())
        if n_real == 0 and n_pred == 0:
            continue
        same = np.sign(pred.lfc[p]) == np.sign(real.lfc[p])
        k = int((pred_sig & same).sum())
        out[p] = k / max(n_pred, n_real)
    return out


def direction_reach(
    pred: DETable, real: DETable, *, target_gene: np.ndarray | None = None
) -> np.ndarray:
    """``R_p = k*/n_real`` — the deepest prefix of the submission's own ranking
    over the reference's significant genes whose directional purity clears 0.9."""
    n_pert = real.lfc.shape[0]
    tg = np.full(n_pert, -1) if target_gene is None else np.asarray(target_gene)
    out = np.full(n_pert, np.nan)
    for p in range(n_pert):
        pool = np.flatnonzero(_drop_target(real.significant[p], tg[p]))
        n_real = pool.size
        if n_real == 0:
            continue
        order = np.lexsort(
            (
                pool,
                -np.abs(pred.lfc[p][pool]),
                pred.pval[p][pool],
                pred.p_adj[p][pool],
                ~pred.significant[p][pool],
            )
        )
        ranked = pool[order]
        keep = real.adjudicable[p][ranked]
        ranked = ranked[keep]
        if ranked.size == 0:
            out[p] = 0.0
            continue
        match = np.sign(pred.lfc[p][ranked]) == np.sign(real.lfc[p][ranked])
        purity = np.cumsum(match) / np.arange(1, ranked.size + 1)
        ok = np.flatnonzero(purity >= REACH_PURITY_FLOOR)
        out[p] = 0.0 if ok.size == 0 else (ok[-1] + 1) / n_real
    return out


def sig_jaccard(
    pred: DETable, real: DETable, *, target_gene: np.ndarray | None = None
) -> np.ndarray:
    """``J_p`` over the significant sets; an empty union is 1."""
    n_pert = real.lfc.shape[0]
    tg = np.full(n_pert, -1) if target_gene is None else np.asarray(target_gene)
    out = np.empty(n_pert)
    for p in range(n_pert):
        a = _drop_target(real.significant[p], tg[p])
        b = _drop_target(pred.significant[p], tg[p])
        union = int((a | b).sum())
        out[p] = 1.0 if union == 0 else int((a & b).sum()) / union
    return out


def lfc_nmae(pred: DETable, real: DETable, *, target_gene: np.ndarray | None = None) -> np.ndarray:
    """Normalised mean absolute fold-change error over the reference's own gate.

    NaN where the gate holds fewer than ``MIN_GATE_SIZE`` genes.
    """
    n_pert = real.lfc.shape[0]
    tg = np.full(n_pert, -1) if target_gene is None else np.asarray(target_gene)
    out = np.full(n_pert, np.nan)
    for p in range(n_pert):
        gate = _drop_target(real.significant[p] & np.isfinite(real.lfc[p]), tg[p])
        if int(gate.sum()) < MIN_GATE_SIZE:
            continue
        truth = real.lfc[p][gate]
        guess = np.nan_to_num(pred.lfc[p][gate], nan=0.0, posinf=0.0, neginf=0.0)
        denom = np.abs(truth).sum()
        if denom <= 0 or not np.isfinite(denom):
            continue
        out[p] = float(np.abs(guess - truth).sum() / denom)
    return out


def scale_score(value: float, baseline: float, replicate: float) -> float:
    """``s = (u - b) / (r - b)``: 0 is the mean-response baseline, 1 a replicate."""
    span = replicate - baseline
    if span == 0:
        raise ValueError("replicate and baseline anchors coincide; score is undefined")
    return float((value - baseline) / span)
