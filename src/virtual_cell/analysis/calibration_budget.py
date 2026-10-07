"""N1 calibration budget: disjoint cell splits, few-shot estimators, γ⊥ metrics.

Preregistered in ``reports/n1_n4_protocol.md``. Every definition here implements
that document. The four-context β/γ are never read; the target context is touched
only through the anchors' fit-part deltas (prediction) and the disjoint evaluation
parts (metrics).
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from scipy import sparse

from virtual_cell.analysis.robustness import aggregate, group_aggregate
from virtual_cell.data.io import DataIntegrityError

LAMBDA_GRID = (0.01, 0.1, 1.0, 10.0, 100.0)
SMALL_K = 5  # below this, the largest λ is used without selection
EPS = 1e-12


# --------------------------------------------------------------------------
# 1. disjoint cell splits (protocol §1)
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class SplitMeans:
    """Part means for one context and one repeat.

    ``pert`` is ``(3, P, G)`` for parts F, E1, E2; ``ctrl`` is ``(3, G)`` for
    CF, CE1, CE2; ``n_pert`` is ``(3, P)`` cell counts; ``n_ctrl`` is ``(3,)``.
    """

    pert: np.ndarray
    ctrl: np.ndarray
    n_pert: np.ndarray
    n_ctrl: np.ndarray


def _three_way(idx: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = len(idx)
    h, q = n // 2, n // 4
    return idx[:h], idx[h : h + q], idx[h + q : h + 2 * q]


def disjoint_split_means(
    X_pert: sparse.csr_matrix,
    group: np.ndarray,
    X_ctrl: sparse.csr_matrix,
    *,
    n_perturbations: int,
    rng: np.random.Generator,
) -> SplitMeans:
    """Split perturbed cells (per perturbation) and control cells into F/E1/E2.

    Controls are permuted and split **before** any delta exists, so no target
    quantity can share control noise across the fit and evaluation sides.
    Disjointness is asserted, not assumed.
    """
    order = np.argsort(group, kind="stable")
    gs = group[order]
    bounds = np.searchsorted(gs, np.arange(n_perturbations + 1))

    parts: list[list[np.ndarray]] = [[], [], []]
    groups: list[list[np.ndarray]] = [[], [], []]
    for p in range(n_perturbations):
        idx = order[bounds[p] : bounds[p + 1]].copy()
        if len(idx) < 4:
            raise DataIntegrityError(f"Perturbation row {p} has {len(idx)} cells; need >= 4.")
        rng.shuffle(idx)
        for j, part in enumerate(_three_way(idx)):
            parts[j].append(part)
            groups[j].append(np.full(len(part), p))

    perm = rng.permutation(X_ctrl.shape[0])
    ctrl_parts = _three_way(perm)
    if min(len(c) for c in ctrl_parts) == 0:
        raise DataIntegrityError("Too few control cells for a three-way split.")

    # disjointness (protocol F4)
    all_pert = [np.concatenate(p) for p in parts]
    for a in range(3):
        for b in range(a + 1, 3):
            if np.intersect1d(all_pert[a], all_pert[b]).size:
                raise AssertionError("perturbed-cell parts overlap")
            if np.intersect1d(ctrl_parts[a], ctrl_parts[b]).size:
                raise AssertionError("control-cell parts overlap")

    pert = np.empty((3, n_perturbations, X_pert.shape[1]))
    n_pert = np.empty((3, n_perturbations), dtype=np.int64)
    for j in range(3):
        sel = all_pert[j]
        grp = np.concatenate(groups[j])
        pert[j] = group_aggregate(X_pert[sel], grp, n_perturbations, "mean_log")
        n_pert[j] = np.bincount(grp, minlength=n_perturbations)
    ctrl = np.vstack([aggregate(X_ctrl[np.sort(c)], "mean_log") for c in ctrl_parts])
    return SplitMeans(
        pert=pert,
        ctrl=ctrl,
        n_pert=n_pert,
        n_ctrl=np.array([len(c) for c in ctrl_parts]),
    )


# --------------------------------------------------------------------------
# 2. estimators (protocol §2.3)
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class SourceView:
    """Everything source-side an estimator may read for one outer fold."""

    S: np.ndarray  # (3, P, G) raw source responses
    S_c: np.ndarray  # (3, P, G) each source centred over the full panel
    A: np.ndarray  # (P, G) source mean
    t_A: np.ndarray  # (G,) source template
    scale: float  # frozen fit_scale
    Z: np.ndarray  # (P, 3G) normalised kernel features


def make_source_view(D: np.ndarray, sources: Sequence[int], scale: float) -> SourceView:
    S = np.asarray(D, dtype=np.float64)[list(sources)]
    S_c = S - S.mean(axis=1, keepdims=True)
    A = S.mean(axis=0)
    norms = np.linalg.norm(S_c, axis=2, keepdims=True)
    Z = (S_c / np.maximum(norms, EPS)).transpose(1, 0, 2).reshape(S.shape[1], -1) / np.sqrt(
        S.shape[0]
    )
    return SourceView(S=S, S_c=S_c, A=A, t_A=A.mean(axis=0), scale=float(scale), Z=Z)


def e0(sv: SourceView, rows: np.ndarray) -> np.ndarray:
    return sv.A[rows]


def e0s(sv: SourceView, rows: np.ndarray) -> np.ndarray:
    return sv.t_A + sv.scale * (sv.A[rows] - sv.t_A)


def e1(sv: SourceView, rows: np.ndarray, K: np.ndarray, Y_K: np.ndarray) -> np.ndarray:
    alpha = (Y_K - e0s(sv, K)).mean(axis=0)
    return e0s(sv, rows) + alpha


def _e2_params(sv: SourceView, K: np.ndarray, Y_K: np.ndarray) -> tuple[np.ndarray, float]:
    At = sv.A[K] - sv.t_A
    if len(K) < 2:
        s2 = sv.scale
    else:
        Atc = At - At.mean(axis=0)
        Yc = Y_K - Y_K.mean(axis=0)
        den = float(np.sum(Atc * Atc))
        s2 = float(np.sum(Yc * Atc)) / den if den > EPS else sv.scale
    c = Y_K.mean(axis=0) - s2 * At.mean(axis=0)
    return c, s2


def e2(sv: SourceView, rows: np.ndarray, K: np.ndarray, Y_K: np.ndarray) -> np.ndarray:
    c, s2 = _e2_params(sv, K, Y_K)
    return c + s2 * (sv.A[rows] - sv.t_A)


def e3(
    sv: SourceView, rows: np.ndarray, K: np.ndarray, Y_K: np.ndarray
) -> tuple[np.ndarray, float]:
    """Gene-wise ridge of target on the 3 centred sources, shrunk toward s/3 each."""
    k = len(K)
    n_s = sv.S_c.shape[0]
    X = sv.S_c[:, K, :].transpose(2, 1, 0)  # (G, k, S)
    w_prior = np.full(n_s, sv.scale / n_s)
    z = Y_K.T - X @ w_prior  # (G, k) residual against the prior
    Xt = np.concatenate([np.ones((X.shape[0], k, 1)), X], axis=2)  # (G, k, S+1)
    Xc = X - X.mean(axis=1, keepdims=True)
    gram_diag = np.einsum("gks,gks->g", Xc, Xc) / n_s  # mean diagonal of centred Gram
    XtX = np.einsum("gki,gkj->gij", Xt, Xt)
    Xtz = np.einsum("gki,gk->gi", Xt, z)
    pen = np.zeros(n_s + 1)
    pen[1:] = 1.0

    def solve(lam_rel: float):
        lam = lam_rel * np.maximum(gram_diag, EPS)
        M = XtX + lam[:, None, None] * np.diag(pen)[None]
        Minv = np.linalg.inv(M)
        coef = np.einsum("gij,gj->gi", Minv, Xtz)
        return coef, Minv

    if k < SMALL_K:
        lam_sel = LAMBDA_GRID[-1]
        coef, _ = solve(lam_sel)
    else:
        best = (np.inf, None, None)
        for lam_rel in LAMBDA_GRID:
            coef, Minv = solve(lam_rel)
            fitted = np.einsum("gki,gi->gk", Xt, coef)
            h = np.einsum("gki,gij,gkj->gk", Xt, Minv, Xt)
            loo = (z - fitted) / np.maximum(1.0 - h, 1e-6)
            err = float(np.sum(loo * loo))
            if err < best[0]:
                best = (err, lam_rel, coef)
        _, lam_sel, coef = best
    w0 = coef[:, 0]
    w = coef[:, 1:] + w_prior  # (G, S)
    pred = w0[None, :] + np.einsum("spg,gs->pg", sv.S_c[:, rows, :], w)
    return pred, float(lam_sel)


def e4(
    sv: SourceView, rows: np.ndarray, K: np.ndarray, Y_K: np.ndarray
) -> tuple[np.ndarray, float]:
    """E2 plus kernel ridge of the anchors' E2 residuals over source similarity."""
    k = len(K)
    R = Y_K - e2(sv, K, K, Y_K)
    Kmat = sv.Z[K] @ sv.Z[K].T
    Kt = sv.Z[rows] @ sv.Z[K].T
    if k < SMALL_K:
        lam_sel = LAMBDA_GRID[-1]
    else:
        best = (np.inf, None)
        for lam in LAMBDA_GRID:
            Ginv = np.linalg.inv(Kmat + lam * np.eye(k))
            alpha = Ginv @ R
            loo = alpha / np.diag(Ginv)[:, None]
            err = float(np.sum(loo * loo))
            if err < best[0]:
                best = (err, lam)
        lam_sel = best[1]
    R_hat = Kt @ np.linalg.solve(Kmat + lam_sel * np.eye(k), R)
    return e2(sv, rows, K, Y_K) + R_hat, float(lam_sel)


ESTIMATORS = ("E0", "E0s", "E1", "E2", "E3", "E4")


def predict_all(
    sv: SourceView, rows: np.ndarray, K: np.ndarray | None, Y_K: np.ndarray | None
) -> tuple[dict[str, np.ndarray], dict[str, float]]:
    """All estimators for one draw. ``K is None`` means k = 0."""
    preds = {"E0": e0(sv, rows), "E0s": e0s(sv, rows)}
    lams: dict[str, float] = {}
    if K is None or len(K) == 0:
        return preds, lams
    preds["E1"] = e1(sv, rows, K, Y_K)
    preds["E2"] = e2(sv, rows, K, Y_K)
    preds["E3"], lams["E3"] = e3(sv, rows, K, Y_K)
    preds["E4"], lams["E4"] = e4(sv, rows, K, Y_K)
    return preds, lams


# --------------------------------------------------------------------------
# 3. metrics (protocol §2.4–2.5)
# --------------------------------------------------------------------------

METRIC_TERMS = ("M0", "M1", "M2", "M3", "M4")


@dataclass(frozen=True)
class EvalContext:
    """Evaluation-only quantities for one (target, repeat) on the test set."""

    e1: np.ndarray  # raw (T, G)
    e2: np.ndarray
    e1c: np.ndarray  # centred over T
    e2c: np.ndarray
    yc: np.ndarray
    Ac: np.ndarray  # source mean centred over T
    u: np.ndarray  # unit conserved directions (T, G)


def _centre(x: np.ndarray) -> np.ndarray:
    return x - x.mean(axis=0, keepdims=True)


def make_eval_context(e1: np.ndarray, e2: np.ndarray, A_T: np.ndarray) -> EvalContext:
    e1c, e2c = _centre(e1), _centre(e2)
    Ac = _centre(A_T)
    u = Ac / np.maximum(np.linalg.norm(Ac, axis=1, keepdims=True), EPS)
    return EvalContext(e1=e1, e2=e2, e1c=e1c, e2c=e2c, yc=(e1c + e2c) / 2, Ac=Ac, u=u)


def _par(x: np.ndarray, u: np.ndarray) -> np.ndarray:
    return np.einsum("pg,pg->p", x, u)[:, None] * u


def _rowwise_r(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    ac, bc = _centre(a.T).T, _centre(b.T).T
    num = np.einsum("pg,pg->p", ac, bc)
    den = np.linalg.norm(ac, axis=1) * np.linalg.norm(bc, axis=1)
    out = np.where(den > EPS, num / np.maximum(den, EPS), 0.0)
    return out


def metric_terms(ev: EvalContext, P: np.ndarray) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Per-perturbation (numerator, denominator) for M0–M4; ratio = 1 − Σnum/Σden."""
    Pc = _centre(P)
    out: dict[str, tuple[np.ndarray, np.ndarray]] = {}

    def nd(a1, a2, p):
        return np.einsum("pg,pg->p", a1 - p, a2 - p), np.einsum("pg,pg->p", a1, a2)

    out["M0"] = nd(ev.e1, ev.e2, P)
    num1, den1 = nd(ev.e1c, ev.e2c, Pc)
    out["M1"] = (num1, den1)
    e1p, e2p, Pp = _par(ev.e1c, ev.u), _par(ev.e2c, ev.u), _par(Pc, ev.u)
    out["M2"] = nd(e1p, e2p, Pp)
    out["M3"] = nd(ev.e1c - e1p, ev.e2c - e2p, Pc - Pp)
    # F6: decomposition identity
    for i in (0, 1):
        tot = out["M2"][i].sum() + out["M3"][i].sum()
        ref = out["M1"][i].sum()
        if abs(tot - ref) > 1e-8 * max(1.0, abs(ref)):
            raise AssertionError("M2 + M3 does not reproduce M1 (F6)")
    # M4: residual of plain source-mean transfer
    num4 = num1
    den4 = np.einsum("pg,pg->p", ev.e1c - ev.Ac, ev.e2c - ev.Ac)
    out["M4"] = (num4, den4)
    return out


def ratio(num: np.ndarray, den: np.ndarray) -> float:
    d = float(np.sum(den))
    return 1.0 - float(np.sum(num)) / d if abs(d) > EPS else float("nan")


def correlation_terms(ev: EvalContext, P: np.ndarray) -> dict[str, float]:
    Pc = _centre(P)
    Pp = _par(Pc, ev.u)
    yp = _par(ev.yc, ev.u)
    r_full = _rowwise_r(Pc, ev.yc)
    P_orth = Pc - Pp
    has_orth = np.linalg.norm(P_orth, axis=1) > 1e-10 * np.maximum(np.linalg.norm(Pc, axis=1), EPS)
    r_orth = np.where(has_orth, _rowwise_r(P_orth, ev.yc - yp), 0.0)
    return {"M1r": float(np.median(r_full)), "M3r": float(np.median(r_orth))}


def bootstrap_ratio(
    num: np.ndarray, den: np.ndarray, *, n_boot: int, rng: np.random.Generator
) -> tuple[float, float]:
    """95 % interval of ``1 − Σnum/Σden`` over resampled perturbations."""
    n = len(num)
    idx = rng.integers(0, n, size=(n_boot, n))
    vals = 1.0 - num[idx].sum(axis=1) / den[idx].sum(axis=1)
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))
