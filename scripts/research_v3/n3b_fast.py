"""Memory-light, algebraically identical implementations of E0s, E2, E4 and metrics for N3-B.

Implementation detail only (n3_protocol.md §3 defines *what* is computed). Each
function reproduces ``virtual_cell.analysis.calibration_budget`` exactly up to
floating-point rounding; ``verify_n3b_fast.py`` checks this on real data before
N3-B is run.

* A subset is held as its consensus mean ``A``, template ``t_A``, frozen scale
  and the perturbation kernel Gram ``ZZᵀ`` (instead of full source arrays).
* E4's kernel ridge uses one eigendecomposition of the anchor Gram for all λ.
* Metrics use cached evaluation dot products: for any prediction ``P``,
  Σ⟨e1−P, e2−P⟩ = Σ⟨e1, e2⟩ − ⟨P, e1+e2⟩ + ‖P‖², and likewise after centring and
  after removing the per-perturbation consensus component.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from virtual_cell.analysis import calibration_budget as cb

EPS = cb.EPS


@dataclass(frozen=True)
class LightView:
    A: np.ndarray
    t_A: np.ndarray
    scale: float
    gram: np.ndarray  # (P, P) = Z Zᵀ


def light_view(sv: cb.SourceView) -> LightView:
    return LightView(A=sv.A, t_A=sv.t_A, scale=sv.scale, gram=sv.Z @ sv.Z.T)


def e0s(lv: LightView, rows: np.ndarray) -> np.ndarray:
    return lv.t_A + lv.scale * (lv.A[rows] - lv.t_A)


def _e2_params(lv: LightView, K: np.ndarray, Y_K: np.ndarray) -> tuple[np.ndarray, float]:
    At = lv.A[K] - lv.t_A
    if len(K) < 2:
        s2 = lv.scale
    else:
        Atc = At - At.mean(axis=0)
        Yc = Y_K - Y_K.mean(axis=0)
        den = float(np.sum(Atc * Atc))
        s2 = float(np.sum(Yc * Atc)) / den if den > EPS else lv.scale
    return Y_K.mean(axis=0) - s2 * At.mean(axis=0), s2


def e2(lv: LightView, rows: np.ndarray, K: np.ndarray, Y_K: np.ndarray) -> np.ndarray:
    c, s2 = _e2_params(lv, K, Y_K)
    return c + s2 * (lv.A[rows] - lv.t_A)


def e4(lv: LightView, rows: np.ndarray, K: np.ndarray, Y_K: np.ndarray) -> tuple[np.ndarray, float]:
    k = len(K)
    R = Y_K - e2(lv, K, K, Y_K)
    Kmat = lv.gram[np.ix_(K, K)]
    Kt = lv.gram[np.ix_(rows, K)]
    lam_vals, Q = np.linalg.eigh(Kmat)
    B = Q.T @ R  # (k, G)
    if k < cb.SMALL_K:
        lam_sel = cb.LAMBDA_GRID[-1]
    else:
        M = B @ B.T
        best = (np.inf, None)
        for lam in cb.LAMBDA_GRID:
            d = 1.0 / (lam_vals + lam)
            W = Q * d  # Q diag(d)
            diag_ginv = np.einsum("ij,ij->i", W, Q)
            alpha_sq = np.einsum("ij,jk,ik->i", W, M, W)
            err = float(np.sum(alpha_sq / diag_ginv**2))
            if err < best[0]:
                best = (err, lam)
        lam_sel = best[1]
    R_hat = (Kt @ (Q * (1.0 / (lam_vals + lam_sel)))) @ B
    return e2(lv, rows, K, Y_K) + R_hat, float(lam_sel)


@dataclass(frozen=True)
class FastEval:
    s_raw: np.ndarray
    c0: float
    sc: np.ndarray
    c1: float
    s_perp: dict
    c3: dict
    u: dict


def fast_eval(e1: np.ndarray, e2_: np.ndarray, consensus: dict[str, np.ndarray]) -> FastEval:
    """``consensus`` maps a name to the test-set source mean whose direction defines γ⊥."""
    e1c, e2c = e1 - e1.mean(axis=0), e2_ - e2_.mean(axis=0)
    sc = e1c + e2c
    c1 = float(np.sum(e1c * e2c))
    s_perp, c3, us = {}, {}, {}
    for name, A_T in consensus.items():
        Ac = A_T - A_T.mean(axis=0)
        u = Ac / np.maximum(np.linalg.norm(Ac, axis=1, keepdims=True), EPS)
        a1, a2 = np.einsum("pg,pg->p", e1c, u), np.einsum("pg,pg->p", e2c, u)
        s_perp[name] = sc - (a1 + a2)[:, None] * u
        c3[name] = c1 - float(np.sum(a1 * a2))
        us[name] = u
    return FastEval(
        s_raw=e1 + e2_, c0=float(np.sum(e1 * e2_)), sc=sc, c1=c1, s_perp=s_perp, c3=c3, u=us
    )


def fast_metrics(fe: FastEval, P: np.ndarray) -> dict[str, float]:
    Pc = P - P.mean(axis=0)
    pp = float(np.sum(Pc * Pc))
    out = {
        "M0": float((-np.sum(P * fe.s_raw) + np.sum(P * P)) / fe.c0) * -1.0,
        "M1": float((np.sum(Pc * fe.sc) - pp) / fe.c1),
    }
    for name, u in fe.u.items():
        pu = np.einsum("pg,pg->p", Pc, u)
        num_minus_c3 = -float(np.sum(Pc * fe.s_perp[name])) + pp - float(np.sum(pu * pu))
        out[f"M3_{name}"] = -num_minus_c3 / fe.c3[name]
    return out
