"""N1 calibration-budget checks F1, F2, F5, F6 (reports/n1_n4_protocol.md §2.9)."""

from __future__ import annotations

import numpy as np
import pytest
from scipy import sparse

from virtual_cell.analysis import calibration_budget as cb

P, G, F = 300, 120, 5
SOURCES = [1, 2, 3]
TARGET = 0


def _world(rng, gamma_scale=1.0, noise=0.05):
    f = rng.normal(size=(P, F))
    W = rng.normal(size=(F, G)) / np.sqrt(F)
    V = rng.normal(size=(F, G)) / np.sqrt(F)
    beta = f @ W
    gamma = gamma_scale * np.tanh(f) @ V  # perturbation-structured, nonlinear in f
    D = np.stack([beta + gamma] + [beta + noise * rng.normal(size=(P, G)) for _ in SOURCES])
    target_latent = beta + gamma
    fit = target_latent + noise * rng.normal(size=(P, G))
    e1 = target_latent + 2 * noise * rng.normal(size=(P, G))
    e2 = target_latent + 2 * noise * rng.normal(size=(P, G))
    return D, fit, e1, e2


def _predict(D, fit, K, T):
    sv = cb.make_source_view(D, SOURCES, scale=0.9)
    return cb.predict_all(sv, T, K, fit[K] if K is not None else None)[0]


def test_f1_non_anchor_and_target_rows_never_read():
    rng = np.random.default_rng(0)
    D, fit, _, _ = _world(rng)
    T, K = np.arange(200, 300), np.arange(0, 40)
    ref = _predict(D, fit, K, T)
    D2, fit2 = D.copy(), fit.copy()
    D2[TARGET] = rng.normal(size=(P, G))  # canonical target row
    mask = np.ones(P, bool)
    mask[K] = False
    fit2[mask] = rng.normal(size=(mask.sum(), G))  # every non-anchor target row
    new = _predict(D2, fit2, K, T)
    for name in ref:
        np.testing.assert_array_equal(ref[name], new[name])


def test_f2_zero_shot_is_source_only():
    rng = np.random.default_rng(1)
    D, fit, _, _ = _world(rng)
    T = np.arange(200, 300)
    ref = _predict(D, fit, None, T)
    D2 = D.copy()
    D2[TARGET] = rng.normal(size=(P, G))
    new = _predict(D2, rng.normal(size=fit.shape), None, T)
    assert set(ref) == {"E0", "E0s"}
    for name in ref:
        np.testing.assert_array_equal(ref[name], new[name])


def _m3(D, fit, e1, e2, K, T, name):
    sv = cb.make_source_view(D, SOURCES, scale=0.9)
    ev = cb.make_eval_context(e1[T], e2[T], sv.A[T])
    P_ = cb.predict_all(sv, T, K, fit[K])[0][name]
    return cb.ratio(*cb.metric_terms(ev, P_)["M3"])


def test_f5_planted_gamma_is_recovered():
    rng = np.random.default_rng(2)
    D, fit, e1, e2 = _world(rng, gamma_scale=1.0)
    T, K = np.arange(200, 300), np.arange(0, 100)
    assert _m3(D, fit, e1, e2, K, T, "E4") > 0.2


@pytest.mark.parametrize("name", ["E3", "E4"])
def test_f5_no_gamma_no_gain(name):
    vals = []
    for seed in range(5):
        rng = np.random.default_rng(10 + seed)
        D, fit, e1, e2 = _world(rng, gamma_scale=0.0)
        T, K = np.arange(200, 300), np.arange(0, 100)
        vals.append(_m3(D, fit, e1, e2, K, T, name))
    assert np.mean(vals) <= 0.02


def test_zero_shot_and_scale_predictors_score_exactly_zero_on_gamma_perp():
    rng = np.random.default_rng(3)
    D, fit, e1, e2 = _world(rng)
    T, K = np.arange(200, 300), np.arange(0, 30)
    sv = cb.make_source_view(D, SOURCES, scale=0.9)
    ev = cb.make_eval_context(e1[T], e2[T], sv.A[T])
    preds = cb.predict_all(sv, T, K, fit[K])[0]
    for name in ("E0", "E0s", "E1", "E2"):
        assert abs(cb.ratio(*cb.metric_terms(ev, preds[name])["M3"])) < 1e-9


def test_template_cannot_move_template_removed_metrics():
    rng = np.random.default_rng(4)
    D, fit, e1, e2 = _world(rng)
    T = np.arange(200, 300)
    sv = cb.make_source_view(D, SOURCES, scale=0.9)
    ev = cb.make_eval_context(e1[T], e2[T], sv.A[T])
    P_ = cb.e0s(sv, T)
    shifted = P_ + rng.normal(size=G)
    a, b = cb.metric_terms(ev, P_), cb.metric_terms(ev, shifted)
    for m in ("M1", "M2", "M3", "M4"):
        assert cb.ratio(*a[m]) == pytest.approx(cb.ratio(*b[m]), abs=1e-10)


def test_split_parts_are_disjoint_and_sized():
    rng = np.random.default_rng(5)
    X = sparse.csr_matrix(rng.random((330, 6)))
    group = np.repeat(np.arange(3), 110)
    C = sparse.csr_matrix(rng.random((101, 6)))
    sm = cb.disjoint_split_means(X, group, C, n_perturbations=3, rng=rng)
    assert (sm.n_pert[0] == 55).all() and (sm.n_pert[1] == 27).all()
    assert sm.n_ctrl.tolist() == [50, 25, 25]
