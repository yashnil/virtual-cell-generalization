"""The v2 clean target: beta really is removed, and nothing else changed.

The first four tests are the scientific justification for v2 existing at all —
they prove on planted balanced data that the unshrunk centred residual is a clean
multiple of gamma while the scaled one is not.
"""

from __future__ import annotations

import numpy as np
import pytest

from virtual_cell.modelling import pathway_gamma_v2 as v2
from virtual_cell.modelling import pathway_residual as pr

N_C, N_P, N_K = 4, 50, 11


@pytest.fixture
def planted():
    rng = np.random.default_rng(77)
    mu = rng.normal(size=N_K)
    alpha = rng.normal(size=(N_C, N_K))
    alpha -= alpha.mean(axis=0)
    beta = rng.normal(size=(N_P, N_K)) * 2.0
    beta -= beta.mean(axis=0)
    gamma = rng.normal(size=(N_C, N_P, N_K)) * 0.6
    gamma -= gamma.mean(axis=0, keepdims=True)
    gamma -= gamma.mean(axis=1, keepdims=True)
    Y = mu + alpha[:, None, :] + beta[None, :, :] + gamma
    return Y, mu, alpha, beta, gamma


@pytest.fixture
def world(planted):
    Y = planted[0]
    rng = np.random.default_rng(78)
    stats = pr.GeneLevelStats(
        pair_corr=rng.uniform(-1, 1, size=(N_C, N_C, N_P)),
        reliability=rng.uniform(0, 1, size=(N_C, N_P)),
        cells=rng.integers(30, 400, size=(N_C, N_P)).astype(float),
        magnitude=rng.uniform(1, 8, size=(N_C, N_P)),
        gene_basal=rng.normal(size=(N_C, N_P)),
    )
    extras = {
        "basal_pathway": rng.normal(size=(N_C, N_K)),
        "in_pathway": rng.integers(0, 2, size=(N_P, N_K)).astype(float),
        "pathway_size": rng.integers(10, 150, size=N_K).astype(float),
    }

    def weight_fn(target, sources):
        base = np.arange(1, len(sources) + 1, dtype=float) + target
        return base / base.sum()

    return Y, stats, extras, weight_fn


# --- the algebra that justifies v2 ----------------------------------------


def test_unshrunk_centred_residual_is_exactly_four_thirds_gamma(planted):
    """The whole point of v2: with all three other contexts as sources."""
    Y, _, _, _, gamma = planted
    for target in range(N_C):
        sources = [c for c in range(N_C) if c != target]
        np.testing.assert_allclose(
            v2.clean_gamma_target(Y, target, sources),
            (4.0 / 3.0) * gamma[target],
            atol=1e-10,
        )


def test_clean_target_removes_beta_exactly(planted):
    """Scaling beta up must not change the clean target at all."""
    Y, mu, alpha, beta, gamma = planted
    target, sources = 0, [1, 2, 3]
    before = v2.clean_gamma_target(Y, target, sources)
    Y_big_beta = mu + alpha[:, None, :] + (50.0 * beta)[None, :, :] + gamma
    after = v2.clean_gamma_target(Y_big_beta, target, sources)
    np.testing.assert_allclose(before, after, atol=1e-9)


def test_scaled_residual_still_contains_beta_when_s_is_not_one(planted):
    """The v1 target does move when beta moves — the contamination is real."""
    Y, mu, alpha, beta, gamma = planted
    target, sources = 0, [1, 2, 3]
    s = 0.55
    v1_target = pr.centre_predictions(Y[target] - pr.baseline(Y, sources, s))
    Y_big_beta = mu + alpha[:, None, :] + (3.0 * beta)[None, :, :] + gamma
    v1_after = pr.centre_predictions(Y_big_beta[target] - pr.baseline(Y_big_beta, sources, s))
    assert not np.allclose(v1_target, v1_after, atol=1e-6)
    # and the difference is exactly the extra beta
    np.testing.assert_allclose(v1_after - v1_target, (1 - s) * 2.0 * beta, atol=1e-9)


def test_scaled_residual_equals_clean_target_only_when_s_is_one(planted):
    Y, _, _, _, gamma = planted
    target, sources = 2, [0, 1, 3]
    np.testing.assert_allclose(
        pr.centre_predictions(Y[target] - pr.baseline(Y, sources, 1.0)),
        v2.clean_gamma_target(Y, target, sources),
        atol=1e-10,
    )


def test_two_source_clean_target_is_beta_free_but_not_a_multiple_of_gamma(planted):
    """Inner pseudo-targets: beta still cancels, scale differs. Documented."""
    Y, _, _, beta, gamma = planted
    target, sources = 0, [1, 2]  # third context (3) excluded
    got = v2.clean_gamma_target(Y, target, sources)
    expected = gamma[target] - 0.5 * (gamma[1] + gamma[2])
    np.testing.assert_allclose(got, expected, atol=1e-10)
    assert not np.allclose(got, (4.0 / 3.0) * gamma[target], atol=1e-3)


def test_clean_target_is_centred_by_construction(planted):
    Y = planted[0]
    np.testing.assert_allclose(v2.clean_gamma_target(Y, 1, [0, 2, 3]).mean(axis=0), 0.0, atol=1e-12)


# --- theory coefficient ----------------------------------------------------


def test_lambda_theory_formula():
    assert v2.lambda_theory(0.0) == pytest.approx(0.75)
    assert v2.lambda_theory(1.0) == pytest.approx(1.0)
    assert v2.lambda_theory(0.55) == pytest.approx(0.8875)


def test_lambda_theory_exactly_corrects_the_gamma_deficit(planted):
    """lambda_theory removes the ENTIRE gamma deficit, leaving only alpha and beta.

    ``Y - B = (4/3)alpha + (1-s)beta + (1 + s/3)gamma`` and
    ``lambda_theory * (4/3) = 1 + s/3`` exactly, so the gamma term cancels and
    the remainder is the template plus the deliberately shrunk conserved effect.
    """
    Y, mu, alpha, beta, gamma = planted
    target, sources = 3, [0, 1, 2]
    s = 0.6
    B = pr.baseline(Y, sources, s)
    assert v2.lambda_theory(s) * (4.0 / 3.0) == pytest.approx(1.0 + s / 3.0)
    corrected = B + v2.lambda_theory(s) * (4.0 / 3.0) * gamma[target]
    expected = (4.0 / 3.0) * alpha[target][None, :] + (1 - s) * beta
    np.testing.assert_allclose(Y[target] - corrected, expected, atol=1e-9)
    # and no trace of gamma remains
    residual_centred = pr.centre_predictions(Y[target] - corrected)
    np.testing.assert_allclose(residual_centred, (1 - s) * beta, atol=1e-9)


# --- leakage ---------------------------------------------------------------


def test_outer_prediction_never_reads_the_target_responses(world):
    """LEAKAGE TEST: scramble Y[target]; prediction must be bit-identical."""
    Y, stats, extras, weight_fn = world
    target, sources = 0, [1, 2, 3]
    sel = v2.Selection("M2", 100.0, 0.75)
    kw = dict(stats=stats, weight_fn=weight_fn, **extras)
    out = v2.fit_and_predict_outer(Y, target, sources, sel, **kw)
    tampered = Y.copy()
    tampered[target] = np.random.default_rng(1).normal(size=Y.shape[1:]) * 900
    out2 = v2.fit_and_predict_outer(tampered, target, sources, sel, **kw)
    np.testing.assert_allclose(out["prediction"], out2["prediction"], atol=0, rtol=0)
    np.testing.assert_allclose(out["gamma_hat"], out2["gamma_hat"], atol=0, rtol=0)
    assert out["scale"] == out2["scale"]
    assert out["lambda_theory"] == out2["lambda_theory"]


def test_inner_selection_never_reads_the_outer_target(world):
    """LEAKAGE TEST: lambda and family selection must not move."""
    Y, stats, extras, weight_fn = world
    target, sources = 3, [0, 1, 2]
    kw = dict(stats=stats, weight_fn=weight_fn, **extras)
    sel = v2.select_by_inner_loco(Y, sources, lambdas=(0.0, 0.5, 1.0), alphas=(10.0, 1000.0), **kw)
    tampered = Y.copy()
    tampered[target] = np.random.default_rng(2).normal(size=Y.shape[1:]) * 400
    sel2 = v2.select_by_inner_loco(
        tampered, sources, lambdas=(0.0, 0.5, 1.0), alphas=(10.0, 1000.0), **kw
    )
    assert (sel.family, sel.alpha, sel.lam) == (sel2.family, sel2.alpha, sel2.lam)
    assert sel.inner_scores == sel2.inner_scores


def test_training_targets_use_only_source_contexts(world):
    """Rows for pseudo-target c must ignore the outer target entirely."""
    Y, stats, extras, weight_fn = world
    spec = v2.FoldSpec(1, (2, 3))  # outer target 0 excluded from sources
    kw = dict(stats=stats, weight_fn=weight_fn, **extras)
    X, R, B, s = v2.assemble_rows(Y, spec, **kw)
    tampered = Y.copy()
    tampered[0] = np.random.default_rng(3).normal(size=Y.shape[1:]) * 500
    X2, R2, B2, s2 = v2.assemble_rows(tampered, spec, **kw)
    np.testing.assert_allclose(X, X2, atol=0, rtol=0)
    np.testing.assert_allclose(R, R2, atol=0, rtol=0)
    np.testing.assert_allclose(B, B2, atol=0, rtol=0)
    assert s == s2


def test_lambda_zero_returns_the_baseline_exactly(world):
    Y, stats, extras, weight_fn = world
    out = v2.fit_and_predict_outer(
        Y, 0, [1, 2, 3], v2.Selection("M2", 100.0, 0.0), stats=stats, weight_fn=weight_fn, **extras
    )
    np.testing.assert_allclose(out["prediction"], out["baseline"], atol=0, rtol=0)


def test_features_are_identical_to_v1(world):
    """v2 must change ONLY the target — the design matrix must be byte-identical."""
    Y, stats, extras, weight_fn = world
    spec = v2.FoldSpec(1, (0, 2, 3))
    kw = dict(stats=stats, weight_fn=weight_fn, **extras)
    X_v2, _, _, _ = v2.assemble_rows(Y, spec, **kw)
    X_v1, _, _, _ = pr.assemble_rows(Y, spec, **kw)
    np.testing.assert_allclose(X_v2, X_v1, atol=0, rtol=0)


def test_v2_target_differs_from_v1_target(world):
    """...and the target must actually differ, or v2 is pointless."""
    Y, stats, extras, weight_fn = world
    spec = v2.FoldSpec(1, (0, 2, 3))
    kw = dict(stats=stats, weight_fn=weight_fn, **extras)
    _, R_v2, _, _ = v2.assemble_rows(Y, spec, **kw)
    _, R_v1, _, _ = pr.assemble_rows(Y, spec, **kw)
    assert not np.allclose(R_v2, R_v1)


# --- bootstrap -------------------------------------------------------------


def test_bootstrap_delta_is_zero_for_identical_predictions():
    rng = np.random.default_rng(4)
    Y = rng.normal(size=(200, 12))
    P = 0.5 * Y + rng.normal(scale=0.3, size=Y.shape)
    out = v2.bootstrap_delta(Y, P, P, n_boot=200)
    assert out["delta_median_pearson"] == pytest.approx(0.0, abs=1e-12)
    assert out["lo"] <= 0.0 <= out["hi"]
    assert out["frac_improved"] == 0.0


def test_bootstrap_delta_detects_a_real_improvement():
    rng = np.random.default_rng(5)
    Y = rng.normal(size=(400, 15))
    bad = rng.normal(size=Y.shape)
    good = 0.8 * Y + 0.2 * rng.normal(size=Y.shape)
    out = v2.bootstrap_delta(Y, bad, good, n_boot=400)
    assert out["delta_median_pearson"] > 0.5
    assert out["lo"] > 0
    assert out["frac_improved"] > 0.95


def test_bootstrap_delta_reports_sample_size():
    rng = np.random.default_rng(6)
    Y = rng.normal(size=(50, 9))
    out = v2.bootstrap_delta(Y, Y, Y, n_boot=50)
    assert out["n"] == 50
