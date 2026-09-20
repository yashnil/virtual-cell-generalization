"""Nested LOCO, feature construction, scaling, residual algebra and leakage.

The leakage tests are the load-bearing ones: each overwrites the outer target's
responses with noise and requires the pipeline to be bit-identical.
"""

from __future__ import annotations

import numpy as np
import pytest

from virtual_cell.modelling import pathway_residual as pr

N_C, N_P, N_K = 4, 60, 12


@pytest.fixture
def world():
    """A planted four-context pathway world plus all inference-available inputs."""
    rng = np.random.default_rng(31)
    mu = rng.normal(size=N_K)
    alpha = rng.normal(size=(N_C, N_K))
    alpha -= alpha.mean(axis=0)
    beta = rng.normal(size=(N_P, N_K))
    beta -= beta.mean(axis=0)
    gamma = rng.normal(size=(N_C, N_P, N_K)) * 0.5
    gamma -= gamma.mean(axis=0, keepdims=True)
    gamma -= gamma.mean(axis=1, keepdims=True)
    Y = mu + alpha[:, None, :] + beta[None, :, :] + gamma
    stats = pr.GeneLevelStats(
        pair_corr=rng.uniform(-1, 1, size=(N_C, N_C, N_P)),
        reliability=rng.uniform(0, 1, size=(N_C, N_P)),
        cells=rng.integers(30, 500, size=(N_C, N_P)).astype(float),
        magnitude=rng.uniform(1, 8, size=(N_C, N_P)),
        gene_basal=rng.normal(size=(N_C, N_P)),
    )
    extras = {
        "basal_pathway": rng.normal(size=(N_C, N_K)),
        "in_pathway": rng.integers(0, 2, size=(N_P, N_K)).astype(float),
        "pathway_size": rng.integers(10, 200, size=N_K).astype(float),
    }

    def weight_fn(target, sources):
        w = rng_stable(target, sources)
        return w

    return Y, mu, alpha, beta, gamma, stats, extras, weight_fn


def rng_stable(target, sources):
    """Deterministic basal-like weights that never touch responses."""
    base = np.arange(1, len(sources) + 1, dtype=float) + target
    return base / base.sum()


# --- residual algebra ------------------------------------------------------


def test_residual_matches_the_documented_algebra(world):
    """R = (4/3)alpha + (1-s)beta + (1+s/3)gamma, exactly."""
    Y, mu, alpha, beta, gamma, *_ = world
    for target in range(N_C):
        sources = [c for c in range(N_C) if c != target]
        s = pr.fit_scale(Y, sources)
        R = Y[target] - pr.baseline(Y, sources, s)
        expected = (
            (4.0 / 3.0) * alpha[target][None, :]
            + (1.0 - s) * beta
            + (1.0 + s / 3.0) * gamma[target]
        )
        np.testing.assert_allclose(R, expected, atol=1e-9)


def test_centred_residual_removes_alpha_but_keeps_beta(world):
    """Centring kills the template exactly; the beta contamination survives."""
    Y, mu, alpha, beta, gamma, *_ = world
    target, sources = 0, [1, 2, 3]
    s = pr.fit_scale(Y, sources)
    R = Y[target] - pr.baseline(Y, sources, s)
    Rc = pr.centre_predictions(R)
    expected = (1.0 - s) * beta + (1.0 + s / 3.0) * gamma[target]
    np.testing.assert_allclose(Rc, expected, atol=1e-9)
    # and the surviving beta part is not negligible
    assert np.linalg.norm((1.0 - s) * beta) > 0.1 * np.linalg.norm(Rc)


def test_residual_is_not_gamma(world):
    """Guard against ever calling R the ANOVA gamma."""
    Y, _, _, _, gamma, *_ = world
    target, sources = 2, [0, 1, 3]
    R = Y[target] - pr.baseline(Y, sources, pr.fit_scale(Y, sources))
    assert not np.allclose(pr.centre_predictions(R), gamma[target], atol=1e-3)


def test_scale_of_one_reduces_residual_to_four_thirds_alpha_plus_gamma(world):
    Y, _, alpha, _, gamma, *_ = world
    target, sources = 1, [0, 2, 3]
    R = Y[target] - pr.baseline(Y, sources, 1.0)
    expected = (4.0 / 3.0) * (alpha[target][None, :] + gamma[target])
    np.testing.assert_allclose(R, expected, atol=1e-9)


# --- scale calibration -----------------------------------------------------


def test_fit_scale_uses_only_the_given_contexts(world):
    Y, *_ = world
    sources = [1, 2, 3]
    before = pr.fit_scale(Y, sources)
    tampered = Y.copy()
    tampered[0] = np.random.default_rng(0).normal(size=Y.shape[1:]) * 50
    assert pr.fit_scale(tampered, sources) == before


def test_fit_scale_is_one_for_a_single_context(world):
    Y, *_ = world
    assert pr.fit_scale(Y, [0]) == 1.0


def test_fit_scale_recovers_a_planted_scalar():
    rng = np.random.default_rng(5)
    base = rng.normal(size=(30, 8))
    base -= base.mean(axis=0)
    k = 0.6
    Y = np.stack([k * base + rng.normal(scale=1e-9, size=base.shape) for _ in range(3)])
    Y = np.stack([base, k * base, k * base])
    assert 0.3 < pr.fit_scale(Y, [0, 1, 2]) < 1.3


# --- baseline and centring -------------------------------------------------


def test_baseline_never_reads_the_target(world):
    Y, *_ = world
    sources = [0, 1, 2]
    B = pr.baseline(Y, sources, 0.5)
    tampered = Y.copy()
    tampered[3] = np.random.default_rng(1).normal(size=Y.shape[1:]) * 99
    np.testing.assert_allclose(pr.baseline(tampered, sources, 0.5), B, atol=0, rtol=0)


def test_centring_is_prediction_only_and_idempotent():
    rng = np.random.default_rng(2)
    R = rng.normal(size=(40, 9))
    Rc = pr.centre_predictions(R)
    np.testing.assert_allclose(Rc.mean(axis=0), 0.0, atol=1e-12)
    np.testing.assert_allclose(pr.centre_predictions(Rc), Rc, atol=1e-12)


# --- features --------------------------------------------------------------


def test_feature_matrix_shape_and_names(world):
    Y, _, _, _, _, stats, extras, weight_fn = world
    X = pr.build_features(Y, 0, [1, 2, 3], stats=stats, weights=weight_fn(0, [1, 2, 3]), **extras)
    assert X.shape == (N_P * N_K, len(pr.FEATURE_NAMES))
    assert np.all(np.isfinite(X))


def test_features_never_read_the_target_responses(world):
    """LEAKAGE TEST: scramble Y[target]; every feature must be bit-identical."""
    Y, _, _, _, _, stats, extras, weight_fn = world
    for target in range(N_C):
        sources = [c for c in range(N_C) if c != target]
        w = weight_fn(target, sources)
        before = pr.build_features(Y, target, sources, stats=stats, weights=w, **extras)
        tampered = Y.copy()
        tampered[target] = np.random.default_rng(target).normal(size=Y.shape[1:]) * 1000
        after = pr.build_features(tampered, target, sources, stats=stats, weights=w, **extras)
        np.testing.assert_allclose(before, after, atol=0, rtol=0)


def test_features_do_change_when_a_source_changes(world):
    """Sanity: the leakage test above would pass trivially if features were constant."""
    Y, _, _, _, _, stats, extras, weight_fn = world
    sources = [1, 2, 3]
    w = weight_fn(0, sources)
    before = pr.build_features(Y, 0, sources, stats=stats, weights=w, **extras)
    tampered = Y.copy()
    tampered[1] += 5.0
    after = pr.build_features(tampered, 0, sources, stats=stats, weights=w, **extras)
    assert not np.allclose(before, after)


def test_target_basal_features_are_actually_used(world):
    """Basal target information is permitted and must reach the design matrix."""
    Y, _, _, _, _, stats, extras, weight_fn = world
    sources = [1, 2, 3]
    w = weight_fn(0, sources)
    before = pr.build_features(Y, 0, sources, stats=stats, weights=w, **extras)
    bumped = dict(extras)
    bp = extras["basal_pathway"].copy()
    bp[0] += 3.0
    bumped["basal_pathway"] = bp
    after = pr.build_features(Y, 0, sources, stats=stats, weights=w, **bumped)
    assert not np.allclose(before, after)


# --- standardisation -------------------------------------------------------


def test_standardiser_fits_on_training_rows_only():
    rng = np.random.default_rng(3)
    train = rng.normal(loc=5.0, scale=2.0, size=(500, 4))
    test = rng.normal(loc=50.0, scale=20.0, size=(100, 4))
    sc = pr.Standardiser().fit(train)
    np.testing.assert_allclose(sc.transform(train).mean(axis=0), 0.0, atol=1e-10)
    np.testing.assert_allclose(sc.transform(train).std(axis=0), 1.0, atol=1e-10)
    # test rows are NOT re-centred; their shift must survive
    assert np.all(sc.transform(test).mean(axis=0) > 5.0)


def test_standardiser_survives_constant_columns():
    X = np.column_stack([np.ones(50), np.arange(50.0)])
    out = pr.Standardiser().fit(X).transform(X)
    assert np.all(np.isfinite(out))


def test_standardiser_refuses_transform_before_fit():
    with pytest.raises(RuntimeError, match="before fit"):
        pr.Standardiser().transform(np.zeros((2, 2)))


# --- models ----------------------------------------------------------------


def test_ridge_recovers_a_planted_linear_signal_at_low_penalty():
    rng = np.random.default_rng(4)
    X = rng.normal(size=(2000, 5))
    beta = np.array([1.5, -2.0, 0.0, 0.5, 3.0])
    y = 7.0 + X @ beta + rng.normal(scale=0.05, size=2000)
    coef = pr.fit_ridge(X, y, alpha=1e-6)
    np.testing.assert_allclose(coef[0], 7.0, atol=0.02)
    np.testing.assert_allclose(coef[1:], beta, atol=0.02)


def test_ridge_shrinks_towards_zero_as_alpha_grows():
    rng = np.random.default_rng(6)
    X = rng.normal(size=(500, 4))
    y = X @ np.array([2.0, 2.0, 2.0, 2.0]) + rng.normal(size=500)
    small = np.linalg.norm(pr.fit_ridge(X, y, 1.0)[1:])
    large = np.linalg.norm(pr.fit_ridge(X, y, 1e6)[1:])
    assert large < small / 10


def test_ridge_intercept_is_not_penalised():
    rng = np.random.default_rng(7)
    X = rng.normal(size=(400, 3))
    y = 25.0 + rng.normal(scale=0.01, size=400)
    coef = pr.fit_ridge(X, y, alpha=1e8)
    assert coef[0] == pytest.approx(25.0, abs=0.05)


def test_deterministic_fit_recovers_a_planted_scalar():
    rng = np.random.default_rng(8)
    x = rng.normal(size=800)
    assert pr.fit_deterministic(x, -1.75 * x) == pytest.approx(-1.75, abs=1e-9)


def test_m0_predicts_exactly_zero(world):
    Y, _, _, _, _, stats, extras, weight_fn = world
    X = pr.build_features(Y, 0, [1, 2, 3], stats=stats, weights=weight_fn(0, [1, 2, 3]), **extras)
    out = pr._predict_residual("M0", X, None, (N_P, N_K))
    assert np.all(out == 0.0)


# --- nested LOCO -----------------------------------------------------------


def test_inner_selection_never_touches_the_outer_target(world):
    """LEAKAGE TEST: the whole selection must be bit-identical under a scrambled target."""
    Y, _, _, _, _, stats, extras, weight_fn = world
    target, sources = 0, [1, 2, 3]
    kw = dict(stats=stats, weight_fn=weight_fn, **extras)
    sel = pr.select_by_inner_loco(Y, sources, lambdas=(0.0, 0.5, 1.0), alphas=(10.0, 1000.0), **kw)
    tampered = Y.copy()
    tampered[target] = np.random.default_rng(9).normal(size=Y.shape[1:]) * 500
    sel2 = pr.select_by_inner_loco(
        tampered, sources, lambdas=(0.0, 0.5, 1.0), alphas=(10.0, 1000.0), **kw
    )
    assert (sel.family, sel.alpha, sel.lam) == (sel2.family, sel2.alpha, sel2.lam)
    assert sel.inner_scores == sel2.inner_scores


def test_outer_prediction_uses_the_target_only_through_basal(world):
    """LEAKAGE TEST: scrambling Y[target] must not change the prediction at all."""
    Y, _, _, _, _, stats, extras, weight_fn = world
    target, sources = 3, [0, 1, 2]
    sel = pr.Selection(family="M2", alpha=100.0, lam=0.5)
    kw = dict(stats=stats, weight_fn=weight_fn, **extras)
    out = pr.fit_and_predict_outer(Y, target, sources, sel, **kw)
    tampered = Y.copy()
    tampered[target] = np.random.default_rng(10).normal(size=Y.shape[1:]) * 777
    out2 = pr.fit_and_predict_outer(tampered, target, sources, sel, **kw)
    np.testing.assert_allclose(out["prediction"], out2["prediction"], atol=0, rtol=0)
    np.testing.assert_allclose(out["baseline"], out2["baseline"], atol=0, rtol=0)


def test_lambda_zero_returns_exactly_the_baseline(world):
    Y, _, _, _, _, stats, extras, weight_fn = world
    sel = pr.Selection(family="M2", alpha=100.0, lam=0.0)
    out = pr.fit_and_predict_outer(Y, 0, [1, 2, 3], sel, stats=stats, weight_fn=weight_fn, **extras)
    np.testing.assert_allclose(out["prediction"], out["baseline"], atol=0, rtol=0)


def test_selection_can_choose_lambda_zero_when_correction_is_useless(world):
    """The model must have an honest way to decline to correct."""
    Y, _, _, _, _, stats, extras, weight_fn = world
    rng = np.random.default_rng(11)
    noise_stats = pr.GeneLevelStats(
        pair_corr=rng.normal(size=stats.pair_corr.shape),
        reliability=rng.normal(size=stats.reliability.shape),
        cells=rng.normal(size=stats.cells.shape),
        magnitude=rng.normal(size=stats.magnitude.shape),
        gene_basal=rng.normal(size=stats.gene_basal.shape),
    )
    sel = pr.select_by_inner_loco(
        Y,
        [1, 2, 3],
        stats=noise_stats,
        weight_fn=weight_fn,
        families=("M0", "M2"),
        alphas=(1e8,),
        lambdas=(0.0, 1.0),
        **extras,
    )
    assert sel.lam in (0.0, 1.0)
    assert sel.family in ("M0", "M2")


def test_train_rows_come_from_every_source_pseudo_target(world):
    Y, _, _, _, _, stats, extras, weight_fn = world
    sel = pr.Selection(family="M2", alpha=100.0, lam=1.0)
    out = pr.fit_and_predict_outer(Y, 0, [1, 2, 3], sel, stats=stats, weight_fn=weight_fn, **extras)
    assert out["n_train_rows"] == 3 * N_P * N_K


def test_assemble_rows_returns_consistent_pieces(world):
    Y, *_rest, stats, extras, weight_fn = world
    spec = pr.FoldSpec(1, (0, 2, 3))
    X, R, B, scale = pr.assemble_rows(Y, spec, stats=stats, weight_fn=weight_fn, **extras)
    assert X.shape[0] == R.size == N_P * N_K
    np.testing.assert_allclose(R.reshape(N_P, N_K), Y[1] - B, atol=1e-12)
    assert 0.0 < scale < 2.0


# --- metrics ---------------------------------------------------------------


def test_metrics_are_perfect_for_a_perfect_prediction():
    rng = np.random.default_rng(12)
    Y = rng.normal(size=(30, 10))
    m = pr.prediction_metrics(Y, Y)
    assert m["pearson"] == pytest.approx(1.0)
    assert m["cosine"] == pytest.approx(1.0)
    assert m["energy_explained"] == pytest.approx(1.0)
    assert m["mse"] == pytest.approx(0.0, abs=1e-18)


def test_zero_prediction_explains_no_energy():
    rng = np.random.default_rng(13)
    Y = rng.normal(size=(20, 8))
    assert pr.prediction_metrics(Y, np.zeros_like(Y))["energy_explained"] == pytest.approx(0.0)


def test_metrics_reject_shape_mismatch():
    with pytest.raises(ValueError, match="shape mismatch"):
        pr.prediction_metrics(np.zeros((3, 4)), np.zeros((3, 5)))


def test_rowwise_pearson_matches_corrcoef():
    rng = np.random.default_rng(14)
    a, b = rng.normal(size=(11, 20)), rng.normal(size=(11, 20))
    np.testing.assert_allclose(
        pr.rowwise_pearson(a, b),
        [np.corrcoef(a[i], b[i])[0, 1] for i in range(11)],
        atol=1e-12,
    )
