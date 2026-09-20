"""Reliability-aware targets, selective prediction, and outer-target isolation.

Every reliability claim is validated against synthetic data with a known latent
signal, known noise and a known true D, rather than asserted.
"""

from __future__ import annotations

import numpy as np
import pytest

from virtual_cell.modelling import transferability as tf

N_C, N_P, N_G = 4, 80, 60


def _latent_world(n_p=400, n_g=200, noise_sd=1.0, seed=0, signal_sd=1.0):
    rng = np.random.default_rng(seed)
    latent = rng.normal(scale=signal_sd, size=(n_p, n_g))
    h1 = latent + rng.normal(scale=noise_sd, size=latent.shape)
    h2 = latent + rng.normal(scale=noise_sd, size=latent.shape)
    return latent, h1, h2, rng


# --- unbiasedness ----------------------------------------------------------


def test_signal_energy_is_unbiased_for_the_latent_norm():
    for noise_sd in (0.5, 1.0, 2.0):
        latent, h1, h2, _ = _latent_world(noise_sd=noise_sd, seed=1)
        truth = np.einsum("pg,pg->p", latent, latent)
        assert np.mean(tf.signal_energy(h1, h2)) == pytest.approx(np.mean(truth), rel=0.05)


def test_naive_energy_is_biased_upward_but_cross_half_is_not():
    latent, h1, h2, _ = _latent_world(noise_sd=1.5, seed=2)
    truth = float(np.mean(np.einsum("pg,pg->p", latent, latent)))
    naive = float(np.mean(np.einsum("pg,pg->p", h1, h1)))
    cross = float(np.mean(tf.signal_energy(h1, h2)))
    assert naive > truth * 1.2
    assert cross == pytest.approx(truth, rel=0.05)


def test_residual_energy_is_unbiased_for_the_true_squared_error():
    latent, h1, h2, rng = _latent_world(noise_sd=1.0, seed=3)
    pred = 0.6 * latent + 0.1 * rng.normal(size=latent.shape)
    truth = float(np.mean(np.einsum("pg,pg->p", latent - pred, latent - pred)))
    assert float(np.mean(tf.residual_energy(h1, h2, pred))) == pytest.approx(truth, rel=0.05)


def test_residual_energy_is_near_zero_for_a_perfect_prediction():
    latent, h1, h2, _ = _latent_world(noise_sd=1.0, seed=4)
    got = float(np.mean(tf.residual_energy(h1, h2, latent)))
    scale = float(np.mean(tf.signal_energy(h1, h2)))
    assert abs(got) < 0.05 * scale


def test_energy_helpers_reject_shape_mismatch():
    with pytest.raises(ValueError, match="shape mismatch"):
        tf.signal_energy(np.zeros((3, 4)), np.zeros((3, 5)))
    with pytest.raises(ValueError, match="shape mismatch"):
        tf.residual_energy(np.zeros((3, 4)), np.zeros((3, 4)), np.zeros((3, 5)))


# --- D behaviour -----------------------------------------------------------


@pytest.mark.parametrize("share,expected", [(0.5, 0.25), (0.75, 0.0625), (0.0, 1.0)])
def test_normalised_d_recovers_planted_values(share, expected):
    latent, h1, h2, _ = _latent_world(n_p=800, noise_sd=0.8, seed=5)
    d, stable = tf.normalised_d(h1, h2, share * latent, min_signal=-np.inf)
    assert np.nanmean(d) == pytest.approx(expected, abs=0.05)


def test_d_equals_one_for_a_zero_prediction():
    latent, h1, h2, _ = _latent_world(n_p=600, noise_sd=0.8, seed=6)
    d, _ = tf.normalised_d(h1, h2, np.zeros_like(latent), min_signal=-np.inf)
    assert np.nanmean(d) == pytest.approx(1.0, abs=0.05)


def test_d_exceeds_one_when_worse_than_predicting_zero():
    """The D>1 regime is real and must be representable, not clipped."""
    latent, h1, h2, _ = _latent_world(n_p=600, noise_sd=0.5, seed=7)
    d, _ = tf.normalised_d(h1, h2, -1.5 * latent, min_signal=-np.inf)
    assert np.nanmedian(d) > 1.0


def test_d_is_never_clipped():
    latent, h1, h2, _ = _latent_world(n_p=300, noise_sd=0.5, seed=8)
    d, _ = tf.normalised_d(h1, h2, -3.0 * latent, min_signal=-np.inf)
    assert np.nanmax(d) > 2.0


def test_near_zero_signal_is_masked_and_would_otherwise_be_unstable():
    """A mixed population: strong-signal rows stay, null rows are excluded.

    The null rows have no reproducible response but the source-only predictor
    still predicts *something* for them - which is exactly the real situation.
    Their denominator approaches zero while their numerator does not, so D
    explodes. Masking must remove them and leave the strong rows bounded.
    """
    rng = np.random.default_rng(9)
    n_g = 100
    latent = np.vstack([rng.normal(scale=1.0, size=(250, n_g)), np.zeros((250, n_g))])
    # a prediction that is non-zero for every row, including the null ones
    pred = 0.5 * latent + 0.5 * rng.normal(scale=1.0, size=latent.shape)
    h1 = latent + rng.normal(scale=1.0, size=latent.shape)
    h2 = latent + rng.normal(scale=1.0, size=latent.shape)

    d_unmasked, _ = tf.normalised_d(h1, h2, pred, min_signal=-np.inf)
    d_masked, stable = tf.normalised_d(h1, h2, pred, min_signal=30.0)

    assert stable[:250].mean() > 0.9  # strong rows kept
    assert stable[250:].mean() < 0.1  # null rows dropped
    assert np.isnan(d_masked[~stable]).all()
    # masked values are bounded; unmasked ones are not
    assert np.nanmax(np.abs(d_masked[stable])) < 10.0
    assert np.nanmax(np.abs(d_unmasked)) > 10 * np.nanmax(np.abs(d_masked[stable]))


def test_stability_rule_is_derived_and_monotone():
    rule = tf.derive_stability_rule(n_genes=100, noise_sd=1.0, n_rep=150, seed=3)
    assert np.isfinite(rule["min_signal_energy"])
    curve = rule["curve"]
    first = np.mean([c["within_tolerance"] for c in curve[:5]])
    last = np.mean([c["within_tolerance"] for c in curve[-5:]])
    assert last > first  # stronger signal -> more reliable D


# --- frozen source agreement ----------------------------------------------


@pytest.fixture
def stats():
    rng = np.random.default_rng(11)
    pair = rng.uniform(-1, 1, size=(N_C, N_C, N_P))
    for a in range(N_C):
        for b in range(N_C):
            pair[b, a] = pair[a, b]
    return tf.SourceStats(
        pair_corr=pair,
        reliability=rng.uniform(0, 1, size=(N_C, N_P)),
        cells=rng.integers(30, 500, size=(N_C, N_P)).astype(float),
        magnitude=rng.uniform(1, 8, size=(N_C, N_P)),
        gene_basal=rng.normal(size=(N_C, N_P)),
        basal_sim=rng.uniform(0.85, 0.95, size=(N_C, N_C)),
    )


def test_source_agreement_is_the_mean_of_pairwise_correlations(stats):
    src = [1, 2, 3]
    expected = np.nanmean(
        np.vstack([stats.pair_corr[1, 2], stats.pair_corr[1, 3], stats.pair_corr[2, 3]]), axis=0
    )
    np.testing.assert_allclose(tf.source_agreement(stats, src), expected, atol=1e-12)


def test_source_agreement_ignores_the_target(stats):
    a = tf.source_agreement(stats, [1, 2, 3])
    b = tf.source_agreement(stats, [3, 2, 1])
    np.testing.assert_allclose(a, b, atol=1e-12)


def test_source_agreement_needs_two_sources(stats):
    with pytest.raises(ValueError, match="at least two"):
        tf.source_agreement(stats, [0])


# --- leakage ---------------------------------------------------------------


@pytest.fixture
def world(stats):
    rng = np.random.default_rng(12)
    Y = rng.normal(size=(N_C, N_P, N_G))
    halves = rng.normal(size=(4, 2, N_C, N_P, N_G))
    return Y, halves, stats


def test_features_never_read_the_target_responses(world):
    """LEAKAGE TEST: scramble Y[target] and the target's halves."""
    Y, halves, stats = world
    rng = np.random.default_rng(13)
    for target in range(N_C):
        src = [c for c in range(N_C) if c != target]
        before = tf.build_confidence_features(Y, halves, target, src, stats=stats)
        Yt = Y.copy()
        Yt[target] = rng.normal(size=Y.shape[1:]) * 500
        Ht = halves.copy()
        Ht[:, :, target] = rng.normal(size=halves[:, :, target].shape) * 500
        after = tf.build_confidence_features(Yt, Ht, target, src, stats=stats)
        np.testing.assert_allclose(before, after, atol=0, rtol=0)


def test_features_do_change_when_a_source_changes(world):
    Y, halves, stats = world
    before = tf.build_confidence_features(Y, halves, 0, [1, 2, 3], stats=stats)
    Yt = Y.copy()
    Yt[1] += 5.0
    after = tf.build_confidence_features(Yt, halves, 0, [1, 2, 3], stats=stats)
    assert not np.allclose(before, after)


def test_feature_matrix_shape_and_finiteness(world):
    Y, halves, stats = world
    X = tf.build_confidence_features(Y, halves, 0, [1, 2, 3], stats=stats)
    assert X.shape == (N_P, len(tf.FEATURE_NAMES))
    assert np.all(np.isfinite(X))


def test_every_ablation_uses_only_declared_features():
    for name, feats in tf.ABLATIONS.items():
        assert set(feats) <= set(tf.FEATURE_NAMES), name
        assert "source_agreement" in feats, name


# --- models ----------------------------------------------------------------


def test_standardiser_fits_on_training_rows_only():
    rng = np.random.default_rng(14)
    train = rng.normal(loc=3.0, scale=2.0, size=(400, 5))
    test = rng.normal(loc=40.0, scale=2.0, size=(50, 5))
    sc = tf.Standardiser().fit(train)
    np.testing.assert_allclose(sc.transform(train).mean(axis=0), 0.0, atol=1e-10)
    assert np.all(sc.transform(test).mean(axis=0) > 5.0)


def test_standardiser_refuses_before_fit():
    with pytest.raises(RuntimeError, match="before fit"):
        tf.Standardiser().transform(np.zeros((2, 2)))


def test_ridge_recovers_a_planted_signal():
    rng = np.random.default_rng(15)
    X = rng.normal(size=(1500, 4))
    beta = np.array([1.0, -2.0, 0.5, 0.0])
    y = 3.0 + X @ beta + rng.normal(scale=0.05, size=1500)
    coef = tf.fit_ridge(X, y, alpha=1e-6)
    np.testing.assert_allclose(coef[1:], beta, atol=0.02)
    np.testing.assert_allclose(coef[0], 3.0, atol=0.02)


def test_monotone_calibration_is_monotone():
    rng = np.random.default_rng(16)
    x = rng.uniform(0, 1, size=500)
    y = x**2 + rng.normal(scale=0.05, size=500)
    iso = tf.fit_monotone(x, y)
    grid = np.linspace(0, 1, 50)
    pred = iso.predict(grid)
    assert np.all(np.diff(pred) >= -1e-9)


# --- selective prediction --------------------------------------------------


def test_risk_coverage_decreases_for_a_perfect_confidence_score():
    rng = np.random.default_rng(17)
    risk = rng.exponential(size=600)
    conf = -risk  # perfect ordering
    rc = tf.risk_coverage(conf, risk)
    means = [r["mean_risk"] for r in rc]
    assert all(means[i] >= means[i + 1] for i in range(len(means) - 1))


def test_risk_coverage_is_flat_for_a_useless_confidence_score():
    rng = np.random.default_rng(18)
    risk = rng.exponential(size=2000)
    conf = rng.normal(size=2000)  # independent of risk
    rc = tf.risk_coverage(conf, risk)
    means = [r["mean_risk"] for r in rc]
    assert abs(means[0] - means[-1]) < 0.3 * means[0]


def test_aurc_prefers_a_better_confidence_score():
    rng = np.random.default_rng(19)
    risk = rng.exponential(size=800)
    good = tf.area_under_risk_coverage(-risk, risk)
    bad = tf.area_under_risk_coverage(rng.normal(size=800), risk)
    assert good < bad


def test_aurc_is_nan_for_empty_input():
    assert np.isnan(tf.area_under_risk_coverage(np.array([]), np.array([])))


def test_prioritisation_captures_more_than_random_with_a_good_score():
    rng = np.random.default_rng(20)
    risk = rng.exponential(size=1000)
    out = tf.prioritisation_capture(-risk, risk)
    for row in out:
        assert row["captured_fraction"] > row["random_expectation"]


def test_prioritisation_matches_random_for_a_useless_score():
    rng = np.random.default_rng(21)
    risk = rng.exponential(size=4000)
    out = tf.prioritisation_capture(rng.normal(size=4000), risk)
    for row in out:
        assert abs(row["captured_fraction"] - row["random_expectation"]) < 0.06


def test_prioritisation_clips_negative_risk_only_for_accounting():
    risk = np.array([-1.0, 2.0, 3.0, 5.0])
    out = tf.prioritisation_capture(np.array([4.0, 3.0, 2.0, 1.0]), risk, budgets=(0.5,))
    # worst-confidence half is [5.0, 3.0]; total clipped = 10.0
    assert out[0]["captured_fraction"] == pytest.approx(0.8)


# --- calibration -----------------------------------------------------------


def test_calibration_curve_is_monotone_for_a_faithful_score():
    rng = np.random.default_rng(22)
    conf = rng.uniform(0, 1, size=1200)
    quality = conf + rng.normal(scale=0.05, size=1200)
    curve = tf.calibration_curve(conf, quality, n_bins=5)
    means = [b["mean_quality"] for b in curve]
    assert all(means[i] < means[i + 1] for i in range(len(means) - 1))


def test_calibration_curve_handles_small_input():
    assert tf.calibration_curve(np.array([1.0, 2.0]), np.array([1.0, 2.0]), n_bins=5) == []
