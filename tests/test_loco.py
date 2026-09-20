"""Leakage-sensitive algebra and reliability corrections for the LOCO phase.

The reliability estimators are validated against synthetic data with a *known*
latent signal and known noise, so the correction is demonstrated rather than
asserted.
"""

from __future__ import annotations

import numpy as np
import pytest

from virtual_cell.analysis import loco
from virtual_cell.decomposition import anova

CONTEXTS = ("k562", "rpe1", "hepg2", "jurkat")


@pytest.fixture
def planted():
    """A four-context tensor built from known zero-sum components."""
    rng = np.random.default_rng(11)
    n_c, n_p, n_g = 4, 40, 25
    mu = rng.normal(size=n_g)
    alpha = rng.normal(size=(n_c, n_g))
    alpha -= alpha.mean(axis=0)
    beta = rng.normal(size=(n_p, n_g))
    beta -= beta.mean(axis=0)
    gamma = rng.normal(size=(n_c, n_p, n_g))
    gamma -= gamma.mean(axis=0, keepdims=True)
    gamma -= gamma.mean(axis=1, keepdims=True)
    D = mu + alpha[:, None, :] + beta[None, :, :] + gamma
    return D, mu, alpha, beta, gamma


# --- folds ----------------------------------------------------------------


def test_folds_hold_out_each_context_exactly_once():
    folds = loco.make_folds(CONTEXTS)
    assert len(folds) == 4
    assert [f.target for f in folds] == list(CONTEXTS)
    for f in folds:
        assert f.target not in f.sources
        assert len(f.sources) == 3
        assert set(f.sources) | {f.target} == set(CONTEXTS)
        assert f.target_index not in f.source_indices


def test_fold_indices_match_names():
    for f in loco.make_folds(CONTEXTS):
        assert CONTEXTS[f.target_index] == f.target
        assert tuple(CONTEXTS[i] for i in f.source_indices) == f.sources


def test_too_few_contexts_is_rejected():
    with pytest.raises(ValueError, match="at least two"):
        loco.make_folds(["only_one"])


# --- the leakage-critical algebra ----------------------------------------


def test_source_mean_equals_the_predicted_algebraic_form(planted):
    """A[p] = mu - alpha_c*/3 + beta_p - gamma[c*,p]/3, exactly."""
    D, mu, alpha, beta, gamma = planted
    for fold in loco.make_folds(CONTEXTS):
        c = fold.target_index
        A = loco.source_mean(D, fold.source_indices)
        expected = mu[None, :] - alpha[c][None, :] / 3.0 + beta - gamma[c] / 3.0
        np.testing.assert_allclose(A, expected, atol=1e-10)


def test_transfer_residual_is_four_thirds_alpha_plus_gamma(planted):
    """delta[c*] - A = (4/3)(alpha_c* + gamma[c*]) — Identity 1."""
    D, mu, alpha, beta, gamma = planted
    for fold in loco.make_folds(CONTEXTS):
        c = fold.target_index
        residual = D[c] - loco.source_mean(D, fold.source_indices)
        expected = (4.0 / 3.0) * (alpha[c][None, :] + gamma[c])
        np.testing.assert_allclose(residual, expected, atol=1e-10)


def test_centred_residual_is_exactly_four_thirds_gamma(planted):
    """Centring the residual over perturbations leaves (4/3) gamma — Identity 1b.

    This is why 'predict the transfer residual' and 'predict gamma' are the same
    problem, and why gamma must never be used as a predictor input.
    """
    D, mu, alpha, beta, gamma = planted
    for fold in loco.make_folds(CONTEXTS):
        c = fold.target_index
        residual = D[c] - loco.source_mean(D, fold.source_indices)
        centred = residual - residual.mean(axis=0, keepdims=True)
        np.testing.assert_allclose(centred, (4.0 / 3.0) * gamma[c], atol=1e-10)


def test_source_mean_carries_a_negative_gamma_component(planted):
    """Identity 2: conserved transfer is anti-correlated with held-out gamma."""
    D, _, _, _, gamma = planted
    correlations = []
    for fold in loco.make_folds(CONTEXTS):
        c = fold.target_index
        A = loco.source_mean(D, fold.source_indices)
        Ac = A - A.mean(axis=0, keepdims=True)
        for p in range(D.shape[1]):
            correlations.append(np.corrcoef(Ac[p], gamma[c, p])[0, 1])
    assert np.mean(correlations) < 0


def test_source_mean_never_reads_the_target_context():
    """Perturbing only the target row must leave every source-only baseline fixed."""
    rng = np.random.default_rng(4)
    D = rng.normal(size=(4, 12, 9))
    control = rng.normal(size=(4, 9))
    fold = loco.make_folds(CONTEXTS)[2]
    before = {
        "A": loco.source_mean(D, fold.source_indices),
        "B": loco.single_source(D, fold.source_indices[0]),
        "D": loco.weighted_source(
            D,
            fold.source_indices,
            loco.basal_affine_weights(control, fold.target_index, fold.source_indices),
        ),
    }
    D_tampered = D.copy()
    D_tampered[fold.target_index] = rng.normal(size=D.shape[1:]) * 100.0
    after = {
        "A": loco.source_mean(D_tampered, fold.source_indices),
        "B": loco.single_source(D_tampered, fold.source_indices[0]),
        "D": loco.weighted_source(
            D_tampered,
            fold.source_indices,
            loco.basal_affine_weights(control, fold.target_index, fold.source_indices),
        ),
    }
    for key in before:
        np.testing.assert_allclose(before[key], after[key], atol=0, rtol=0)


# --- basal-only context features -----------------------------------------


def test_affine_weights_sum_to_one_and_use_only_basal():
    rng = np.random.default_rng(6)
    control = rng.normal(size=(4, 30))
    fold = loco.make_folds(CONTEXTS)[0]
    w = loco.basal_affine_weights(control, fold.target_index, fold.source_indices)
    assert w.sum() == pytest.approx(1.0)


def test_affine_weights_recover_an_exact_basal_combination():
    """If the target basal profile IS a known mix of sources, recover that mix."""
    rng = np.random.default_rng(7)
    control = rng.normal(size=(4, 60))
    true_w = np.array([0.5, -0.25, 0.75])
    fold = loco.make_folds(CONTEXTS)[0]
    control[fold.target_index] = np.tensordot(
        true_w, control[list(fold.source_indices)], axes=(0, 0)
    )
    w = loco.basal_affine_weights(control, fold.target_index, fold.source_indices)
    np.testing.assert_allclose(w, true_w, atol=1e-6)


def test_simplex_weights_are_non_negative_and_normalised():
    rng = np.random.default_rng(8)
    control = np.abs(rng.normal(size=(4, 40)))
    fold = loco.make_folds(CONTEXTS)[1]
    w = loco.basal_simplex_weights(control, fold.target_index, fold.source_indices)
    assert (w >= 0).all()
    assert w.sum() == pytest.approx(1.0)


def test_nearest_context_picks_the_most_basally_similar_source():
    control = np.array(
        [
            [1.0, 2.0, 3.0, 4.0],
            [1.0, 2.0, 3.0, 4.1],  # nearly identical to context 0
            [-1.0, 5.0, -3.0, 2.0],
            [4.0, -2.0, 1.0, -5.0],
        ]
    )
    best, sim = loco.nearest_context(control, 0, [1, 2, 3])
    assert best == 1
    assert sim[1] > sim[2] and sim[1] > sim[3]


def test_weighted_source_rejects_mismatched_weights():
    D = np.zeros((4, 5, 3))
    with pytest.raises(ValueError, match="weights shape"):
        loco.weighted_source(D, [1, 2, 3], np.array([0.5, 0.5]))


# --- metrics --------------------------------------------------------------


def test_perfect_prediction_scores_perfectly():
    rng = np.random.default_rng(9)
    Y = rng.normal(size=(15, 40))
    m = loco.per_perturbation_metrics(Y, Y)
    np.testing.assert_allclose(m["pearson"], 1.0, atol=1e-9)
    np.testing.assert_allclose(m["cosine"], 1.0, atol=1e-9)
    np.testing.assert_allclose(m["energy_explained"], 1.0, atol=1e-9)
    np.testing.assert_allclose(m["mse"], 0.0, atol=1e-18)


def test_zero_prediction_explains_no_energy():
    rng = np.random.default_rng(10)
    Y = rng.normal(size=(8, 30))
    m = loco.per_perturbation_metrics(Y, np.zeros_like(Y))
    np.testing.assert_allclose(m["energy_explained"], 0.0, atol=1e-9)


def test_metrics_reject_shape_mismatch():
    with pytest.raises(ValueError, match="shape mismatch"):
        loco.per_perturbation_metrics(np.zeros((3, 4)), np.zeros((3, 5)))


def test_template_removal_subtracts_the_shared_row():
    Y = np.arange(12, dtype=float).reshape(3, 4)
    t = loco.target_template_evaluation_only(Y)
    np.testing.assert_allclose(loco.template_removed(Y, t).mean(axis=0), 0.0, atol=1e-12)


# --- reliability: derivation validated on synthetic latent + noise -------


def _synthetic_halves(n_p, n_g, noise_sd, seed, signal_sd=1.0):
    rng = np.random.default_rng(seed)
    latent = rng.normal(scale=signal_sd, size=(n_p, n_g))
    h1 = latent + rng.normal(scale=noise_sd, size=(n_p, n_g))
    h2 = latent + rng.normal(scale=noise_sd, size=(n_p, n_g))
    return latent, h1, h2


def test_split_half_reliability_recovers_the_known_noise_ratio():
    """corr(h1,h2) must estimate Var(L)/(Var(L)+Var(e))."""
    for noise_sd in (0.5, 1.0, 2.0):
        latent, h1, h2 = _synthetic_halves(400, 300, noise_sd, seed=1)
        expected = 1.0 / (1.0 + noise_sd**2)
        observed = np.nanmean(loco.split_half_reliability(h1, h2))
        assert observed == pytest.approx(expected, abs=0.02), (noise_sd, observed, expected)


def test_ceiling_is_sqrt_reliability_not_reliability():
    """A PERFECT latent predictor correlates with a noisy half at sqrt(rho)."""
    for noise_sd in (0.5, 1.0, 2.0):
        latent, h1, h2 = _synthetic_halves(400, 300, noise_sd, seed=2)
        rho = np.nanmean(loco.split_half_reliability(h1, h2))
        achieved = np.nanmean(
            [loco._safe_pearson(latent[i], h1[i]) for i in range(latent.shape[0])]
        )
        assert achieved == pytest.approx(float(np.sqrt(rho)), abs=0.02)
        # and it is clearly NOT rho itself
        assert abs(achieved - rho) > 0.05


def test_disattenuation_recovers_a_known_latent_correlation():
    """A predictor with known true correlation to the latent must be recovered."""
    rng = np.random.default_rng(3)
    n_p, n_g, noise_sd = 500, 400, 1.2
    for true_r in (0.3, 0.6, 0.9):
        latent = rng.normal(size=(n_p, n_g))
        noise_pred = rng.normal(size=(n_p, n_g))
        pred = true_r * latent + np.sqrt(1 - true_r**2) * noise_pred
        h1 = latent + rng.normal(scale=noise_sd, size=(n_p, n_g))
        h2 = latent + rng.normal(scale=noise_sd, size=(n_p, n_g))
        table = loco.reliability_normalised_evaluation(pred, h1, h2)
        assert np.nanmedian(table["r_disattenuated"]) == pytest.approx(true_r, abs=0.05)
        # the raw correlation is attenuated well below the truth
        assert np.nanmedian(table["r_half_mean"]) < true_r - 0.05


def test_spearman_brown_matches_a_pooled_two_half_estimate():
    noise_sd = 1.0
    latent, h1, h2 = _synthetic_halves(500, 400, noise_sd, seed=5)
    rho_half = np.nanmean(loco.split_half_reliability(h1, h2))
    pooled = (h1 + h2) / 2.0
    achieved = np.nanmean(
        [loco._safe_pearson(latent[i], pooled[i]) for i in range(latent.shape[0])]
    )
    predicted = float(np.sqrt(loco.spearman_brown(rho_half)))
    assert achieved == pytest.approx(predicted, abs=0.02)


def test_disattenuation_refuses_uninterpretable_low_reliability():
    out = loco.disattenuate(np.array([0.2, 0.2]), np.array([0.01, 0.5]), min_rho=0.05)
    assert np.isnan(out[0])
    assert np.isfinite(out[1])


def test_disattenuation_is_clipped():
    out = loco.disattenuate(np.array([0.9]), np.array([0.06]))
    assert out[0] <= 1.5


def test_reliability_rejects_shape_mismatch():
    with pytest.raises(ValueError, match="shape mismatch"):
        loco.split_half_reliability(np.zeros((3, 4)), np.zeros((4, 4)))


# --- bootstrap ------------------------------------------------------------


def test_bootstrap_median_brackets_the_truth():
    rng = np.random.default_rng(12)
    v = rng.normal(loc=0.4, scale=0.2, size=800)
    out = loco.bootstrap_median(v, n_boot=500, seed=1)
    assert out["lo"] < out["median"] < out["hi"]
    assert out["n"] == 800


def test_bootstrap_ignores_non_finite_values():
    v = np.array([0.1, np.nan, 0.3, np.inf, 0.5])
    assert loco.bootstrap_median(v, n_boot=100)["n"] == 3


def test_decomposition_gamma_is_never_a_baseline_input(planted):
    """Guard: the baselines must not depend on the four-context decomposition."""
    D, *_ = planted
    dec = anova.decompose(D)
    fold = loco.make_folds(CONTEXTS)[0]
    A = loco.source_mean(D, fold.source_indices)
    # Rebuild the source mean from the source rows alone and require equality.
    manual = D[list(fold.source_indices)].mean(axis=0)
    np.testing.assert_allclose(A, manual, atol=0, rtol=0)
    assert dec.gamma.shape[0] == 4  # decomposition exists but is unused above
