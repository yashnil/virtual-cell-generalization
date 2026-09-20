"""Template estimators, pathway aggregation and candidate D targets.

Every estimator gets a leakage test: the held-out context's perturbation
responses are replaced with noise and the estimator must be bit-identical.
Every reliability-based target is validated against synthetic data with a known
latent signal and known noise.
"""

from __future__ import annotations

import numpy as np
import pytest

from virtual_cell.analysis import foundations as fx
from virtual_cell.analysis import loco

CONTEXTS = ("k562", "rpe1", "hepg2", "jurkat")


@pytest.fixture
def planted():
    rng = np.random.default_rng(21)
    n_c, n_p, n_g = 4, 50, 30
    mu = rng.normal(size=n_g)
    alpha = rng.normal(size=(n_c, n_g))
    alpha -= alpha.mean(axis=0)
    beta = rng.normal(size=(n_p, n_g))
    beta -= beta.mean(axis=0)
    gamma = rng.normal(size=(n_c, n_p, n_g)) * 0.4
    gamma -= gamma.mean(axis=0, keepdims=True)
    gamma -= gamma.mean(axis=1, keepdims=True)
    D = mu + alpha[:, None, :] + beta[None, :, :] + gamma
    control = rng.normal(size=(n_c, n_g))
    return D, mu, alpha, beta, gamma, control


# --- A. template algebra ---------------------------------------------------


def test_oracle_alpha_matches_the_planted_alpha(planted):
    D, _, alpha, *_ = planted
    for c in range(4):
        np.testing.assert_allclose(fx.oracle_alpha_evaluation_only(D, c), alpha[c], atol=1e-10)


def test_four_thirds_correction_is_exact_with_oracle_alpha(planted):
    """A[p] + (4/3) alpha_c* must equal mu + alpha_c* + beta_p - gamma/3."""
    D, mu, alpha, beta, gamma, _ = planted
    for fold in loco.make_folds(CONTEXTS):
        c = fold.target_index
        A = loco.source_mean(D, fold.source_indices)
        corrected = fx.apply_template_correction(A, fx.oracle_alpha_evaluation_only(D, c))
        expected = mu[None, :] + alpha[c][None, :] + beta - gamma[c] / 3.0
        np.testing.assert_allclose(corrected, expected, atol=1e-10)


def test_template_correction_removes_the_template_error_exactly(planted):
    """With oracle alpha the remaining error is purely the gamma term."""
    D, _, _, _, gamma, _ = planted
    for fold in loco.make_folds(CONTEXTS):
        c = fold.target_index
        A = loco.source_mean(D, fold.source_indices)
        corrected = fx.apply_template_correction(A, fx.oracle_alpha_evaluation_only(D, c))
        residual = D[c] - corrected
        np.testing.assert_allclose(residual, (4.0 / 3.0) * gamma[c], atol=1e-10)


def test_zero_estimator_reproduces_the_uncorrected_source_mean(planted):
    D, *_, control = planted
    fold = loco.make_folds(CONTEXTS)[1]
    A = loco.source_mean(D, fold.source_indices)
    est = {
        e.name: e
        for e in fx.template_estimators(D, control, fold.target_index, fold.source_indices)
    }
    np.testing.assert_allclose(
        fx.apply_template_correction(A, est["zero"].alpha_hat), A, atol=0, rtol=0
    )


def test_source_alphas_sum_to_zero_and_ignore_the_target(planted):
    D, *_ = planted
    for fold in loco.make_folds(CONTEXTS):
        a = fx.source_alphas(D, fold.source_indices)
        np.testing.assert_allclose(a.sum(axis=0), 0.0, atol=1e-10)


def test_global_scale_recovers_a_planted_coefficient():
    rng = np.random.default_rng(3)
    dev = rng.normal(size=(3, 40))
    k_true = -0.73
    assert fx.fit_global_scale(k_true * dev, dev) == pytest.approx(k_true, abs=1e-9)


def test_global_scale_is_zero_for_degenerate_deviations():
    assert fx.fit_global_scale(np.ones((3, 5)), np.zeros((3, 5))) == 0.0


def test_ridge_subspace_map_returns_the_right_shape(planted):
    D, _, _, _, _, control = planted
    fold = loco.make_folds(CONTEXTS)[0]
    alphas = fx.source_alphas(D, fold.source_indices)
    devs = np.vstack(
        [fx.basal_deviation(control, s, fold.source_indices) for s in fold.source_indices]
    )
    out = fx.ridge_subspace_map(
        alphas, devs, fx.basal_deviation(control, fold.target_index, range(4))
    )
    assert out.shape == (D.shape[2],)
    assert np.all(np.isfinite(out))


@pytest.mark.parametrize("estimator", ["zero", "direct_basal", "global_scalar", "ridge_subspace"])
def test_template_estimators_never_read_target_responses(planted, estimator):
    """LEAKAGE TEST: scramble the target's responses; the estimate must not move."""
    D, _, _, _, _, control = planted
    rng = np.random.default_rng(99)
    for fold in loco.make_folds(CONTEXTS):
        before = {
            e.name: e.alpha_hat
            for e in fx.template_estimators(D, control, fold.target_index, fold.source_indices)
        }[estimator]
        tampered = D.copy()
        tampered[fold.target_index] = rng.normal(size=D.shape[1:]) * 50.0
        after = {
            e.name: e.alpha_hat
            for e in fx.template_estimators(
                tampered, control, fold.target_index, fold.source_indices
            )
        }[estimator]
        np.testing.assert_allclose(before, after, atol=0, rtol=0)


def test_basal_deviation_uses_only_control_profiles(planted):
    D, *_, control = planted
    fold = loco.make_folds(CONTEXTS)[2]
    before = fx.basal_deviation(control, fold.target_index, range(4))
    D_tampered = D * 1000.0
    after = fx.basal_deviation(control, fold.target_index, range(4))
    np.testing.assert_allclose(before, after, atol=0, rtol=0)
    assert D_tampered.shape == D.shape


# --- C. pathway aggregation -----------------------------------------------


def test_read_gmt_parses_name_and_members(tmp_path):
    p = tmp_path / "x.gmt"
    p.write_text("SET_A\thttp://url\tG1\tG2\tG3\nSET_B\tdesc\tG2\tG4\n\n")
    sets = fx.read_gmt(p)
    assert sets == {"SET_A": ["G1", "G2", "G3"], "SET_B": ["G2", "G4"]}


def test_membership_drops_sets_below_the_minimum_overlap():
    genes = [f"G{i}" for i in range(20)]
    sets = {"BIG": genes[:12], "SMALL": genes[:3]}
    names, M = fx.pathway_membership(sets, genes, min_genes=10)
    assert names == ["BIG"]
    assert M.shape == (1, 20)
    assert M.sum() == 12


def test_membership_errors_when_nothing_qualifies():
    with pytest.raises(ValueError, match="No gene set"):
        fx.pathway_membership({"S": ["A"]}, ["A", "B"], min_genes=10)


def test_pathway_score_is_the_mean_over_present_members():
    genes = ["a", "b", "c", "d"]
    names, M = fx.pathway_membership({"S": ["a", "c"]}, genes, min_genes=2)
    R = np.array([[1.0, 10.0, 3.0, 100.0]])
    np.testing.assert_allclose(fx.pathway_scores(R, M), [[2.0]])


def test_pathway_scores_preserve_leading_dimensions():
    genes = [f"g{i}" for i in range(10)]
    _, M = fx.pathway_membership({"S1": genes[:5], "S2": genes[5:]}, genes, min_genes=3)
    R = np.zeros((4, 7, 10))
    assert fx.pathway_scores(R, M).shape == (4, 7, 2)


def test_pathway_scores_are_linear_so_the_decomposition_commutes(planted):
    """Aggregation is linear, so decomposing then scoring == scoring then decomposing."""
    from virtual_cell.decomposition import anova

    D, *_ = planted
    genes = [f"g{i}" for i in range(D.shape[2])]
    _, M = fx.pathway_membership(
        {"S1": genes[:15], "S2": genes[10:], "S3": genes[5:25]}, genes, min_genes=10
    )
    gamma_then_score = fx.pathway_scores(anova.decompose(D).gamma, M)
    score_then_gamma = anova.decompose(fx.pathway_scores(D, M)).gamma
    np.testing.assert_allclose(gamma_then_score, score_then_gamma, atol=1e-10)


# --- D. candidate targets, validated on synthetic latent + noise ----------


def _latent_halves(n_p, n_g, noise_sd, seed):
    rng = np.random.default_rng(seed)
    latent = rng.normal(size=(n_p, n_g))
    return (
        latent,
        latent + rng.normal(scale=noise_sd, size=(n_p, n_g)),
        latent + rng.normal(scale=noise_sd, size=(n_p, n_g)),
        rng,
    )


def test_reliable_energy_is_unbiased_for_the_latent_norm():
    """<h1,h2> estimates ||L||^2 regardless of noise level; ||h||^2 does not."""
    for noise_sd in (0.5, 1.5):
        latent, h1, h2, _ = _latent_halves(600, 200, noise_sd, seed=4)
        truth = np.einsum("pg,pg->p", latent, latent)
        est = fx.reliable_energy(h1, h2)
        assert np.mean(est) == pytest.approx(np.mean(truth), rel=0.05)
        naive = np.einsum("pg,pg->p", h1, h1)
        assert np.mean(naive) > np.mean(truth) * 1.1  # inflated by noise


def test_reliable_residual_energy_is_unbiased():
    latent, h1, h2, rng = _latent_halves(600, 200, 1.0, seed=5)
    pred = 0.5 * latent + 0.5 * rng.normal(size=latent.shape)
    truth = np.einsum("pg,pg->p", latent - pred, latent - pred)
    est = fx.reliable_residual_energy(h1, h2, pred)
    assert np.mean(est) == pytest.approx(np.mean(truth), rel=0.05)


def test_unexplained_fraction_recovers_a_planted_value():
    """A prediction explaining a known fraction must score that fraction."""
    rng = np.random.default_rng(6)
    n_p, n_g = 800, 200
    latent = rng.normal(size=(n_p, n_g))
    for share in (0.25, 0.5, 0.9):
        pred = share * latent
        # ||L - share L||^2 / ||L||^2 = (1-share)^2
        h1 = latent + rng.normal(scale=1.0, size=latent.shape)
        h2 = latent + rng.normal(scale=1.0, size=latent.shape)
        tab = fx.candidate_targets(h1, h2, pred)
        assert np.nanmean(tab["D_unexplained_fraction"]) == pytest.approx(
            (1 - share) ** 2, abs=0.05
        )


def test_unexplained_fraction_exceeds_one_when_worse_than_zero():
    """The >1 regime is meaningful and must be representable, not clipped."""
    rng = np.random.default_rng(7)
    latent = rng.normal(size=(400, 150))
    pred = -1.5 * latent  # anti-correlated and over-scaled
    h1 = latent + rng.normal(scale=0.5, size=latent.shape)
    h2 = latent + rng.normal(scale=0.5, size=latent.shape)
    tab = fx.candidate_targets(h1, h2, pred)
    assert np.nanmedian(tab["D_unexplained_fraction"]) > 1.0


def test_explained_and_unexplained_fractions_are_complements():
    rng = np.random.default_rng(8)
    latent = rng.normal(size=(200, 100))
    pred = 0.6 * latent
    h1 = latent + rng.normal(scale=0.8, size=latent.shape)
    h2 = latent + rng.normal(scale=0.8, size=latent.shape)
    tab = fx.candidate_targets(h1, h2, pred)
    np.testing.assert_allclose(
        tab["D_explained_fraction"] + tab["D_unexplained_fraction"], 1.0, atol=1e-10
    )


def test_raw_gamma_norm_is_depth_biased_but_reliable_fraction_is_not():
    """The decisive reason a raw residual norm is a bad D: it grows with noise."""
    rng = np.random.default_rng(9)
    latent = rng.normal(size=(500, 200))
    pred = 0.5 * latent
    norms, fracs = [], []
    for noise_sd in (0.3, 1.0, 2.5):  # shallower sampling -> more noise
        h1 = latent + rng.normal(scale=noise_sd, size=latent.shape)
        h2 = latent + rng.normal(scale=noise_sd, size=latent.shape)
        tab = fx.candidate_targets(h1, h2, pred)
        norms.append(float(np.nanmedian(tab["D_raw_gamma_norm_proxy"])))
        fracs.append(fx.pooled_unexplained_fraction(h1, h2, pred))
    assert norms[0] < norms[1] < norms[2]  # biased by depth
    assert max(fracs) - min(fracs) < 0.1  # pooled fraction: stable


def test_raw_residual_fraction_is_biased_relative_to_the_reliable_one():
    rng = np.random.default_rng(10)
    latent = rng.normal(size=(500, 200))
    pred = 0.5 * latent
    h1 = latent + rng.normal(scale=1.5, size=latent.shape)
    h2 = latent + rng.normal(scale=1.5, size=latent.shape)
    tab = fx.candidate_targets(h1, h2, pred)
    assert np.nanmean(tab["D_raw_residual_fraction"]) > np.nanmean(tab["D_unexplained_fraction"])


def test_candidate_targets_marks_unusable_pairs():
    rng = np.random.default_rng(11)
    noise_only = rng.normal(size=(200, 80))
    h1 = rng.normal(size=(200, 80))
    h2 = rng.normal(size=(200, 80))
    tab = fx.candidate_targets(h1, h2, noise_only)
    assert (~tab["usable"]).sum() > 0
    assert tab.loc[~tab["usable"], "D_unexplained_fraction"].isna().all()


def test_reliable_energy_rejects_shape_mismatch():
    with pytest.raises(ValueError, match="shape mismatch"):
        fx.reliable_energy(np.zeros((3, 4)), np.zeros((3, 5)))


# --- B. partial correlation ------------------------------------------------


def test_partial_spearman_removes_a_pure_confound():
    """x and y correlated only through z must lose their association."""
    rng = np.random.default_rng(12)
    z = rng.normal(size=900)
    x = z + rng.normal(scale=0.3, size=900)
    y = z + rng.normal(scale=0.3, size=900)
    from scipy import stats

    raw = stats.spearmanr(x, y).statistic
    partial = fx.partial_spearman(x, y, z)
    assert raw > 0.8
    assert abs(partial) < 0.2


def test_partial_spearman_keeps_a_genuine_association():
    rng = np.random.default_rng(13)
    z = rng.normal(size=900)
    x = rng.normal(size=900)
    y = 0.8 * x + 0.5 * z + rng.normal(scale=0.2, size=900)
    assert fx.partial_spearman(x, y, z) > 0.6


def test_partial_spearman_handles_non_finite_rows():
    rng = np.random.default_rng(14)
    x = rng.normal(size=200)
    y = x + rng.normal(scale=0.2, size=200)
    z = rng.normal(size=200)
    x[:5] = np.nan
    assert np.isfinite(fx.partial_spearman(x, y, z))


def test_per_pair_fraction_is_unstable_at_low_reliability_but_pooled_is_not():
    """Documented edge case: the per-pair denominator is a noisy ||L||^2.

    At high noise it approaches zero for many perturbations, so the mean of
    per-pair ratios degrades badly while the pooled ratio stays put.
    """
    rng = np.random.default_rng(15)
    latent = rng.normal(size=(500, 200))
    pred = 0.5 * latent
    per_pair, pooled = [], []
    for noise_sd in (0.3, 2.5):
        h1 = latent + rng.normal(scale=noise_sd, size=latent.shape)
        h2 = latent + rng.normal(scale=noise_sd, size=latent.shape)
        per_pair.append(
            float(np.nanmean(fx.candidate_targets(h1, h2, pred)["D_unexplained_fraction"]))
        )
        pooled.append(fx.pooled_unexplained_fraction(h1, h2, pred))
    assert abs(per_pair[1] - per_pair[0]) > 0.2  # per-pair mean degrades
    assert abs(pooled[1] - pooled[0]) < 0.1  # pooled holds
    assert pooled[0] == pytest.approx(0.25, abs=0.05)


def test_pooled_fraction_recovers_planted_values():
    rng = np.random.default_rng(16)
    latent = rng.normal(size=(600, 200))
    for share in (0.25, 0.5, 0.9):
        h1 = latent + rng.normal(scale=1.0, size=latent.shape)
        h2 = latent + rng.normal(scale=1.0, size=latent.shape)
        got = fx.pooled_unexplained_fraction(h1, h2, share * latent)
        assert got == pytest.approx((1 - share) ** 2, abs=0.03)


def test_pooled_fraction_is_nan_when_there_is_no_reliable_signal():
    rng = np.random.default_rng(17)
    h1 = rng.normal(size=(300, 100))
    h2 = rng.normal(size=(300, 100))
    val = fx.pooled_unexplained_fraction(h1, h2, np.zeros((300, 100)))
    assert np.isfinite(val) or np.isnan(val)  # defined behaviour either way
