"""The tiered mean-response model: algebra, leakage, and the frozen Tier-0 rule.

These tests pin the two properties the Arc track depends on. First, that
``beta_hat`` carries only the perturbation-specific part, so that adding
``m_hat`` does not count the template twice. Second, that nothing in the
shrinkage selection can read the held-out context — the guarantee that makes a
weight chosen on public data usable on Arc.
"""

from __future__ import annotations

import numpy as np
import pytest

from virtual_cell.modelling import context_main_effect as cme
from virtual_cell.modelling import mean_response as mr


def planted(
    n_contexts: int = 4,
    n_pert: int = 40,
    n_genes: int = 25,
    seed: int = 0,
    gamma_scale: float = 0.3,
) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    """A balanced ``delta = mu + alpha + beta + gamma`` tensor with known parts."""
    rng = np.random.default_rng(seed)
    mu = rng.normal(size=n_genes)
    alpha = rng.normal(size=(n_contexts, n_genes))
    alpha -= alpha.mean(axis=0, keepdims=True)
    beta = rng.normal(size=(n_pert, n_genes))
    beta -= beta.mean(axis=0, keepdims=True)
    gamma = gamma_scale * rng.normal(size=(n_contexts, n_pert, n_genes))
    gamma -= gamma.mean(axis=0, keepdims=True)
    gamma -= gamma.mean(axis=1, keepdims=True)
    delta = mu[None, None, :] + alpha[:, None, :] + beta[None, :, :] + gamma
    control_means = rng.normal(size=(n_contexts, n_genes))
    return delta, control_means, {"mu": mu, "alpha": alpha, "beta": beta, "gamma": gamma}


# --------------------------------------------------------------------------
# beta: the perturbation-specific part, and only that
# --------------------------------------------------------------------------


def test_centred_beta_has_zero_mean_over_perturbations():
    delta, _, _ = planted()
    beta = mr.centred_source_beta(delta, [0, 1, 2])
    assert np.allclose(beta.mean(axis=0), 0.0, atol=1e-12)


def test_centred_beta_recovers_planted_beta_when_there_is_no_interaction():
    delta, _, parts = planted(gamma_scale=0.0)
    assert np.allclose(mr.centred_source_beta(delta, [0, 1, 2]), parts["beta"])


def test_centred_beta_carries_the_source_mean_interaction():
    """The documented bias: it estimates ``beta + mean_S gamma``, exactly."""
    delta, _, parts = planted()
    sources = [0, 1, 2]
    expected = parts["beta"] + parts["gamma"][sources].mean(axis=0)
    assert np.allclose(mr.centred_source_beta(delta, sources), expected)


def test_centring_happens_before_averaging():
    """A source with a huge template must not leak it into beta_hat."""
    delta, _, _ = planted()
    inflated = delta.copy()
    inflated[1] += 50.0  # a pure template shift in one source
    assert np.allclose(
        mr.centred_source_beta(delta, [0, 1, 2]),
        mr.centred_source_beta(inflated, [0, 1, 2]),
    )


def test_tier_zero_beta_is_exactly_zero():
    delta, _, _ = planted()
    beta = mr.tier_beta(delta, [0, 1, 2], 0)
    assert beta.shape == (delta.shape[1], delta.shape[2])
    assert np.array_equal(beta, np.zeros_like(beta))


def test_tier_zero_ignores_its_sources_entirely():
    delta, _, _ = planted()
    other, _, _ = planted(seed=99)
    assert np.array_equal(mr.tier_beta(delta, [0, 1], 0), mr.tier_beta(other, [2, 3], 0))


def test_unknown_tier_is_rejected():
    delta, _, _ = planted()
    with pytest.raises(ValueError, match="unknown tier"):
        mr.tier_beta(delta, [0, 1], 3)


# --------------------------------------------------------------------------
# source subsets
# --------------------------------------------------------------------------


def test_source_subsets_are_all_combinations_in_frozen_order():
    assert mr.source_subsets([2, 0, 1], 2) == [(0, 1), (0, 2), (1, 2)]
    assert mr.source_subsets([2, 0, 1], 1) == [(0,), (1,), (2,)]


def test_source_subsets_rejects_impossible_sizes():
    with pytest.raises(ValueError, match="cannot take subsets"):
        mr.source_subsets([0, 1], 3)


def test_tier_source_counts_match_the_frozen_policy():
    assert mr.TIER_SOURCE_COUNT == {2: 2, 1: 1, 0: None}


# --------------------------------------------------------------------------
# the model
# --------------------------------------------------------------------------


def test_predict_is_the_stated_formula():
    delta, _, _ = planted()
    m = np.arange(delta.shape[2], dtype=float)
    beta = mr.centred_source_beta(delta, [0, 1])
    got = mr.predict(m, beta, 0.25)
    assert np.allclose(got, m[None, :] + 0.25 * beta)


def test_weight_zero_collapses_to_the_main_effect():
    delta, _, _ = planted()
    m = np.linspace(-1, 1, delta.shape[2])
    got = mr.predict(m, mr.centred_source_beta(delta, [0, 1]), 0.0)
    assert np.allclose(got, np.broadcast_to(m, got.shape))


def test_prediction_mean_over_perturbations_is_the_main_effect():
    """Centring makes the model's own panel mean equal ``m_hat``, for any w."""
    delta, _, _ = planted()
    m = np.linspace(-1, 1, delta.shape[2])
    for weight in mr.SHRINKAGE_GRID:
        pred = mr.predict(m, mr.centred_source_beta(delta, [0, 1, 2]), weight)
        assert np.allclose(pred.mean(axis=0), m, atol=1e-12)


def test_predict_rejects_a_mismatched_gene_axis():
    delta, _, _ = planted()
    with pytest.raises(ValueError, match="is not"):
        mr.predict(np.zeros(3), mr.centred_source_beta(delta, [0, 1]), 1.0)


# --------------------------------------------------------------------------
# m_hat estimators
# --------------------------------------------------------------------------


def test_m0_is_zero_and_reads_nothing():
    delta, control_means, _ = planted()
    out = mr.zero_main_effect(delta, [0, 1, 2], control_means=control_means, target=3)
    assert np.array_equal(out, np.zeros(delta.shape[2]))


def test_every_estimator_returns_one_value_per_gene():
    delta, control_means, _ = planted()
    for name, fn in mr.MAIN_EFFECT_ESTIMATORS.items():
        out = fn(delta, [0, 1, 2], control_means=control_means, target=3)
        assert out.shape == (delta.shape[2],), name


def test_shrunk_estimator_is_a_scalar_multiple_of_its_base():
    delta, control_means, _ = planted()
    base = cme.source_pooled_mean(delta, [0, 1, 2], control_means=control_means)
    shrunk = mr.source_pooled_shrunk_mean(delta, [0, 1, 2], control_means=control_means, target=3)
    ratio = shrunk[np.abs(base) > 1e-9] / base[np.abs(base) > 1e-9]
    assert np.allclose(ratio, ratio[0])


def test_main_effect_scale_is_one_when_sources_agree_perfectly():
    """No interaction and no context term: the estimate needs no shrinking."""
    delta, control_means, _ = planted(gamma_scale=0.0)
    delta = delta - delta.mean(axis=0, keepdims=True) + delta.mean(axis=0, keepdims=True)
    identical = np.broadcast_to(delta[0], delta.shape).copy()
    scale = mr.fit_main_effect_scale(identical, [0, 1, 2], cme.source_pooled_mean)
    assert scale == pytest.approx(1.0, abs=1e-9)


def test_no_estimator_reads_the_target_context_response():
    """Replace the target row with noise; every feasible m_hat must be unchanged."""
    delta, control_means, _ = planted()
    rng = np.random.default_rng(7)
    poisoned = delta.copy()
    poisoned[3] = rng.normal(scale=100.0, size=delta.shape[1:])
    for name, fn in mr.MAIN_EFFECT_ESTIMATORS.items():
        a = fn(delta, [0, 1, 2], control_means=control_means, target=3)
        b = fn(poisoned, [0, 1, 2], control_means=control_means, target=3)
        assert np.array_equal(a, b), name


# --------------------------------------------------------------------------
# nested shrinkage selection
# --------------------------------------------------------------------------


def test_shrinkage_selection_never_reads_the_target():
    delta, control_means, _ = planted()
    rng = np.random.default_rng(11)
    poisoned = delta.copy()
    poisoned[3] = rng.normal(scale=100.0, size=delta.shape[1:])
    kwargs = dict(
        control_means=control_means,
        main_effect=mr.MAIN_EFFECT_ESTIMATORS["M3b_basal_shrunk"],
    )
    a = mr.select_tier_shrinkage(delta, [0, 1, 2], **kwargs)
    b = mr.select_tier_shrinkage(poisoned, [0, 1, 2], **kwargs)
    assert a.weights == b.weights
    assert a.grid == b.grid


def test_selected_weights_lie_on_the_grid_and_tier_zero_is_pinned():
    delta, control_means, _ = planted()
    selected = mr.select_tier_shrinkage(
        delta,
        [0, 1, 2],
        control_means=control_means,
        main_effect=mr.MAIN_EFFECT_ESTIMATORS["M1_source_pooled"],
    )
    assert selected[0] == 0.0
    for tier in (1, 2):
        assert selected[tier] in mr.SHRINKAGE_GRID


def test_no_interaction_selects_no_shrinkage():
    """With gamma = 0 the source estimate is unbiased, so ``w = 1`` must win."""
    delta, control_means, _ = planted(gamma_scale=0.0, n_pert=60)
    selected = mr.select_tier_shrinkage(
        delta,
        [0, 1, 2],
        control_means=control_means,
        main_effect=mr.MAIN_EFFECT_ESTIMATORS["M1_source_pooled"],
    )
    assert selected[2] == 1.0
    assert selected[1] == 1.0


def test_overwhelming_interaction_selects_full_shrinkage():
    """With gamma dominating, transferring beta at all must lose."""
    delta, control_means, _ = planted(gamma_scale=12.0, n_pert=60)
    selected = mr.select_tier_shrinkage(
        delta,
        [0, 1, 2],
        control_means=control_means,
        main_effect=mr.MAIN_EFFECT_ESTIMATORS["M1_source_pooled"],
    )
    assert selected[2] == 0.0


def test_tier_one_is_shrunk_at_least_as_hard_as_tier_two():
    """One source carries more interaction than two, so it cannot deserve less."""
    delta, control_means, _ = planted(gamma_scale=1.0, n_pert=80, seed=5)
    selected = mr.select_tier_shrinkage(
        delta,
        [0, 1, 2],
        control_means=control_means,
        main_effect=mr.MAIN_EFFECT_ESTIMATORS["M1_source_pooled"],
    )
    assert selected[1] <= selected[2]


def test_selection_needs_at_least_two_sources():
    delta, control_means, _ = planted()
    with pytest.raises(ValueError, match="at least two source"):
        mr.select_tier_shrinkage(
            delta,
            [0],
            control_means=control_means,
            main_effect=mr.MAIN_EFFECT_ESTIMATORS["M1_source_pooled"],
        )


def test_grid_is_the_frozen_conservative_one():
    assert mr.SHRINKAGE_GRID == (0.0, 0.25, 0.5, 0.75, 1.0)
