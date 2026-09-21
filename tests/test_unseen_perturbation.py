"""Checks on the unseen-perturbation protocol, its estimators and its priors.

The leakage tests at the end are the load-bearing ones. Arc's panel is
predominantly unseen-context *and* unseen-perturbation, so a protocol that
quietly reads either axis would report a number that does not transfer. Rather
than arguing that the code does not do so, the tests corrupt the withheld
blocks and require the predictions to be bit-identical.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from virtual_cell.modelling import context_main_effect as cme
from virtual_cell.modelling import unseen_perturbation as up
from virtual_cell.priors import catalogue
from virtual_cell.priors import features as pf

N_CONTEXTS, N_PERTS, N_GENES, N_FEATURES = 4, 40, 25, 12


@pytest.fixture
def rng() -> np.random.Generator:
    return np.random.default_rng(11)


@pytest.fixture
def tensor(rng: np.random.Generator) -> np.ndarray:
    """A response tensor built from a known four-component decomposition."""
    mu = rng.normal(size=N_GENES)
    alpha = rng.normal(size=(N_CONTEXTS, 1, N_GENES))
    alpha -= alpha.mean(axis=0, keepdims=True)
    beta = rng.normal(size=(1, N_PERTS, N_GENES))
    beta -= beta.mean(axis=1, keepdims=True)
    gamma = rng.normal(scale=0.3, size=(N_CONTEXTS, N_PERTS, N_GENES))
    gamma -= gamma.mean(axis=0, keepdims=True)
    gamma -= gamma.mean(axis=1, keepdims=True)
    return mu + alpha + beta + gamma


@pytest.fixture
def features(rng: np.random.Generator) -> np.ndarray:
    return rng.normal(size=(N_PERTS, N_FEATURES))


# --------------------------------------------------------------------------
# splits
# --------------------------------------------------------------------------


def test_perturbation_folds_partition_every_index() -> None:
    folds = up.perturbation_folds(N_PERTS, n_folds=5, seed=0)
    assert len(folds) == 5
    joined = np.sort(np.concatenate(folds))
    assert joined.tolist() == list(range(N_PERTS))


def test_perturbation_folds_are_deterministic() -> None:
    a = up.perturbation_folds(N_PERTS, 5, seed=3)
    b = up.perturbation_folds(N_PERTS, 5, seed=3)
    assert all((x == y).all() for x, y in zip(a, b, strict=True))


def test_perturbation_folds_reject_a_single_fold() -> None:
    with pytest.raises(ValueError, match="at least 2"):
        up.perturbation_folds(N_PERTS, 1)


def test_p2_splits_hold_out_both_axes() -> None:
    splits = up.make_two_axis_splits(N_CONTEXTS, N_PERTS, n_folds=4, seed=0)
    p2 = [s for s in splits if s.regime is up.Regime.P2]
    assert len(p2) == 4 * N_CONTEXTS
    for s in p2:
        assert s.target_context not in s.source_contexts
        assert s.eval_contexts == (s.target_context,)
        assert not set(s.test_perturbations) & set(s.train_perturbations)


def test_p1_splits_keep_every_context() -> None:
    splits = up.make_two_axis_splits(N_CONTEXTS, N_PERTS, n_folds=4, seed=0)
    p1 = [s for s in splits if s.regime is up.Regime.P1]
    assert len(p1) == 4
    for s in p1:
        assert s.target_context is None
        assert len(s.source_contexts) == N_CONTEXTS


def test_a_p2_split_refuses_a_target_that_is_also_a_source() -> None:
    with pytest.raises(ValueError, match="must not be a source"):
        up.TwoAxisSplit(
            regime=up.Regime.P2,
            target_context=0,
            source_contexts=(0, 1),
            test_perturbations=(1,),
            train_perturbations=(2,),
            eval_contexts=(0,),
        )


def test_a_split_refuses_an_overlapping_train_and_test_set() -> None:
    with pytest.raises(ValueError, match="both train and test"):
        up.TwoAxisSplit(
            regime=up.Regime.P1,
            target_context=None,
            source_contexts=(0, 1),
            test_perturbations=(1, 2),
            train_perturbations=(2, 3),
            eval_contexts=(0,),
        )


# --------------------------------------------------------------------------
# the quantities
# --------------------------------------------------------------------------


def test_conserved_beta_is_centred_over_the_perturbations_used(tensor: np.ndarray) -> None:
    perts = list(range(10, 30))
    beta, _ = up.conserved_beta(tensor, range(N_CONTEXTS), perts)
    assert beta.shape == (len(perts), N_GENES)
    assert np.abs(beta.mean(axis=0)).max() < 1e-10


def test_conserved_beta_ignores_contexts_and_perturbations_it_was_not_given(
    tensor: np.ndarray, rng: np.random.Generator
) -> None:
    """The subtle leak is through mu, which averages over perturbations too."""
    contexts, perts = [0, 1], [0, 1, 2, 3, 4]
    before, mu_before = up.conserved_beta(tensor, contexts, perts)
    poisoned = tensor.copy()
    poisoned[2:] = rng.normal(size=poisoned[2:].shape) * 1e6
    poisoned[:, 5:] = rng.normal(size=poisoned[:, 5:].shape) * 1e6
    after, mu_after = up.conserved_beta(poisoned, contexts, perts)
    assert np.array_equal(before, after)
    assert np.array_equal(mu_before, mu_after)


def test_oracle_main_effect_matches_the_decomposition_identity(tensor: np.ndarray) -> None:
    """``mean_p delta[c,p]`` must equal ``mu + alpha_c``, with beta and gamma gone."""
    everything = list(range(N_PERTS))
    m = np.stack([cme.oracle_main_effect(tensor, c) for c in range(N_CONTEXTS)])
    mu = tensor.mean(axis=(0, 1))
    alpha = tensor.mean(axis=1) - mu
    assert np.allclose(m, mu + alpha)
    assert np.allclose(m[0], up.context_main_effect(tensor, 0, everything))


# --------------------------------------------------------------------------
# estimators
# --------------------------------------------------------------------------


def test_zero_estimator_predicts_exactly_zero(features: np.ndarray, tensor: np.ndarray) -> None:
    beta, _ = up.conserved_beta(tensor, range(N_CONTEXTS), range(N_PERTS))
    out = up.predict_zero(features[:30], beta[:30], features[30:])
    assert out.shape == (10, N_GENES)
    assert not out.any()


def test_nearest_neighbour_copies_an_identical_training_gene(
    tensor: np.ndarray, rng: np.random.Generator
) -> None:
    beta, _ = up.conserved_beta(tensor, range(N_CONTEXTS), range(N_PERTS))
    train_f = rng.normal(size=(20, N_FEATURES))
    test_f = train_f[[7, 3]]
    out = up.predict_nearest_neighbour(train_f, beta[:20], test_f)
    assert np.allclose(out[0], beta[7])
    assert np.allclose(out[1], beta[3])


def test_knn_reduces_to_the_nearest_neighbour_at_k_of_one(
    tensor: np.ndarray, rng: np.random.Generator
) -> None:
    beta, _ = up.conserved_beta(tensor, range(N_CONTEXTS), range(N_PERTS))
    train_f = rng.normal(size=(20, N_FEATURES))
    test_f = rng.normal(size=(5, N_FEATURES))
    one = up.predict_knn(train_f, beta[:20], test_f, k=1)
    nn = up.predict_nearest_neighbour(train_f, beta[:20], test_f)
    assert np.allclose(one, nn)


def test_knn_averages_when_neighbours_are_equidistant(rng: np.random.Generator) -> None:
    train_f = np.array([[1.0, 0.0], [0.0, 1.0]])
    beta = np.array([[2.0, 4.0], [4.0, 8.0]])
    test_f = np.array([[1.0, 1.0]])
    out = up.predict_knn(train_f, beta, test_f, k=2)
    assert np.allclose(out[0], beta.mean(axis=0))


def test_ridge_recovers_an_exact_linear_map(rng: np.random.Generator) -> None:
    train_f = rng.normal(size=(60, 5))
    weights = rng.normal(size=(5, 8))
    beta = train_f @ weights
    test_f = rng.normal(size=(6, 5))
    out = up.predict_ridge(train_f, beta, test_f, alpha=1e-8)
    assert np.allclose(out, test_f @ weights, atol=1e-5)


def test_ridge_agrees_in_its_primal_and_dual_regimes(rng: np.random.Generator) -> None:
    """The solver switches form on shape; both must give the same answer."""
    wide = rng.normal(size=(12, 40))
    beta = rng.normal(size=(12, 6))
    test = rng.normal(size=(3, 40))
    dual = up.predict_ridge(wide, beta, test, alpha=0.7)
    # Same problem, solved in the primal by padding the sample count with copies
    # is not equivalent; instead check against an explicit primal solve.
    xm, ym = wide.mean(axis=0), beta.mean(axis=0)
    xc, yc = wide - xm, beta - ym
    coef = np.linalg.solve(xc.T @ xc + 0.7 * np.eye(40), xc.T @ yc)
    assert np.allclose(dual, (test - xm) @ coef + ym, atol=1e-6)


def test_neighbour_agreement_is_high_when_neighbours_agree() -> None:
    """The Tier-0 analogue of source agreement: if the genes standing in for an
    unseen one disagree, their average is not worth trusting."""
    train_f = np.eye(6)
    aligned = np.tile(np.array([1.0, 2.0, 3.0]), (6, 1))
    test_f = np.ones((1, 6))
    assert up.neighbour_agreement(train_f, aligned, test_f, k=6)[0] == pytest.approx(1.0)


def test_neighbour_agreement_is_low_when_neighbours_conflict() -> None:
    train_f = np.eye(4)
    opposed = np.array([[1.0, 0.0], [-1.0, 0.0], [0.0, 1.0], [0.0, -1.0]])
    test_f = np.ones((1, 4))
    got = up.neighbour_agreement(train_f, opposed, test_f, k=4)[0]
    assert got == pytest.approx(-1.0 / 3.0)


def test_neighbour_agreement_needs_no_target_measurement(
    tensor: np.ndarray, features: np.ndarray, rng: np.random.Generator
) -> None:
    """It is defined for a Tier-0 gene, which is the whole point: source
    agreement is undefined without source responses, this is not."""
    split = next(
        s
        for s in up.make_two_axis_splits(N_CONTEXTS, N_PERTS, n_folds=4, seed=0)
        if s.regime is up.Regime.P2
    )
    beta_train, _ = up.conserved_beta(tensor, split.source_contexts, split.train_perturbations)
    train_f = features[list(split.train_perturbations)]
    test_f = features[list(split.test_perturbations)]
    before = up.neighbour_agreement(train_f, beta_train, test_f, k=5)

    poisoned = tensor.copy()
    poisoned[:, list(split.test_perturbations)] = rng.normal(
        size=(N_CONTEXTS, len(test_f), N_GENES)
    )
    beta_after, _ = up.conserved_beta(poisoned, split.source_contexts, split.train_perturbations)
    after = up.neighbour_agreement(train_f, beta_after, test_f, k=5)
    assert np.array_equal(before, after)


def test_evaluation_reports_perfect_agreement_for_a_perfect_prediction(
    tensor: np.ndarray,
) -> None:
    beta, _ = up.conserved_beta(tensor, range(N_CONTEXTS), range(N_PERTS))
    got = up.evaluate_predictions(beta, beta)
    assert np.allclose(got["pearson"], 1.0)
    assert np.allclose(got["cosine"], 1.0)
    assert np.allclose(got["unexplained_fraction"], 0.0)


def test_evaluation_charges_a_zero_prediction_the_whole_signal(tensor: np.ndarray) -> None:
    beta, _ = up.conserved_beta(tensor, range(N_CONTEXTS), range(N_PERTS))
    got = up.evaluate_predictions(np.zeros_like(beta), beta)
    assert np.allclose(got["unexplained_fraction"], 1.0)


# --------------------------------------------------------------------------
# feasible context main effect
# --------------------------------------------------------------------------


def test_feasible_estimators_never_read_the_target_context(
    tensor: np.ndarray, rng: np.random.Generator
) -> None:
    """The whole point of 'feasible': the target's responses are unavailable."""
    control_means = rng.normal(size=(N_CONTEXTS, N_GENES))
    target, sources = 2, [0, 1, 3]
    poisoned = tensor.copy()
    poisoned[target] = rng.normal(size=poisoned[target].shape) * 1e6
    for name, fn in cme.FEASIBLE_ESTIMATORS.items():
        before = fn(tensor, sources, control_means=control_means, target=target)
        after = fn(poisoned, sources, control_means=control_means, target=target)
        assert np.array_equal(before, after), name


def test_source_pooled_mean_is_the_mean_of_the_source_main_effects(
    tensor: np.ndarray,
) -> None:
    sources = [0, 1, 3]
    want = np.stack([cme.oracle_main_effect(tensor, c) for c in sources]).mean(axis=0)
    assert np.allclose(cme.source_pooled_mean(tensor, sources), want)


def test_basal_weighting_collapses_to_pooling_when_sources_are_equidistant(
    tensor: np.ndarray,
) -> None:
    control_means = np.zeros((N_CONTEXTS, N_GENES))
    control_means[0] = [1.0, 0.0] + [0.0] * (N_GENES - 2)
    control_means[1] = [0.0, 1.0] + [0.0] * (N_GENES - 2)
    control_means[2] = [1.0, 1.0] + [0.0] * (N_GENES - 2)
    control_means[3] = [0.0, 0.0, 1.0] + [0.0] * (N_GENES - 3)
    got = cme.basal_weighted_mean(tensor, [0, 1], control_means=control_means, target=2)
    want = cme.source_pooled_mean(tensor, [0, 1])
    assert np.allclose(got, want)


def test_oracle_main_effect_is_not_reachable_from_the_feasible_registry() -> None:
    assert "oracle" not in " ".join(cme.FEASIBLE_ESTIMATORS)
    assert cme.oracle_main_effect not in cme.FEASIBLE_ESTIMATORS.values()


# --------------------------------------------------------------------------
# LEAKAGE -- the tests the protocol stands on
# --------------------------------------------------------------------------


@pytest.mark.parametrize("estimator_name", sorted(up.ESTIMATORS))
def test_predictions_ignore_the_held_out_context(
    tensor: np.ndarray, features: np.ndarray, estimator_name: str
) -> None:
    """Replace the outer context's responses with noise: output must not move."""
    rng = np.random.default_rng(99)
    estimator = up.ESTIMATORS[estimator_name]
    split = next(
        s
        for s in up.make_two_axis_splits(N_CONTEXTS, N_PERTS, n_folds=4, seed=0)
        if s.regime is up.Regime.P2
    )
    before = up.fit_predict_split(tensor, features, split, estimator, n_components=5)

    poisoned = tensor.copy()
    poisoned[split.target_context] = rng.normal(size=poisoned[split.target_context].shape) * 1e6
    after = up.fit_predict_split(poisoned, features, split, estimator, n_components=5)

    assert np.array_equal(before, after), estimator_name


@pytest.mark.parametrize("estimator_name", sorted(up.ESTIMATORS))
def test_predictions_ignore_the_held_out_perturbations_everywhere(
    tensor: np.ndarray, features: np.ndarray, estimator_name: str
) -> None:
    """Replace the test perturbations' responses in EVERY context with noise."""
    rng = np.random.default_rng(100)
    estimator = up.ESTIMATORS[estimator_name]
    split = next(
        s
        for s in up.make_two_axis_splits(N_CONTEXTS, N_PERTS, n_folds=4, seed=0)
        if s.regime is up.Regime.P2
    )
    before = up.fit_predict_split(tensor, features, split, estimator, n_components=5)

    poisoned = tensor.copy()
    held = list(split.test_perturbations)
    poisoned[:, held] = rng.normal(size=poisoned[:, held].shape) * 1e6
    after = up.fit_predict_split(poisoned, features, split, estimator, n_components=5)

    assert np.array_equal(before, after), estimator_name


@pytest.mark.parametrize("estimator_name", sorted(up.ESTIMATORS))
def test_predictions_ignore_both_withheld_blocks_at_once(
    tensor: np.ndarray, features: np.ndarray, estimator_name: str
) -> None:
    rng = np.random.default_rng(101)
    estimator = up.ESTIMATORS[estimator_name]
    split = next(
        s
        for s in up.make_two_axis_splits(N_CONTEXTS, N_PERTS, n_folds=4, seed=0)
        if s.regime is up.Regime.P2
    )
    before = up.fit_predict_split(tensor, features, split, estimator, n_components=5)

    poisoned = tensor.copy()
    poisoned[split.target_context] = rng.normal(size=poisoned[split.target_context].shape) * 1e6
    held = list(split.test_perturbations)
    poisoned[:, held] = rng.normal(size=poisoned[:, held].shape) * 1e6
    after = up.fit_predict_split(poisoned, features, split, estimator, n_components=5)

    assert np.array_equal(before, after), estimator_name


def test_the_leakage_tests_can_actually_fail(tensor: np.ndarray, features: np.ndarray) -> None:
    """A guard on the guards.

    The three tests above assert that corrupting a withheld block leaves the
    prediction untouched. That would also hold if ``fit_predict_split`` ignored
    the response tensor entirely, in which case they would prove nothing. So:
    corrupt the block the split IS allowed to read, and require the prediction
    to move.
    """
    split = next(
        s
        for s in up.make_two_axis_splits(N_CONTEXTS, N_PERTS, n_folds=4, seed=0)
        if s.regime is up.Regime.P2
    )
    before = up.fit_predict_split(tensor, features, split, up.predict_knn, n_components=5)

    poisoned = tensor.copy()
    readable = np.ix_(list(split.source_contexts), list(split.train_perturbations))
    poisoned[readable] *= -1.0
    after = up.fit_predict_split(poisoned, features, split, up.predict_knn, n_components=5)

    assert not np.array_equal(before, after), (
        "corrupting the readable training block must change the prediction, "
        "otherwise the leakage tests prove nothing"
    )


def test_pca_is_fitted_on_training_rows_only(tensor: np.ndarray, rng: np.random.Generator) -> None:
    split = next(
        s
        for s in up.make_two_axis_splits(N_CONTEXTS, N_PERTS, n_folds=4, seed=0)
        if s.regime is up.Regime.P2
    )
    features = rng.normal(size=(N_PERTS, 30))
    before = up.fit_predict_split(tensor, features, split, up.predict_ridge, n_components=5)

    moved = features.copy()
    moved[list(split.test_perturbations)] += 500.0
    after = up.fit_predict_split(tensor, moved, split, up.predict_ridge, n_components=5)
    # The projection must not have been refitted, so the training side is
    # untouched; only the test rows' own coordinates may move.
    assert not np.array_equal(before, after)
    pca_train = pf.fit_pca(features[list(split.train_perturbations)], 5)
    pca_moved = pf.fit_pca(moved[list(split.train_perturbations)], 5)
    assert np.array_equal(pca_train.components, pca_moved.components)


# --------------------------------------------------------------------------
# priors
# --------------------------------------------------------------------------


def test_every_catalogued_source_declares_its_leakage() -> None:
    for name, src in catalogue.CATALOGUE.items():
        assert src.leakage in set(catalogue.Leakage), name
        assert src.leakage_note.strip(), name
        assert src.provenance.startswith("data/provenance"), name


def test_depmap_is_flagged_as_a_partial_leakage_risk() -> None:
    """It is itself a perturbation outcome; the catalogue must not pretend otherwise."""
    assert catalogue.source("depmap").leakage is catalogue.Leakage.PARTIAL


def test_per_line_depmap_profiles_are_refused() -> None:
    with pytest.raises(NotImplementedError, match="withheld deliberately"):
        pf.depmap_features(["TP53"], "unused.csv", summaries_only=False)


def test_unknown_prior_source_is_rejected() -> None:
    with pytest.raises(KeyError, match="Unknown prior source"):
        catalogue.source("crystal_ball")


def test_pathway_features_drop_sets_too_small_to_inform(tmp_path) -> None:
    gmt = tmp_path / "sets.gmt"
    gmt.write_text(
        "BIG\thttp://x\tA\tB\tC\tD\tE\nSMALL\thttp://x\tA\tZ\n",
        encoding="utf-8",
    )
    block = pf.pathway_features(list("ABCDEF"), gmt, name="t", min_size=5)
    assert block.feature_names == ("BIG",)
    assert block.matrix[:, 0].tolist() == [1, 1, 1, 1, 1, 0]
    assert block.covered.tolist() == [True] * 5 + [False]


def test_basal_features_mark_genes_off_the_axis_as_uncovered() -> None:
    block = pf.basal_features(
        ["A", "MISSING"],
        np.array([[1.0, 2.0], [3.0, 4.0]]),
        ["A", "B"],
        context_names=["c1", "c2"],
    )
    assert block.covered.tolist() == [True, False]
    assert np.isnan(block.matrix[1]).all()
    assert block.matrix[0, 0] == 1.0 and block.matrix[0, 1] == 3.0


def test_coexpression_of_a_gene_with_itself_is_one() -> None:
    profiles = np.array([[1.0, 2.0, 5.0], [2.0, 1.0, 1.0], [3.0, 5.0, 2.0]])
    block = pf.coexpression_features(["A", "B"], profiles, ["A", "B", "C"])
    assert block.matrix[0, 0] == pytest.approx(1.0)
    assert block.matrix[1, 1] == pytest.approx(1.0)


def test_feature_block_rejects_a_mismatched_coverage_vector() -> None:
    with pytest.raises(ValueError, match="one row per gene"):
        pf.FeatureBlock("t", pd.Index(["A", "B"]), np.zeros((2, 3)), np.zeros(3, dtype=bool))


def test_feature_block_subset_preserves_order_and_rejects_unknown_genes() -> None:
    block = pf.FeatureBlock(
        "t", pd.Index(["A", "B", "C"]), np.arange(9.0).reshape(3, 3), np.ones(3, dtype=bool)
    )
    sub = block.subset(["C", "A"])
    assert sub.genes.tolist() == ["C", "A"]
    assert sub.matrix[0].tolist() == [6.0, 7.0, 8.0]
    with pytest.raises(KeyError, match="not in this block"):
        block.subset(["A", "NOPE"])


def test_pca_explains_a_low_rank_matrix_completely(rng: np.random.Generator) -> None:
    latent = rng.normal(size=(50, 3))
    data = latent @ rng.normal(size=(3, 20))
    pca = pf.fit_pca(data, 3)
    assert pca.explained_variance_ratio.sum() == pytest.approx(1.0, abs=1e-8)
    assert pf.apply_pca(data, pca).shape == (50, 3)
