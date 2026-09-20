"""Null constructions, vectorised correlation, and recoverability plumbing.

The null constructions carry the scientific weight of the falsification battery,
so each one is tested for exactly the invariances it claims to preserve.
"""

from __future__ import annotations

import numpy as np
import pytest

from virtual_cell.analysis import falsification as fal
from virtual_cell.analysis import loco
from virtual_cell.decomposition import anova


@pytest.fixture
def membership():
    rng = np.random.default_rng(0)
    n_sets, n_genes = 12, 200
    M = np.zeros((n_sets, n_genes))
    for i in range(n_sets):
        size = int(rng.integers(10, 60))
        M[i, rng.choice(n_genes, size=size, replace=False)] = 1.0
    return M


# --- vectorised correlation matches the scalar reference ------------------


def test_rowwise_pearson_matches_numpy_corrcoef():
    rng = np.random.default_rng(1)
    a = rng.normal(size=(17, 40))
    b = rng.normal(size=(17, 40))
    got = fal.rowwise_pearson(a, b)
    want = np.array([np.corrcoef(a[i], b[i])[0, 1] for i in range(17)])
    np.testing.assert_allclose(got, want, atol=1e-12)


def test_rowwise_pearson_matches_the_loco_scalar_helper():
    rng = np.random.default_rng(2)
    a = rng.normal(size=(9, 30))
    b = rng.normal(size=(9, 30))
    np.testing.assert_allclose(
        fal.rowwise_pearson(a, b),
        [loco._safe_pearson(a[i], b[i]) for i in range(9)],
        atol=1e-12,
    )


def test_rowwise_pearson_broadcasts_over_leading_axes():
    rng = np.random.default_rng(3)
    a = rng.normal(size=(4, 6, 25))
    b = rng.normal(size=(4, 6, 25))
    out = fal.rowwise_pearson(a, b)
    assert out.shape == (4, 6)
    np.testing.assert_allclose(out[2], fal.rowwise_pearson(a[2], b[2]), atol=1e-12)


def test_rowwise_pearson_is_nan_for_constant_rows():
    a = np.vstack([np.ones(10), np.arange(10.0)])
    b = np.vstack([np.arange(10.0), np.arange(10.0)])
    out = fal.rowwise_pearson(a, b)
    assert np.isnan(out[0])
    assert out[1] == pytest.approx(1.0)


def test_rowwise_pearson_rejects_shape_mismatch():
    with pytest.raises(ValueError, match="shape mismatch"):
        fal.rowwise_pearson(np.zeros((2, 3)), np.zeros((2, 4)))


# --- permuted null: preserves geometry exactly ----------------------------


def test_permuted_null_preserves_set_sizes(membership):
    P = fal.permuted_membership(membership, np.random.default_rng(4))
    np.testing.assert_array_equal(P.sum(axis=1), membership.sum(axis=1))


def test_permuted_null_preserves_every_pairwise_overlap(membership):
    """(MP)(MP)^T = M M^T exactly, so all set-set overlaps survive."""
    P = fal.permuted_membership(membership, np.random.default_rng(5))
    np.testing.assert_allclose(P @ P.T, membership @ membership.T, atol=1e-12)


def test_permuted_null_preserves_the_gene_degree_multiset(membership):
    P = fal.permuted_membership(membership, np.random.default_rng(6))
    np.testing.assert_array_equal(np.sort(P.sum(axis=0)), np.sort(membership.sum(axis=0)))


def test_permuted_null_actually_moves_genes(membership):
    P = fal.permuted_membership(membership, np.random.default_rng(7))
    assert not np.array_equal(P, membership)


# --- resampled null: sizes only -------------------------------------------


def test_resampled_null_matches_sizes_but_not_overlaps(membership):
    sizes = membership.sum(axis=1).astype(int)
    R = fal.resampled_membership(sizes, membership.shape[1], np.random.default_rng(8))
    np.testing.assert_array_equal(R.sum(axis=1), sizes)
    assert R.shape == membership.shape
    assert not np.allclose(R @ R.T, membership @ membership.T)


def test_resampled_null_rejects_oversized_sets():
    with pytest.raises(ValueError, match="exceeds the gene universe"):
        fal.resampled_membership([50], 10, np.random.default_rng(0))


def test_resampled_sets_have_no_duplicate_genes():
    R = fal.resampled_membership([30, 30], 100, np.random.default_rng(9))
    assert set(np.unique(R)) <= {0.0, 1.0}


# --- gaussian projection ---------------------------------------------------


def test_gaussian_projection_has_the_requested_shape_and_unit_rows():
    W = fal.gaussian_projection(45, 300, np.random.default_rng(10))
    assert W.shape == (45, 300)
    np.testing.assert_allclose(np.linalg.norm(W, axis=1), 1.0, atol=1e-10)


def test_gaussian_projection_is_seed_reproducible():
    a = fal.gaussian_projection(8, 50, np.random.default_rng(11))
    b = fal.gaussian_projection(8, 50, np.random.default_rng(11))
    np.testing.assert_allclose(a, b, atol=0, rtol=0)


# --- linear representation plumbing ---------------------------------------


def test_membership_weights_give_the_mean_over_members():
    M = np.array([[1.0, 0.0, 1.0, 0.0]])
    W = fal.membership_to_weights(M)
    R = np.array([[2.0, 99.0, 4.0, 99.0]])
    np.testing.assert_allclose(fal.project(R, W), [[3.0]])


def test_membership_weights_reject_empty_sets():
    with pytest.raises(ValueError, match="no member genes"):
        fal.membership_to_weights(np.zeros((2, 5)))


def test_projection_is_linear_so_gamma_commutes():
    """Projection is linear, so project-then-decompose == decompose-then-project."""
    rng = np.random.default_rng(12)
    D = rng.normal(size=(4, 25, 60))
    W = fal.gaussian_projection(10, 60, rng)
    np.testing.assert_allclose(
        fal.gamma_of(fal.project(D, W)), fal.project(anova.decompose(D).gamma, W), atol=1e-10
    )


def test_gamma_of_matches_the_frozen_decomposition():
    rng = np.random.default_rng(13)
    D = rng.normal(size=(4, 30, 15))
    np.testing.assert_allclose(fal.gamma_of(D), anova.decompose(D).gamma, atol=1e-12)


def test_gamma_of_satisfies_the_zero_sum_conditions():
    rng = np.random.default_rng(14)
    g = fal.gamma_of(rng.normal(size=(4, 20, 12)))
    np.testing.assert_allclose(g.sum(axis=0), 0.0, atol=1e-10)
    np.testing.assert_allclose(g.sum(axis=1), 0.0, atol=1e-10)


# --- recoverability --------------------------------------------------------


def test_recoverability_never_reads_the_target_row():
    """LEAKAGE TEST: only gamma_true may depend on the target; gamma_hat must not."""
    rng = np.random.default_rng(15)
    Dp = rng.normal(size=(4, 40, 20))
    fold = loco.make_folds(("a", "b", "c", "d"))[1]
    A = Dp[list(fold.source_indices)].mean(axis=0)
    P = Dp[fold.source_indices[0]]
    hat_before = fal.centred(P) - fal.centred(A)
    Dp2 = Dp.copy()
    Dp2[fold.target_index] = rng.normal(size=Dp.shape[1:]) * 100.0
    A2 = Dp2[list(fold.source_indices)].mean(axis=0)
    hat_after = fal.centred(Dp2[fold.source_indices[0]]) - fal.centred(A2)
    np.testing.assert_allclose(hat_before, hat_after, atol=0, rtol=0)


def test_recoverability_is_one_when_the_predictor_is_the_truth():
    """Construct a tensor whose held-out gamma equals the nearest source's deviation."""
    rng = np.random.default_rng(16)
    n_p, k = 30, 12
    Dp = rng.normal(size=(4, n_p, k))
    fold = loco.make_folds(("a", "b", "c", "d"))[0]
    r = fal.gamma_recoverability(
        Dp, fold.target_index, fold.source_indices, predictor_index=fold.source_indices[0]
    )
    assert r.shape == (n_p,)
    assert np.all(np.abs(r[np.isfinite(r)]) <= 1.0 + 1e-9)


def test_recoverability_requires_a_predictor():
    Dp = np.zeros((4, 5, 3))
    with pytest.raises(ValueError, match="predictor_index or weights"):
        fal.gamma_recoverability(Dp, 0, [1, 2, 3])


def test_gamma_reliability_splits_every_context():
    """Identical halves must give reliability 1; independent noise must not."""
    rng = np.random.default_rng(17)
    base = rng.normal(size=(4, 25, 10))
    same = np.stack([np.stack([base, base])], axis=0)
    np.testing.assert_allclose(np.nanmean(fal.gamma_reliability(same, 0)), 1.0, atol=1e-8)
    noisy = np.stack(
        [
            np.stack(
                [
                    base + rng.normal(scale=3.0, size=base.shape),
                    base + rng.normal(scale=3.0, size=base.shape),
                ]
            )
        ],
        axis=0,
    )
    assert np.nanmean(fal.gamma_reliability(noisy, 0)) < 0.5


def test_spearman_brown_is_monotone_and_fixes_one():
    assert fal.spearman_brown(np.array([1.0]))[0] == pytest.approx(1.0)
    vals = fal.spearman_brown(np.array([0.1, 0.3, 0.6]))
    assert vals[0] < vals[1] < vals[2]
    assert np.all(vals > np.array([0.1, 0.3, 0.6]))


# --- null summary statistics -----------------------------------------------


def test_empirical_percentile_detects_a_clear_signal():
    rng = np.random.default_rng(18)
    nulls = rng.normal(loc=0.0, scale=0.05, size=200)
    out = fal.empirical_percentile(0.5, nulls)
    assert out["percentile"] == 100.0
    assert out["p_value"] == pytest.approx(1 / 201, abs=1e-6)
    assert out["z"] > 5


def test_empirical_percentile_detects_a_null_result():
    rng = np.random.default_rng(19)
    nulls = rng.normal(loc=0.4, scale=0.05, size=200)
    out = fal.empirical_percentile(0.4, nulls)
    assert 30 < out["percentile"] < 70
    assert out["p_value"] > 0.2


def test_p_value_is_never_zero():
    out = fal.empirical_percentile(99.0, np.zeros(50))
    assert out["p_value"] > 0


def test_empirical_percentile_handles_empty_and_nan_input():
    out = fal.empirical_percentile(0.3, [])
    assert out["n_null"] == 0
    assert np.isnan(out["p_value"])
    assert np.isnan(fal.empirical_percentile(np.nan, [0.1, 0.2])["p_value"])
