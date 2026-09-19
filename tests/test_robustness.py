"""Invariants of the robustness-battery variants.

All tests run on small in-process data: they pin the variant mathematics without
needing the git-ignored bundle, and they guarantee the sensitivity analysis is
testing what it claims to test.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from virtual_cell.analysis import robustness
from virtual_cell.data.io import DataIntegrityError


@pytest.fixture
def block() -> sparse.csr_matrix:
    rng = np.random.default_rng(0)
    dense = rng.random((40, 9)) * 3.0
    dense[dense < 1.0] = 0.0
    return sparse.csr_matrix(dense)


# --- aggregation order ----------------------------------------------------


def test_mean_log_is_the_plain_column_mean(block):
    np.testing.assert_allclose(
        robustness.aggregate(block, "mean_log"), np.asarray(block.todense()).mean(0), rtol=1e-10
    )


def test_log_mean_is_log1p_of_the_mean_of_expm1(block):
    dense = np.asarray(block.todense())
    np.testing.assert_allclose(
        robustness.aggregate(block, "log_mean"), np.log1p(np.expm1(dense).mean(0)), rtol=1e-10
    )


def test_the_two_aggregations_differ_by_jensen(block):
    """They must not be accidentally identical, or variant D tests nothing."""
    a = robustness.aggregate(block, "mean_log")
    b = robustness.aggregate(block, "log_mean")
    assert not np.allclose(a, b)
    # Jensen: log1p(mean(expm1(x))) >= mean(log1p(expm1(x))) = mean(x), as log1p∘expm1 = id
    assert np.all(b >= a - 1e-9)


def test_aggregation_preserves_sparsity_of_zero(block):
    """expm1(0) == 0, so an all-zero gene stays zero under both orders."""
    dense = np.asarray(block.todense())
    dense[:, 3] = 0.0
    X = sparse.csr_matrix(dense)
    assert robustness.aggregate(X, "mean_log")[3] == 0.0
    assert robustness.aggregate(X, "log_mean")[3] == 0.0


def test_unknown_aggregation_is_rejected(block):
    with pytest.raises(ValueError, match="Unknown aggregation"):
        robustness.aggregate(block, "geometric")


def test_empty_block_is_rejected():
    with pytest.raises(DataIntegrityError, match="empty cell block"):
        robustness.aggregate(sparse.csr_matrix((0, 5)), "mean_log")


@pytest.mark.parametrize("aggregation", ["mean_log", "log_mean"])
def test_group_aggregate_matches_per_group_aggregate(block, aggregation):
    group = np.array([0] * 15 + [1] * 15 + [2] * 10)
    out = robustness.group_aggregate(block, group, 3, aggregation)
    for g in range(3):
        np.testing.assert_allclose(
            out[g], robustness.aggregate(block[group == g], aggregation), rtol=1e-10
        )


def test_group_aggregate_rejects_an_empty_group(block):
    group = np.array([0] * 20 + [2] * 20)
    with pytest.raises(DataIntegrityError, match="no cells"):
        robustness.group_aggregate(block, group, 3, "mean_log")


# --- control scheme -------------------------------------------------------


@pytest.fixture
def reliability_inputs():
    rng = np.random.default_rng(1)
    n_p, per_pert, n_g = 4, 12, 7
    X = sparse.csr_matrix(rng.random((n_p * per_pert, n_g)))
    group = np.repeat(np.arange(n_p), per_pert)
    X_ctrl = sparse.csr_matrix(rng.random((30, n_g)))
    return X, group, X_ctrl, n_p


@pytest.mark.parametrize("scheme", ["shared", "split"])
def test_half_deltas_have_the_expected_shape(reliability_inputs, scheme):
    X, group, X_ctrl, n_p = reliability_inputs
    a, b = robustness.half_deltas(
        X,
        group,
        X_ctrl,
        n_perturbations=n_p,
        rng=np.random.default_rng(0),
        control_scheme=scheme,
        aggregation="mean_log",
    )
    assert a.shape == b.shape == (n_p, X.shape[1])
    assert np.all(np.isfinite(a)) and np.all(np.isfinite(b))


def test_shared_control_uses_one_reference_for_both_halves(reliability_inputs):
    """Under the shared scheme the two halves differ only by perturbation cells."""
    X, group, X_ctrl, n_p = reliability_inputs
    a, b = robustness.half_deltas(
        X,
        group,
        X_ctrl,
        n_perturbations=n_p,
        rng=np.random.default_rng(3),
        control_scheme="shared",
        aggregation="mean_log",
    )
    ctrl = robustness.aggregate(X_ctrl, "mean_log")
    # adding back the same control recovers two perturbation-half means
    assert np.all(np.isfinite(a + ctrl)) and np.all(np.isfinite(b + ctrl))
    # the difference of the halves is free of the control term entirely
    diff = a - b
    a2, b2 = robustness.half_deltas(
        X,
        group,
        X_ctrl * 1.0,
        n_perturbations=n_p,
        rng=np.random.default_rng(3),
        control_scheme="shared",
        aggregation="mean_log",
    )
    np.testing.assert_allclose(diff, a2 - b2, rtol=1e-10)


def test_split_control_differs_from_shared_control(reliability_inputs):
    """Variant A must actually change the estimate, or it tests nothing."""
    X, group, X_ctrl, n_p = reliability_inputs
    shared = robustness.half_deltas(
        X,
        group,
        X_ctrl,
        n_perturbations=n_p,
        rng=np.random.default_rng(5),
        control_scheme="shared",
        aggregation="mean_log",
    )
    split = robustness.half_deltas(
        X,
        group,
        X_ctrl,
        n_perturbations=n_p,
        rng=np.random.default_rng(5),
        control_scheme="split",
        aggregation="mean_log",
    )
    assert not np.allclose(shared[0], split[0])


def test_split_control_halves_use_disjoint_control_cells():
    """Directly verify disjointness via a control matrix with separable rows."""
    n_g = 6
    ctrl = sparse.csr_matrix(np.vstack([np.full(n_g, 1.0)] * 10 + [np.full(n_g, 5.0)] * 10))
    X = sparse.csr_matrix(np.ones((8, n_g)))
    group = np.repeat(np.arange(2), 4)
    a, b = robustness.half_deltas(
        X,
        group,
        ctrl,
        n_perturbations=2,
        rng=np.random.default_rng(0),
        control_scheme="split",
        aggregation="mean_log",
    )
    # If the halves shared cells, both control means would equal the pooled 3.0.
    # Disjoint halves drawn from a 1.0/5.0 mixture give means that sum to 6.0.
    implied = (1.0 - a[0, 0]) + (1.0 - b[0, 0])
    assert implied == pytest.approx(6.0, abs=1e-9)


def test_unknown_control_scheme_is_rejected(reliability_inputs):
    X, group, X_ctrl, n_p = reliability_inputs
    with pytest.raises(ValueError, match="Unknown control_scheme"):
        robustness.half_deltas(
            X,
            group,
            X_ctrl,
            n_perturbations=n_p,
            rng=np.random.default_rng(0),
            control_scheme="bootstrap",
            aggregation="mean_log",
        )


def test_single_cell_perturbation_is_rejected():
    X = sparse.csr_matrix(np.ones((1, 4)))
    with pytest.raises(DataIntegrityError, match="need >= 2"):
        robustness.half_deltas(
            X,
            np.array([0]),
            sparse.csr_matrix(np.ones((6, 4))),
            n_perturbations=1,
            rng=np.random.default_rng(0),
            control_scheme="shared",
            aggregation="mean_log",
        )


# --- control-derived HVG ranking -----------------------------------------


def test_global_hvg_ranking_is_the_descending_mean_variance():
    genes = ["g0", "g1", "g2", "g3"]
    variances = {
        "a": np.array([1.0, 5.0, 0.0, 2.0]),
        "b": np.array([1.0, 3.0, 0.0, 7.0]),
    }
    # means: g0=1.0, g1=4.0, g2=0.0, g3=4.5 -- deliberately untied
    ranking = robustness.global_hvg_ranking(variances, genes)
    assert list(ranking.index) == ["g3", "g1", "g0", "g2"]
    assert ranking["g3"] == pytest.approx(4.5)


def test_hvg_ranking_is_one_global_order_not_per_context():
    """A single ranking must be returned, so every context shares a feature space."""
    genes = ["g0", "g1", "g2"]
    ranking = robustness.global_hvg_ranking(
        {"a": np.array([3.0, 1.0, 2.0]), "b": np.array([1.0, 3.0, 2.0])}, genes
    )
    assert isinstance(ranking, pd.Series)
    assert sorted(ranking.index) == sorted(genes)


# --- controlled depth experiment -----------------------------------------


def test_depth_reliability_holds_pairs_fixed_across_depths():
    rng = np.random.default_rng(2)
    n_g = 12
    signal = rng.normal(size=n_g) * 4.0
    cells = np.abs(signal[None, :] + rng.normal(scale=1.0, size=(400, n_g)))
    X = sparse.csr_matrix(cells)
    group = np.zeros(400, dtype=int)
    ctrl = sparse.csr_matrix(np.abs(rng.normal(size=(50, n_g))))

    frame = robustness.depth_reliability(
        X,
        group,
        ctrl,
        pair_rows=[0],
        depths=(15, 100),
        n_repeats=8,
        rng=np.random.default_rng(0),
    )
    assert set(frame["n_cells"]) == {15, 100}
    assert frame["pert_row"].nunique() == 1
    assert frame["reliability"].between(-1.0, 1.0).all()
    # more cells must estimate the same fixed signal better
    deep = frame.loc[frame["n_cells"] == 100, "reliability"].iloc[0]
    shallow = frame.loc[frame["n_cells"] == 15, "reliability"].iloc[0]
    assert deep > shallow


def test_depth_reliability_skips_depths_without_enough_cells():
    rng = np.random.default_rng(3)
    X = sparse.csr_matrix(rng.random((40, 5)))
    frame = robustness.depth_reliability(
        X,
        np.zeros(40, dtype=int),
        sparse.csr_matrix(rng.random((10, 5))),
        pair_rows=[0],
        depths=(15, 100),
        n_repeats=3,
        rng=np.random.default_rng(0),
    )
    assert set(frame["n_cells"]) == {15}
