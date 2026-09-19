"""Mathematical invariants of the four-component response decomposition.

These tests are data-independent: they pin the algebra and the conventions that
``reports/molina_zhang_reproduction_spec.md`` freezes, so our implementation can
be trusted before the authors' processed artifacts are ever obtained. Where a
ground truth is needed it is *planted* in synthetic tensors rather than copied
from the paper, so nothing here is tuned to a target number.
"""

from __future__ import annotations

import numpy as np
import pytest

from virtual_cell.decomposition import anova

RNG_SEED = 12345
REFERENCE_SEED_FOR_TEST = anova.REFERENCE_SEED


@pytest.fixture
def tensor() -> np.ndarray:
    rng = np.random.default_rng(RNG_SEED)
    return rng.normal(size=(4, 23, 17))


@pytest.fixture
def cl_deltas() -> dict[str, dict[str, np.ndarray]]:
    rng = np.random.default_rng(RNG_SEED)
    perts = [f"P{i:03d}" for i in range(11)]
    return {cl: {p: rng.normal(size=7) for p in perts} for cl in ("k562", "rpe1", "hepg2")}


# --- exact reconstruction -------------------------------------------------


def test_decomposition_reconstructs_the_tensor_exactly(tensor):
    dec = anova.decompose(tensor)
    np.testing.assert_allclose(dec.reconstruct(), tensor, rtol=0, atol=1e-10)


def test_reconstruction_holds_for_degenerate_shapes():
    rng = np.random.default_rng(0)
    for shape in [(1, 1, 1), (2, 1, 5), (1, 6, 5), (7, 2, 1)]:
        D = rng.normal(size=shape)
        np.testing.assert_allclose(anova.decompose(D).reconstruct(), D, atol=1e-10)


# --- zero-sum side conditions ---------------------------------------------


def test_alpha_sums_to_zero_over_cell_lines(tensor):
    dec = anova.decompose(tensor)
    np.testing.assert_allclose(dec.alpha.sum(axis=0), 0.0, atol=1e-10)


def test_beta_sums_to_zero_over_perturbations(tensor):
    dec = anova.decompose(tensor)
    np.testing.assert_allclose(dec.beta.sum(axis=0), 0.0, atol=1e-10)


def test_gamma_sums_to_zero_over_both_axes(tensor):
    dec = anova.decompose(tensor)
    np.testing.assert_allclose(dec.gamma.sum(axis=0), 0.0, atol=1e-10)
    np.testing.assert_allclose(dec.gamma.sum(axis=1), 0.0, atol=1e-10)


# --- balanced-design orthogonality ----------------------------------------


def test_sums_of_squares_partition_the_total_exactly(tensor):
    """The defining property of the reference SS convention."""
    dec = anova.decompose(tensor)
    ss = anova.sums_of_squares(dec)
    total = anova.total_sum_of_squares(tensor)
    assert sum(ss.values()) == pytest.approx(total, rel=1e-10)


def test_components_are_mutually_orthogonal(tensor):
    """Every cross term in the expansion of the total vanishes."""
    dec = anova.decompose(tensor)
    n_c, n_p, _ = tensor.shape
    mu = np.broadcast_to(dec.mu, tensor.shape)
    alpha = np.broadcast_to(dec.alpha[:, None, :], tensor.shape)
    beta = np.broadcast_to(dec.beta[None, :, :], tensor.shape)
    parts = {"mu": mu, "alpha": alpha, "beta": beta, "gamma": dec.gamma}
    names = list(parts)
    scale = n_c * n_p
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            cross = float(np.sum(parts[a] * parts[b])) / scale
            assert cross == pytest.approx(0.0, abs=1e-9), f"{a}.{b} = {cross}"


def test_uncorrected_fractions_sum_to_one(tensor):
    dec = anova.decompose(tensor)
    fractions = anova.uncorrected_fractions(dec, total=anova.total_sum_of_squares(tensor))
    assert sum(fractions.values()) == pytest.approx(1.0, rel=1e-10)


# --- planted-component recovery -------------------------------------------


def test_planted_components_are_recovered_exactly():
    """Build a tensor from known zero-sum components; the estimator must invert it."""
    rng = np.random.default_rng(7)
    n_c, n_p, G = 4, 31, 13
    mu = rng.normal(size=G)
    alpha = rng.normal(size=(n_c, G))
    alpha -= alpha.mean(axis=0)
    beta = rng.normal(size=(n_p, G))
    beta -= beta.mean(axis=0)
    gamma = rng.normal(size=(n_c, n_p, G))
    gamma -= gamma.mean(axis=0, keepdims=True)
    gamma -= gamma.mean(axis=1, keepdims=True)

    D = mu + alpha[:, None, :] + beta[None, :, :] + gamma
    dec = anova.decompose(D)

    np.testing.assert_allclose(dec.mu, mu, atol=1e-10)
    np.testing.assert_allclose(dec.alpha, alpha, atol=1e-10)
    np.testing.assert_allclose(dec.beta, beta, atol=1e-10)
    np.testing.assert_allclose(dec.gamma, gamma, atol=1e-10)


def test_planted_variance_shares_are_recovered():
    """A tensor built with a known SS budget must report that budget back."""
    rng = np.random.default_rng(11)
    n_c, n_p, G = 4, 200, 50

    def centred(shape, axes):
        x = rng.normal(size=shape)
        for ax in axes:
            x = x - x.mean(axis=ax, keepdims=True)
        return x

    mu = rng.normal(size=G) * 3.0
    alpha = centred((n_c, G), [0]) * 2.0
    beta = centred((n_p, G), [0]) * 1.5
    gamma = centred((n_c, n_p, G), [0, 1]) * 1.0

    D = mu + alpha[:, None, :] + beta[None, :, :] + gamma
    dec = anova.decompose(D)
    ss = anova.sums_of_squares(dec)

    expected = {
        "mu": float(np.sum(mu**2)),
        "alpha": float(np.mean(np.sum(alpha**2, axis=1))),
        "beta": float(np.mean(np.sum(beta**2, axis=1))),
        "gamma": float(np.mean(np.sum(gamma**2, axis=2))),
    }
    for key, value in expected.items():
        assert ss[key] == pytest.approx(value, rel=1e-8), key


# --- permutation invariance -----------------------------------------------


def test_variance_shares_are_invariant_to_perturbation_order(tensor):
    rng = np.random.default_rng(3)
    order = rng.permutation(tensor.shape[1])
    base = anova.sums_of_squares(anova.decompose(tensor))
    permuted = anova.sums_of_squares(anova.decompose(tensor[:, order, :]))
    for key in base:
        assert permuted[key] == pytest.approx(base[key], rel=1e-10), key


def test_variance_shares_are_invariant_to_cell_line_order(tensor):
    rng = np.random.default_rng(4)
    order = rng.permutation(tensor.shape[0])
    base = anova.sums_of_squares(anova.decompose(tensor))
    permuted = anova.sums_of_squares(anova.decompose(tensor[order, :, :]))
    for key in base:
        assert permuted[key] == pytest.approx(base[key], rel=1e-10), key


def test_beta_follows_a_perturbation_permutation(tensor):
    rng = np.random.default_rng(5)
    order = rng.permutation(tensor.shape[1])
    base = anova.decompose(tensor)
    permuted = anova.decompose(tensor[:, order, :])
    np.testing.assert_allclose(permuted.beta, base.beta[order], atol=1e-10)


# --- template projection --------------------------------------------------


def test_projection_removes_the_template_direction():
    rng = np.random.default_rng(9)
    D = rng.normal(size=(3, 40, 12)) + 5.0
    resid = anova.project_out_template(D)
    for c in range(D.shape[0]):
        unit = D[c].mean(0) / np.linalg.norm(D[c].mean(0))
        np.testing.assert_allclose(resid[c] @ unit, 0.0, atol=1e-9)


def test_projection_never_increases_the_total(tensor):
    assert (
        anova.total_sum_of_squares(anova.project_out_template(tensor))
        <= anova.total_sum_of_squares(tensor) + 1e-12
    )


# --- split-half noise correction ------------------------------------------


def _synthetic_cell_data(n_cells, noise_scale, seed=0, n_c=4, n_p=25, G=20):
    """Cells drawn around a planted per-(cl, pert) mean, with known noise."""
    rng = np.random.default_rng(seed)
    cell_lines = [f"cl{i}" for i in range(n_c)]
    perts = [f"P{j:03d}" for j in range(n_p)]
    truth = rng.normal(size=(n_c, n_p, G))
    data = {}
    for ci, cl in enumerate(cell_lines):
        pert_cells = {
            p: truth[ci, pj] + rng.normal(scale=noise_scale, size=(n_cells, G))
            for pj, p in enumerate(perts)
        }
        data[cl] = {"ctrl_mean": np.zeros(G), "pert_cells": pert_cells}
    return data, cell_lines, perts, truth


def test_cross_half_signal_is_unbiased_for_a_noise_free_tensor(tensor):
    """With identical halves the cross-product reduces to the plain SS."""
    signal = anova.cross_half_signal(tensor, tensor)
    ss = anova.sums_of_squares(anova.decompose(tensor))
    for key in ss:
        assert signal[key] == pytest.approx(ss[key], rel=1e-10), key


def test_split_half_correction_recovers_a_known_noise_share():
    """Noise share must track the planted noise level, and the shares must sum to 1."""
    data, cls, perts, _ = _synthetic_cell_data(n_cells=200, noise_scale=1.0, seed=2)
    D = np.array(
        [[data[cl]["pert_cells"][p].mean(0) - data[cl]["ctrl_mean"] for p in perts] for cl in cls]
    )
    total = anova.total_sum_of_squares(D)
    signal = anova.split_half_signal(data, cell_lines=cls, perturbations=perts, n_splits=20, seed=1)
    fractions = anova.noise_corrected_fractions(signal, total)

    assert sum(fractions.values()) == pytest.approx(1.0, rel=1e-10)
    assert all(-0.05 < v < 1.05 for v in fractions.values()), fractions
    # 200 cells at unit noise: the pseudobulk mean is precise, so noise is small.
    assert fractions["noise"] < 0.15


def test_noise_share_grows_when_cells_per_perturbation_shrink():
    """Monotonicity is the qualitative behaviour the correction exists to capture."""
    shares = []
    for n_cells in (400, 50, 12):
        data, cls, perts, _ = _synthetic_cell_data(n_cells=n_cells, noise_scale=1.0, seed=3)
        D = np.array(
            [
                [data[cl]["pert_cells"][p].mean(0) - data[cl]["ctrl_mean"] for p in perts]
                for cl in cls
            ]
        )
        signal = anova.split_half_signal(
            data, cell_lines=cls, perturbations=perts, n_splits=10, seed=1
        )
        shares.append(
            anova.noise_corrected_fractions(signal, anova.total_sum_of_squares(D))["noise"]
        )
    assert shares[0] < shares[1] < shares[2], shares


def test_split_half_is_reproducible_under_a_fixed_seed():
    data, cls, perts, _ = _synthetic_cell_data(n_cells=40, noise_scale=1.0, seed=4)
    kwargs = dict(cell_lines=cls, perturbations=perts, n_splits=5, seed=REFERENCE_SEED_FOR_TEST)
    first = anova.split_half_signal(data, **kwargs)
    second = anova.split_half_signal(data, **kwargs)
    assert first == second


def test_split_half_rejects_a_single_cell_perturbation():
    data, cls, perts, _ = _synthetic_cell_data(n_cells=1, noise_scale=1.0, seed=5)
    with pytest.raises(anova.DecompositionError, match="at least 2"):
        anova.split_half_signal(data, cell_lines=cls, perturbations=perts, n_splits=1)


# --- frozen sets and malformed input --------------------------------------


def test_shared_perturbations_are_the_sorted_intersection(cl_deltas):
    del cl_deltas["k562"]["P000"]
    shared = anova.shared_perturbations(cl_deltas)
    assert shared == sorted(shared)
    assert "P000" not in shared
    assert len(shared) == 10


def test_frozen_perturbation_set_is_honoured_in_order(cl_deltas):
    frozen = ["P005", "P001", "P009"]
    D, cls, perts = anova.build_response_tensor(cl_deltas, perturbations=frozen)
    assert perts == frozen
    assert D.shape[1] == 3
    np.testing.assert_allclose(D[cls.index("rpe1"), 0], cl_deltas["rpe1"]["P005"])


def test_frozen_cell_line_order_is_honoured(cl_deltas):
    order = ["hepg2", "k562", "rpe1"]
    D, cls, _ = anova.build_response_tensor(cl_deltas, cell_lines=order)
    assert cls == order
    np.testing.assert_allclose(D[0, 0], cl_deltas["hepg2"][sorted(cl_deltas["hepg2"])[0]])


def test_frozen_response_gene_set_is_enforced(cl_deltas):
    cl_deltas["k562"]["P000"] = np.zeros(6)  # every other vector has length 7
    with pytest.raises(anova.DecompositionError, match="Inconsistent response dimensionality"):
        anova.build_response_tensor(cl_deltas)


def test_unbalanced_design_is_rejected_when_the_set_is_frozen(cl_deltas):
    del cl_deltas["hepg2"]["P003"]
    with pytest.raises(anova.DecompositionError, match="Unbalanced design"):
        anova.build_response_tensor(cl_deltas, perturbations=["P002", "P003"])


def test_missing_cell_line_is_rejected(cl_deltas):
    with pytest.raises(anova.DecompositionError, match="absent from the input"):
        anova.build_response_tensor(cl_deltas, cell_lines=["k562", "jurkat"])


def test_no_shared_perturbations_is_rejected():
    with pytest.raises(anova.DecompositionError, match="share no perturbations"):
        anova.shared_perturbations({"a": {"P1": np.zeros(3)}, "b": {"P2": np.zeros(3)}})


def test_empty_input_is_rejected():
    with pytest.raises(anova.DecompositionError, match="No cell lines"):
        anova.shared_perturbations({})


@pytest.mark.parametrize(
    "bad",
    [
        np.zeros((3, 4)),
        np.zeros((3, 4, 5, 6)),
        np.zeros((0, 4, 5)),
        np.zeros((3, 0, 5)),
    ],
)
def test_malformed_tensor_shapes_are_rejected(bad):
    with pytest.raises(anova.DecompositionError):
        anova.decompose(bad)


@pytest.mark.parametrize("bad_value", [np.nan, np.inf, -np.inf])
def test_non_finite_values_are_rejected(tensor, bad_value):
    tensor = tensor.copy()
    tensor[0, 0, 0] = bad_value
    with pytest.raises(anova.DecompositionError, match="NaN or infinite"):
        anova.decompose(tensor)


def test_mismatched_halves_are_rejected(tensor):
    with pytest.raises(anova.DecompositionError, match="Half shapes differ"):
        anova.cross_half_signal(tensor, tensor[:, :-1, :])


def test_zero_total_is_rejected():
    dec = anova.decompose(np.zeros((2, 3, 4)))
    with pytest.raises(anova.DecompositionError, match="not positive"):
        anova.uncorrected_fractions(dec)


def test_beta_fraction_is_one_when_there_is_no_interaction():
    rng = np.random.default_rng(21)
    beta = rng.normal(size=(9, 6))
    beta -= beta.mean(axis=0)
    D = np.broadcast_to(beta[None, :, :], (4, 9, 6)).copy()
    dec = anova.decompose(D, perturbations=[f"P{i}" for i in range(9)])
    for value in anova.beta_fraction_per_perturbation(dec).values():
        assert value == pytest.approx(1.0, abs=1e-8)


# --- equivalence with the authors' released implementation ----------------


def _reference_fig1f_algebra(D):
    """Verbatim transcription of the reference implementation's algebra.

    Source: xinyizhanglab/perturbation-decomposition @ a152147,
    ``figures/fig1_f.py`` lines 78-90 (decomposition and sums of squares).
    Kept literal, including the Python-level loops, so that any drift between
    our vectorised code and the reference is caught here.
    """
    n_cl, n_p, _ = D.shape
    mu = D.mean(axis=(0, 1))
    alpha = D.mean(axis=1) - mu
    beta = D.mean(axis=0) - mu
    gamma = D - mu[None, None, :] - alpha[:, None, :] - beta[None, :, :]
    ss_total = np.mean([np.sum(D[c, p] ** 2) for c in range(n_cl) for p in range(n_p)])
    ss = {
        "mu": np.sum(mu**2),
        "alpha": np.mean([np.sum(alpha[c] ** 2) for c in range(n_cl)]),
        "beta": np.mean([np.sum(beta[p] ** 2) for p in range(n_p)]),
        "gamma": np.mean([np.sum(gamma[c, p] ** 2) for c in range(n_cl) for p in range(n_p)]),
    }
    return mu, alpha, beta, gamma, float(ss_total), ss


def test_matches_the_reference_implementation_algebra(tensor):
    mu, alpha, beta, gamma, ss_total, ss = _reference_fig1f_algebra(tensor)
    dec = anova.decompose(tensor)

    np.testing.assert_allclose(dec.mu, mu, atol=1e-12)
    np.testing.assert_allclose(dec.alpha, alpha, atol=1e-12)
    np.testing.assert_allclose(dec.beta, beta, atol=1e-12)
    np.testing.assert_allclose(dec.gamma, gamma, atol=1e-12)

    assert anova.total_sum_of_squares(tensor) == pytest.approx(ss_total, rel=1e-12)
    ours = anova.sums_of_squares(dec)
    for key, value in ss.items():
        assert ours[key] == pytest.approx(float(value), rel=1e-12), key


def test_matches_the_reference_projective_template_removal(tensor):
    """Source: ``figures/fig1_f.py`` lines 94-104."""
    D_out = np.zeros_like(tensor)
    for ci in range(tensor.shape[0]):
        T_c = tensor[ci].mean(0)
        Tn_c = T_c / (np.linalg.norm(T_c) + 1e-30)
        for pi in range(tensor.shape[1]):
            D_out[ci, pi] = tensor[ci, pi] - np.dot(tensor[ci, pi], Tn_c) * Tn_c
    np.testing.assert_allclose(anova.project_out_template(tensor), D_out, atol=1e-12)


def test_matches_the_reference_cross_half_signal(tensor):
    """Source: ``figures/fig1_f.py`` lines 131-144 (``anova_cross_sigs``)."""
    rng = np.random.default_rng(99)
    other = tensor + rng.normal(scale=0.3, size=tensor.shape)
    n_cl, n_p, _ = tensor.shape

    mu1 = tensor.mean(axis=(0, 1))
    mu2 = other.mean(axis=(0, 1))
    a1 = tensor.mean(axis=1) - mu1
    a2 = other.mean(axis=1) - mu2
    b1 = tensor.mean(axis=0) - mu1
    b2 = other.mean(axis=0) - mu2
    g1 = tensor - mu1 - a1[:, None, :] - b1[None, :, :]
    g2 = other - mu2 - a2[:, None, :] - b2[None, :, :]
    expected = {
        "mu": np.dot(mu1, mu2),
        "alpha": np.mean([np.dot(a1[c], a2[c]) for c in range(n_cl)]),
        "beta": np.mean([np.dot(b1[p], b2[p]) for p in range(n_p)]),
        "gamma": np.mean([np.dot(g1[c, p], g2[c, p]) for c in range(n_cl) for p in range(n_p)]),
    }
    ours = anova.cross_half_signal(tensor, other)
    for key, value in expected.items():
        assert ours[key] == pytest.approx(float(value), rel=1e-10), key
