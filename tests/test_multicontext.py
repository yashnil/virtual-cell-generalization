"""Checks on the delta-matrix scoring path used for the Feng multi-context study.

This path exists only because Feng ships log fold changes rather than cells, so
the frozen cell-level entry point does not apply. Its correctness claim is that
it is the *same* scoring, and the real proof of that is in the analysis script,
which reproduces the frozen arch1 numbers before touching Feng. These tests pin
the pieces that script cannot: the algebra of the source block, the
positive-control transfer, and the cross-line decomposition.
"""

from __future__ import annotations

import numpy as np
import pytest

from virtual_cell.modelling import multicontext as mc
from virtual_cell.modelling import unseen_perturbation as up
from virtual_cell.modelling.external_benchmark import KNN_K, PCA_COMPONENTS, RIDGE_ALPHA

N_CTX, N_PERT, N_GENE = 4, 30, 18
GENES = [f"G{i}" for i in range(N_GENE)]
PERTS = [f"P{i}" for i in range(N_PERT)]


@pytest.fixture
def rng() -> np.random.Generator:
    return np.random.default_rng(7)


@pytest.fixture
def delta(rng: np.random.Generator) -> np.ndarray:
    return rng.normal(size=(N_CTX, N_PERT, N_GENE))


def test_source_block_beta_is_centred_over_perturbations(delta: np.ndarray) -> None:
    block = mc.build_source_block(delta, GENES, PERTS, GENES)
    assert block.beta.shape == (N_PERT, N_GENE)
    assert np.abs(block.beta.mean(axis=0)).max() < 1e-12


def test_source_block_restricts_to_the_requested_gene_axis(delta: np.ndarray) -> None:
    subset = GENES[2:10]
    block = mc.build_source_block(delta, GENES, PERTS, subset)
    assert block.delta.shape == (N_CTX, N_PERT, len(subset))
    assert np.allclose(block.delta, delta[:, :, 2:10])


def test_source_block_rejects_a_gene_off_the_axis(delta: np.ndarray) -> None:
    with pytest.raises(ValueError, match="present on the source axis"):
        mc.build_source_block(delta, GENES, PERTS, [*GENES[:3], "NOT_A_GENE"])


def test_measured_transfer_is_near_perfect_when_the_target_matches_the_sources(
    delta: np.ndarray,
) -> None:
    """The positive control must actually fire when transfer is easy.

    Not exactly 1: the truth is centred on the 12 evaluated perturbations while
    the transfer is centred on all 30 source perturbations, so the two sides
    carry different offsets. That asymmetry is in the frozen protocol and is
    pinned here rather than papered over.
    """
    block = mc.build_source_block(delta, GENES, PERTS, GENES)
    seen = PERTS[:12]
    observed = delta[:, :12].mean(axis=0)
    got = mc.score_measured_transfer(observed, seen, block)
    assert got["n_seen"] == 12
    assert got["raw_pearson"] > 0.95


def test_measured_transfer_is_near_zero_against_an_unrelated_target(
    delta: np.ndarray, rng: np.random.Generator
) -> None:
    block = mc.build_source_block(delta, GENES, PERTS, GENES)
    seen = PERTS[:12]
    got = mc.score_measured_transfer(rng.normal(size=(12, N_GENE)), seen, block)
    assert abs(got["raw_pearson"]) < 0.5


def test_measured_transfer_rejects_a_perturbation_the_sources_lack(delta: np.ndarray) -> None:
    block = mc.build_source_block(delta, GENES, PERTS, GENES)
    with pytest.raises(ValueError, match="must be a source perturbation"):
        mc.score_measured_transfer(np.zeros((1, N_GENE)), ["NOT_A_PERT"], block)


def test_cross_line_agreement_is_one_when_every_line_agrees(rng: np.random.Generator) -> None:
    shared = rng.normal(size=(5, N_GENE))
    stack = np.stack([shared, shared, shared])
    got = mc.cross_line_agreement(stack)
    assert np.allclose(got["cross_line_agreement"], 1.0)
    assert np.allclose(got["line_specific_fraction"], 0.0, atol=1e-12)


def test_cross_line_agreement_falls_and_line_specific_energy_rises_with_disagreement(
    rng: np.random.Generator,
) -> None:
    shared = rng.normal(size=(5, N_GENE))
    noisy = np.stack([shared + rng.normal(0, 3.0, size=shared.shape) for _ in range(4)])
    got = mc.cross_line_agreement(noisy)
    assert np.nanmedian(got["cross_line_agreement"]) < 0.5
    assert np.nanmedian(got["line_specific_fraction"]) > 0.4


def test_cross_line_agreement_skips_a_perturbation_present_in_one_line(
    rng: np.random.Generator,
) -> None:
    stack = rng.normal(size=(3, 2, N_GENE))
    stack[1:, 0, :] = np.nan  # only line 0 carries perturbation 0
    got = mc.cross_line_agreement(stack)
    assert np.isnan(got["cross_line_agreement"][0])
    assert np.isfinite(got["cross_line_agreement"][1])


def test_score_unseen_reports_every_frozen_estimator(
    delta: np.ndarray, rng: np.random.Generator
) -> None:
    block = mc.build_source_block(delta, GENES, PERTS, GENES)
    test_names = [f"U{i}" for i in range(6)]
    vocabulary = sorted([*PERTS, *test_names])
    matrix = rng.normal(size=(len(vocabulary), 12))
    observed = rng.normal(size=(len(test_names), N_GENE))
    table, detail = mc.score_unseen(observed, test_names, block, matrix, vocabulary)
    assert list(table["estimator"]) == list(up.ESTIMATORS)
    assert (table["n_test"] == len(test_names)).all()
    assert list(detail["perturbation"]) == test_names
    assert {"neighbour_agreement", "support_distance", "knn_pearson"} <= set(detail.columns)


def test_score_unseen_recovers_a_perfectly_informative_representation(
    delta: np.ndarray,
) -> None:
    """If a test gene's features equal a training gene's, k-NN at k=1 copies it."""
    block = mc.build_source_block(delta, GENES, PERTS, GENES)
    vocabulary = [*PERTS, "COPY"]
    matrix = np.eye(len(vocabulary))
    matrix[-1] = matrix[3]  # "COPY" is feature-identical to PERTS[3]
    observed = (block.beta[3] + block.beta[4])[None, :] / 1.0
    _, detail = mc.score_unseen(observed, ["COPY"], block, matrix, vocabulary)
    assert np.isfinite(detail["neighbour_agreement"]).all()


def test_the_frozen_constants_are_imported_not_redefined() -> None:
    """A local copy of these would let the protocol drift silently."""
    import inspect

    source = inspect.getsource(mc)
    assert "KNN_K = " not in source
    assert "RIDGE_ALPHA = " not in source
    assert "PCA_COMPONENTS = " not in source
    assert (KNN_K, RIDGE_ALPHA, PCA_COMPONENTS) == (25, 10.0, 64)
