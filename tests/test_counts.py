"""Count recovery from ``log1p(CP10K)``: exactness, the failure modes, the refusals.

The single-cell half of the Arc benchmark rests entirely on this inversion
being the original count matrix rather than something that merely looks like
one, so the tests plant known counts and demand them back exactly.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy import sparse

from virtual_cell.data.counts import (
    TARGET_SUM,
    recover_counts,
    recovery_diagnostics,
    row_minimum,
)


def normalise(counts: np.ndarray, *, dtype=np.float32) -> sparse.csr_matrix:
    """Apply exactly the transform scPertEval documents, at its storage dtype."""
    counts = np.asarray(counts, dtype=np.float64)
    lib = counts.sum(axis=1, keepdims=True)
    cpm = np.divide(counts * TARGET_SUM, lib, out=np.zeros_like(counts), where=lib > 0)
    return sparse.csr_matrix(np.log1p(cpm).astype(dtype).astype(np.float64))


def planted(rng: np.random.Generator, n_cells: int, n_genes: int, lam: float = 1.5) -> np.ndarray:
    counts = rng.poisson(lam, size=(n_cells, n_genes)).astype(np.int64)
    counts[counts.sum(axis=1) == 0, 0] = 1  # no empty cells
    return counts


# --------------------------------------------------------------------------
# exactness
# --------------------------------------------------------------------------


def test_planted_counts_are_recovered_exactly():
    rng = np.random.default_rng(0)
    counts = planted(rng, 60, 200)
    recovered = recover_counts(normalise(counts))
    assert np.array_equal(recovered.counts.toarray(), counts)


def test_library_sizes_are_recovered_exactly():
    rng = np.random.default_rng(1)
    counts = planted(rng, 40, 150)
    recovered = recover_counts(normalise(counts))
    assert np.array_equal(recovered.library_sizes, counts.sum(axis=1))


def test_recovery_is_exact_at_float64_storage():
    rng = np.random.default_rng(2)
    counts = planted(rng, 30, 120)
    recovered = recover_counts(normalise(counts, dtype=np.float64))
    assert np.array_equal(recovered.counts.toarray(), counts)
    assert recovered.max_integrality_error < 1e-8


def test_deep_cells_are_recovered():
    """Depth is where float32 storage bites; 20,000 UMIs must still invert."""
    rng = np.random.default_rng(3)
    counts = rng.poisson(40, size=(25, 500)).astype(np.int64)
    recovered = recover_counts(normalise(counts))
    assert np.array_equal(recovered.counts.toarray(), counts)


# --------------------------------------------------------------------------
# the quantum: cells whose rarest gene carries more than one count
# --------------------------------------------------------------------------


def test_quantum_is_one_when_a_singleton_exists():
    rng = np.random.default_rng(4)
    counts = planted(rng, 20, 100)
    counts[:, 0] = 1
    assert np.all(recover_counts(normalise(counts)).row_quantum == 1)


def test_a_row_without_a_singleton_returns_the_minimal_solution():
    """The documented ambiguity, pinned so it cannot be forgotten.

    A cell whose counts share a common factor is indistinguishable from the
    cell with those counts divided through: both normalise to the same row.
    The recovery returns the minimal one, and this test states which.
    """
    counts = np.array([[3, 6, 9, 0, 12], [6, 3, 0, 3, 3]], dtype=np.int64)
    recovered = recover_counts(normalise(counts, dtype=np.float64))
    assert np.array_equal(recovered.counts.toarray(), counts // 3)
    assert np.array_equal(recovered.library_sizes, counts.sum(axis=1) // 3)


def test_scaled_and_unscaled_rows_normalise_identically():
    """Why the ambiguity is not a bug that could be fixed with more care."""
    small = np.array([[1, 2, 3, 0, 4]], dtype=np.int64)
    large = small * 7
    assert np.allclose(
        normalise(small, dtype=np.float64).toarray(),
        normalise(large, dtype=np.float64).toarray(),
    )


def test_coprime_rows_are_recovered_exactly_even_without_a_one():
    """Coprimality, not the presence of a 1, is what the recovery needs."""
    counts = np.array([[2, 3, 5, 0, 7], [4, 6, 0, 9, 2]], dtype=np.int64)
    recovered = recover_counts(normalise(counts, dtype=np.float64))
    assert np.array_equal(recovered.counts.toarray(), counts)


# --------------------------------------------------------------------------
# refusals
# --------------------------------------------------------------------------


def test_matrix_normalised_before_gene_filtering_is_refused():
    """Row sums below the target mean the inversion is underdetermined."""
    rng = np.random.default_rng(5)
    counts = planted(rng, 20, 200)
    full = normalise(counts).toarray()
    with pytest.raises(ValueError, match="do not sum"):
        recover_counts(sparse.csr_matrix(full[:, :100]))


def test_non_count_matrix_is_refused():
    rng = np.random.default_rng(6)
    values = rng.gamma(2.0, 1.0, size=(15, 80))
    lib = values.sum(axis=1, keepdims=True)
    x = np.log1p(values * TARGET_SUM / lib)
    with pytest.raises(ValueError, match="not integral"):
        recover_counts(sparse.csr_matrix(x))


def test_negative_stored_value_is_refused():
    x = sparse.csr_matrix(np.array([[1.0, -0.5, 2.0]]))
    with pytest.raises(ValueError, match="strictly positive"):
        recover_counts(x)


# --------------------------------------------------------------------------
# edges and diagnostics
# --------------------------------------------------------------------------


def test_empty_rows_survive_and_are_reported():
    rng = np.random.default_rng(7)
    counts = planted(rng, 10, 60)
    counts[3] = 0
    recovered = recover_counts(normalise(counts))
    assert recovered.n_empty_rows == 1
    assert recovered.library_sizes[3] == 0
    assert np.array_equal(recovered.counts.toarray(), counts)


def test_row_minimum_handles_empty_rows():
    data = np.array([3.0, 1.0, 2.0, 5.0])
    indptr = np.array([0, 3, 3, 4])
    out = row_minimum(data, indptr)
    assert out[0] == 1.0
    assert np.isinf(out[1])
    assert out[2] == 5.0


def test_diagnostics_round_trip_is_tight():
    rng = np.random.default_rng(8)
    counts = planted(rng, 40, 300)
    x = normalise(counts)
    diag = recovery_diagnostics(x, recover_counts(x))
    assert diag["max_round_trip_error"] < 1e-5
    assert diag["median_library_size"] == float(np.median(counts.sum(axis=1)))


def test_sparsity_is_preserved():
    """The inversion must not invent a detected gene or lose one."""
    rng = np.random.default_rng(9)
    counts = planted(rng, 50, 400, lam=0.3)
    recovered = recover_counts(normalise(counts))
    assert np.array_equal((recovered.counts.toarray() > 0), (counts > 0))
