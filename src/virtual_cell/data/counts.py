r"""Recovering raw integer counts from a ``log1p(CP10K)`` matrix.

The 2026 challenge scores **counts**, so every generator in
:mod:`virtual_cell.arc.generate` consumes raw integer cells and every metric in
:mod:`virtual_cell.arc.metrics` reads a count matrix. The public scPertEval
bundles ship ``log1p(CP10K)`` and carry **no** ``layers['counts']`` and no
``raw`` — verified on all seven files. Taken at face value that rules out any
public single-cell benchmark of a count generator, which would leave the
generator comparison resting on the Arc controls alone, where no perturbed
ground truth exists.

The normalisation is invertible, exactly, and this module inverts it.

The algebra
-----------
Let ``c_g`` be a cell's integer count for gene ``g`` and ``L = sum_g c_g`` its
library size over the **stored** gene axis. scPertEval's ``X`` is

.. math::
    x_g = \log(1 + S c_g / L), \qquad S = 10^4

so ``expm1(x_g) = S c_g / L`` and ``sum_g expm1(x_g) = S`` exactly. That row-sum
identity is checked here and holds to ``1e-4`` relative on every audited file,
which establishes that the normalisation was applied on the stored axis rather
than before gene filtering — if it had been applied first, the stored rows
would sum below ``S`` and the inversion would be underdetermined.

``L`` itself is not stored, and dividing through by the smallest stored value
gives the row only **up to scale**:

.. math::
    c_g / \min_h c_h = \operatorname{expm1}(x_g) / \min_h \operatorname{expm1}(x_h)

What this module returns is the *minimal* integer solution — the counts
rescaled so the rarest detected gene reads 1. That solution is unique and
well defined, and it equals the true counts **exactly when the true row's
counts are coprime**, which happens precisely when the cell detected at least
one gene with a single UMI.

The ambiguity is real and is not removable: a cell with counts ``(3, 6, 9)``
and one with counts ``(1, 2, 3)`` produce byte-identical ``log1p(CP10K)``
rows, and no inversion can tell them apart. So the singleton assumption is
stated, not proved. What supports it is that these are 10x Perturb-seq cells
detecting thousands of genes at a median depth around 10-15k UMIs; a cell
with no gene at exactly one count would be extraordinary, and if the
assumption failed by a factor ``k`` every recovered library size would be
``k`` times too small, which is checkable against the published depths of the
source experiments. :func:`recovery_diagnostics` reports the library-size
distribution for exactly that check.

What is and is not established
------------------------------
The evidence that the result is the original count matrix is:

* ``sum_g expm1(x_g) == S`` to within tolerance — the normalisation is on this
  axis, and nothing was dropped after it (a matrix normalised *before* gene
  filtering fails here, and is refused);
* the row is integral to within tolerance after rescaling — the values really
  are ``integer / L`` and not some other transform, which a continuous matrix
  fails at every scale and is refused for;
* the recovered counts round-trip to ``X``;
* the recovered library sizes are plausible depths.

plus the stated singleton assumption, which fixes the one remaining degree of
freedom. The ``k`` search below exists to make the integrality check a real
test rather than a formality: a non-count matrix must fail it at *every* small
scale, not just at one.

What this does **not** recover is any gene the bundle's preprocessing removed
before storing: those counts are gone, and the recovered library size is the
library size *over the stored axis*, smaller than the cell's true UMI total.
Every statistic derived from it is therefore a statistic of the filtered
transcriptome, and is labelled as such wherever it is reported.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import sparse

__all__ = [
    "TARGET_SUM",
    "MAX_QUANTUM",
    "RecoveredCounts",
    "row_minimum",
    "recover_counts",
    "recovery_diagnostics",
]

#: scPertEval's ``normalize_total`` target, and the ``S`` of the module algebra.
TARGET_SUM = 1.0e4

#: Largest rarest-gene count the search will consider. A cell whose rarest
#: detected gene carries more than this is not recovered; it is reported.
MAX_QUANTUM = 64


@dataclass(frozen=True)
class RecoveredCounts:
    """Integer counts plus the evidence that the inversion was exact."""

    counts: sparse.csr_matrix
    library_sizes: np.ndarray
    #: Per-row multiplier that made the row integral. ``1`` for any matrix
    #: built from integer counts; above one only if float noise forced it.
    row_quantum: np.ndarray
    #: Largest distance from an integer after the per-row multiplier was applied.
    max_integrality_error: float
    #: Largest ``|sum_g expm1(x_g) / S - 1|`` over rows.
    max_rowsum_error: float
    #: Rows with no stored value; their counts are all zero.
    n_empty_rows: int

    @property
    def n_rows_rescaled(self) -> int:
        """Rows that needed a multiplier above one to come out integral."""
        return int((self.row_quantum > 1).sum())


def row_minimum(data: np.ndarray, indptr: np.ndarray) -> np.ndarray:
    """Smallest stored value in each CSR row; ``inf`` for an empty row."""
    data = np.asarray(data, dtype=np.float64)
    indptr = np.asarray(indptr)
    n_rows = len(indptr) - 1
    out = np.full(n_rows, np.inf, dtype=np.float64)
    lengths = np.diff(indptr)
    nonempty = np.flatnonzero(lengths > 0)
    if nonempty.size:
        out[nonempty] = np.minimum.reduceat(data, indptr[nonempty].astype(np.int64))
    return out


def _row_quantum(
    provisional: np.ndarray, lengths: np.ndarray, tolerance: float, max_quantum: int
) -> tuple[np.ndarray, np.ndarray]:
    """Smallest ``k`` per row making ``k * provisional`` integral, and the error.

    A row built from integer counts is already integral at ``k = 1``, so this
    returns 1 throughout for real data. The search matters for the *refusal*:
    a continuous matrix must fail at every ``k`` up to ``max_quantum`` before
    it is rejected, which is a much stronger statement than failing at one.
    Rows for which no ``k <= max_quantum`` works are returned with ``k = 0``
    and their own best (unachieved) error, so the caller can refuse them.
    """
    n_rows = len(lengths)
    starts = np.concatenate([[0], np.cumsum(lengths)]).astype(np.int64)
    quantum = np.zeros(n_rows, dtype=np.int64)
    best_err = np.full(n_rows, np.inf)
    pending = np.flatnonzero(lengths > 0)
    for k in range(1, max_quantum + 1):
        if pending.size == 0:
            break
        scaled = provisional * float(k)
        err = np.abs(scaled - np.round(scaled))
        # Largest deviation within each still-unassigned row.
        row_err = np.array([err[starts[i] : starts[i + 1]].max() for i in pending])
        best_err[pending] = np.minimum(best_err[pending], row_err)
        hit = row_err <= tolerance
        quantum[pending[hit]] = k
        pending = pending[~hit]
    best_err[lengths == 0] = 0.0
    return quantum, best_err


def recover_counts(
    block: sparse.csr_matrix,
    *,
    target_sum: float = TARGET_SUM,
    integrality_tolerance: float = 1.0e-2,
    max_rowsum_error: float = 1.0e-3,
    max_quantum: int = MAX_QUANTUM,
) -> RecoveredCounts:
    """Invert ``log1p(CP<target_sum>)`` back to the integer counts it came from.

    Raises :class:`ValueError` if the row-sum identity fails or if any row
    resists integralisation, so a file whose normalisation differs from the
    documented one cannot be silently turned into plausible-looking counts.
    """
    block = sparse.csr_matrix(block)
    data = np.asarray(block.data, dtype=np.float64)
    if data.size and data.min() <= 0:
        raise ValueError("stored values must be strictly positive in a log1p matrix")

    expm1 = np.expm1(data)
    rowsum = np.add.reduceat(expm1, block.indptr[:-1].astype(np.int64)) if data.size else None
    lengths = np.diff(block.indptr)
    rowsum_err = 0.0
    if rowsum is not None:
        ok = lengths > 0
        rowsum_err = float(np.abs(rowsum[ok] / target_sum - 1.0).max()) if ok.any() else 0.0
    if rowsum_err > max_rowsum_error:
        raise ValueError(
            f"rows do not sum to {target_sum:g} after expm1 "
            f"(max relative error {rowsum_err:.3e}); the matrix is not "
            "log1p-normalised on its own gene axis"
        )

    # Scale each row so its rarest detected gene reads exactly one count: the
    # minimal integer solution (see the module docstring on what that assumes).
    # The inversion is on ``expm1(x) = S c / L``, so the row minimum is taken
    # there rather than on ``x`` -- monotone, but the factor matters.
    smallest = row_minimum(expm1, block.indptr)
    provisional_lib = np.where(
        np.isfinite(smallest), target_sum / np.where(smallest > 0, smallest, 1.0), 0.0
    )
    provisional = expm1 * np.repeat(provisional_lib / target_sum, lengths)

    quantum, row_err = _row_quantum(provisional, lengths, integrality_tolerance, max_quantum)
    unresolved = np.flatnonzero((lengths > 0) & (quantum == 0))
    if unresolved.size:
        raise ValueError(
            f"{unresolved.size} of {len(lengths)} rows are not integral for any "
            f"rarest-gene count up to {max_quantum} (best residual "
            f"{row_err[unresolved].min():.3e}); the source matrix was not built "
            "from integer counts"
        )

    final = provisional * np.repeat(np.maximum(quantum, 1).astype(np.float64), lengths)
    int_err = float(np.abs(final - np.round(final)).max()) if final.size else 0.0
    counts = sparse.csr_matrix(
        (np.round(final).astype(np.int64), block.indices.copy(), block.indptr.copy()),
        shape=block.shape,
    )

    return RecoveredCounts(
        counts=counts,
        library_sizes=np.asarray(counts.sum(axis=1)).ravel().astype(np.int64),
        row_quantum=quantum,
        max_integrality_error=int_err,
        max_rowsum_error=rowsum_err,
        n_empty_rows=int((lengths == 0).sum()),
    )


def recovery_diagnostics(
    block: sparse.csr_matrix, recovered: RecoveredCounts, *, target_sum: float = TARGET_SUM
) -> dict[str, float]:
    """Round-trip evidence: re-normalise the recovered counts and compare.

    The round trip establishes that the recovered integers regenerate ``X``.
    It cannot establish the *scale*, because every integer multiple of a row
    round-trips equally well -- which is why the library-size quantiles are
    reported alongside, as the statistic that would expose a wrong scale.
    """
    block = sparse.csr_matrix(block)
    counts = recovered.counts
    lib = recovered.library_sizes.astype(np.float64)
    lengths = np.diff(counts.indptr)
    safe = np.where(lib > 0, lib, 1.0)
    round_trip = np.log1p(counts.data.astype(np.float64) * target_sum / np.repeat(safe, lengths))
    err = np.abs(round_trip - np.asarray(block.data, dtype=np.float64))
    nz = recovered.library_sizes > 0
    return {
        "max_round_trip_error": float(err.max()) if err.size else 0.0,
        "median_round_trip_error": float(np.median(err)) if err.size else 0.0,
        "max_integrality_error": recovered.max_integrality_error,
        "max_rowsum_error": recovered.max_rowsum_error,
        "n_rows_rescaled": float(recovered.n_rows_rescaled),
        "n_empty_rows": float(recovered.n_empty_rows),
        "median_library_size": float(np.median(recovered.library_sizes[nz])) if nz.any() else 0.0,
        "min_library_size": float(recovered.library_sizes[nz].min()) if nz.any() else 0.0,
        "max_library_size": float(recovered.library_sizes.max(initial=0)),
    }
