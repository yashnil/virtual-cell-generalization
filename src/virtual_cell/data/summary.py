"""Per-context summary statistics for single-cell count matrices."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, dataclass

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from virtual_cell.data.io import DEFAULT_CONTEXT_KEY


def library_sizes(adata: ad.AnnData) -> np.ndarray:
    """Total counts per cell (row sums of ``adata.X``) as a 1-D float array."""
    X = adata.X
    if sparse.issparse(X):
        return np.asarray(X.sum(axis=1)).ravel().astype(float)
    return np.asarray(X).sum(axis=1).astype(float)


def sparsity(adata: ad.AnnData) -> float:
    """Fraction of entries in ``adata.X`` that are exactly zero."""
    X = adata.X
    n_total = adata.n_obs * adata.n_vars
    if n_total == 0:
        return float("nan")
    if sparse.issparse(X):
        n_nonzero = int(np.count_nonzero(X.data))
    else:
        n_nonzero = int(np.count_nonzero(np.asarray(X)))
    return 1.0 - n_nonzero / n_total


def genes_detected_per_cell(adata: ad.AnnData) -> np.ndarray:
    """Number of genes with non-zero counts in each cell."""
    X = adata.X
    if sparse.issparse(X):
        X = sparse.csr_matrix(X, copy=True)
        X.eliminate_zeros()  # explicit stored zeros must not count as detected
        return np.diff(X.indptr)
    return np.count_nonzero(np.asarray(X), axis=1)


@dataclass(frozen=True)
class ContextSummary:
    """Summary statistics describing one cellular context."""

    context: str
    n_cells: int
    n_genes: int
    mean_library_size: float
    median_library_size: float
    min_library_size: float
    max_library_size: float
    sparsity: float
    mean_genes_detected: float

    def to_dict(self) -> dict:
        return asdict(self)


def summarize_context(
    adata: ad.AnnData, *, context_key: str | None = DEFAULT_CONTEXT_KEY
) -> ContextSummary:
    """Compute :class:`ContextSummary` for one context.

    The context label is read from ``adata.obs[context_key]``; if the column is
    absent or ``context_key`` is ``None`` the label is reported as ``"unknown"``.
    """
    if context_key is not None and context_key in adata.obs.columns:
        label = str(adata.obs[context_key].iloc[0])
    else:
        label = "unknown"

    lib = library_sizes(adata)
    return ContextSummary(
        context=label,
        n_cells=int(adata.n_obs),
        n_genes=int(adata.n_vars),
        mean_library_size=float(lib.mean()),
        median_library_size=float(np.median(lib)),
        min_library_size=float(lib.min()),
        max_library_size=float(lib.max()),
        sparsity=float(sparsity(adata)),
        mean_genes_detected=float(genes_detected_per_cell(adata).mean()),
    )


def summarize_contexts(
    contexts: Mapping[str, ad.AnnData], *, context_key: str | None = DEFAULT_CONTEXT_KEY
) -> pd.DataFrame:
    """Summarise several contexts into one DataFrame (one row per context)."""
    rows = [summarize_context(a, context_key=context_key).to_dict() for a in contexts.values()]
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.set_index("context")
    return df
