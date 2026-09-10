"""Pseudobulk (mean expression) utilities and cross-context comparisons.

The basal mean-expression profile of a context is the simplest possible
context representation (``z_c`` baseline in the project plan) and the reference
point for all perturbation-effect estimates, so it is worth getting right early.
"""

from __future__ import annotations

from collections.abc import Mapping

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from virtual_cell.data.summary import library_sizes


def normalized_matrix(
    adata: ad.AnnData, *, target_sum: float = 1e4, log1p: bool = True
) -> sparse.csr_matrix:
    """Return a library-size-normalised (and optionally log1p) copy of ``adata.X``.

    Each cell is scaled to ``target_sum`` total counts. Cells with zero counts
    are left as all-zero rows rather than producing NaNs.
    """
    X = sparse.csr_matrix(adata.X, dtype=float)
    lib = library_sizes(adata)
    scale = np.zeros_like(lib)
    nonzero = lib > 0
    scale[nonzero] = target_sum / lib[nonzero]
    X = sparse.diags(scale) @ X
    if log1p:
        X = X.log1p()
    return sparse.csr_matrix(X)


def mean_expression(
    adata: ad.AnnData,
    *,
    normalize: bool = True,
    target_sum: float = 1e4,
    log1p: bool = True,
) -> pd.Series:
    """Per-gene mean expression across all cells (pseudobulk profile).

    With ``normalize=False`` this is the mean of raw counts. With the default
    settings it is the mean of ``log1p(target_sum * x / library_size)``, the
    standard scanpy-style normalised expression.
    """
    X = normalized_matrix(adata, target_sum=target_sum, log1p=log1p) if normalize else adata.X
    means = np.asarray(X.mean(axis=0)).ravel()
    return pd.Series(means, index=adata.var_names, name="mean_expression")


def shared_genes(contexts: Mapping[str, ad.AnnData]) -> pd.Index:
    """Genes present in every context, in the order of the first context."""
    if not contexts:
        raise ValueError("No contexts provided.")
    iterator = iter(contexts.values())
    genes = pd.Index(next(iterator).var_names)
    for adata in iterator:
        genes = genes[genes.isin(adata.var_names)]
    if len(genes) == 0:
        raise ValueError("Contexts share no genes.")
    return genes


def basal_mean_table(
    contexts: Mapping[str, ad.AnnData],
    *,
    normalize: bool = True,
    target_sum: float = 1e4,
    log1p: bool = True,
) -> pd.DataFrame:
    """Genes x contexts table of basal mean expression on the shared gene set."""
    genes = shared_genes(contexts)
    columns = {
        label: mean_expression(
            adata[:, genes], normalize=normalize, target_sum=target_sum, log1p=log1p
        )
        for label, adata in contexts.items()
    }
    return pd.DataFrame(columns, index=genes)


def basal_mean_correlation(table: pd.DataFrame, method: str = "pearson") -> pd.DataFrame:
    """Pairwise correlation between context basal profiles."""
    return table.corr(method=method)


def basal_log_fold_change(table: pd.DataFrame, a: str, b: str, eps: float = 1e-6) -> pd.Series:
    """Per-gene ``log2((mean_a + eps) / (mean_b + eps))`` between two contexts."""
    return np.log2((table[a] + eps) / (table[b] + eps)).rename(f"log2fc_{a}_vs_{b}")
