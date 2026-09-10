import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from virtual_cell.preprocessing.pseudobulk import (
    basal_log_fold_change,
    basal_mean_correlation,
    basal_mean_table,
    mean_expression,
    normalized_matrix,
    shared_genes,
)


def _adata(counts, genes, context):
    counts = np.asarray(counts)
    obs = pd.DataFrame({"context": context}, index=[f"{context}{i}" for i in range(len(counts))])
    return ad.AnnData(X=sparse.csr_matrix(counts), obs=obs, var=pd.DataFrame(index=genes))


def test_normalized_matrix_rows_sum_to_target():
    adata = _adata([[1, 3], [0, 0], [5, 5]], ["g1", "g2"], "A")
    X = normalized_matrix(adata, target_sum=100, log1p=False).toarray()
    np.testing.assert_allclose(X.sum(axis=1), [100, 0, 100])
    Xlog = normalized_matrix(adata, target_sum=100, log1p=True).toarray()
    np.testing.assert_allclose(Xlog, np.log1p(X))


def test_mean_expression_raw_matches_numpy():
    counts = np.array([[1, 3], [0, 0], [5, 5]])
    adata = _adata(counts, ["g1", "g2"], "A")
    m = mean_expression(adata, normalize=False)
    np.testing.assert_allclose(m.to_numpy(), counts.mean(axis=0))
    assert list(m.index) == ["g1", "g2"]


def test_shared_genes_and_table_alignment():
    a = _adata([[1, 2, 3]], ["g1", "g2", "g3"], "A")
    b = _adata([[4, 5]], ["g3", "g1"], "B")
    genes = shared_genes({"A": a, "B": b})
    assert list(genes) == ["g1", "g3"]
    table = basal_mean_table({"A": a, "B": b}, normalize=False)
    assert list(table.index) == ["g1", "g3"]
    assert table.loc["g1", "A"] == 1 and table.loc["g3", "A"] == 3
    assert table.loc["g1", "B"] == 5 and table.loc["g3", "B"] == 4


def test_no_shared_genes_raises():
    a = _adata([[1]], ["g1"], "A")
    b = _adata([[1]], ["g2"], "B")
    with pytest.raises(ValueError, match="share no genes"):
        shared_genes({"A": a, "B": b})


def test_correlation_and_log_fold_change_on_synthetic(synthetic_contexts):
    table = basal_mean_table(synthetic_contexts)
    corr = basal_mean_correlation(table)
    assert corr.shape == (3, 3)
    np.testing.assert_allclose(np.diag(corr), 1.0)
    lfc = basal_log_fold_change(table, "A", "B")
    assert lfc.shape == (table.shape[0],)
    assert np.isfinite(lfc).all()
