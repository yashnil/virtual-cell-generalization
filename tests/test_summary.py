import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from virtual_cell.data.summary import (
    ContextSummary,
    genes_detected_per_cell,
    library_sizes,
    sparsity,
    summarize_context,
    summarize_contexts,
)


@pytest.fixture
def tiny():
    counts = np.array(
        [
            [0, 3, 0, 1],
            [2, 0, 0, 0],
            [0, 0, 0, 0],
            [5, 5, 5, 5],
        ]
    )
    obs = pd.DataFrame({"context": "T"}, index=[f"c{i}" for i in range(4)])
    var = pd.DataFrame(index=[f"g{i}" for i in range(4)])
    return ad.AnnData(X=sparse.csr_matrix(counts), obs=obs, var=var), counts


def test_library_sizes_and_sparsity_match_numpy(tiny):
    adata, counts = tiny
    np.testing.assert_allclose(library_sizes(adata), counts.sum(axis=1))
    assert sparsity(adata) == pytest.approx(1 - np.count_nonzero(counts) / counts.size)
    np.testing.assert_array_equal(genes_detected_per_cell(adata), np.count_nonzero(counts, axis=1))


def test_summarize_context_values(tiny):
    adata, counts = tiny
    s = summarize_context(adata)
    assert isinstance(s, ContextSummary)
    assert s.context == "T"
    assert s.n_cells == 4
    assert s.n_genes == 4
    lib = counts.sum(axis=1)  # [4, 2, 0, 20]
    assert s.mean_library_size == pytest.approx(lib.mean())
    assert s.median_library_size == pytest.approx(np.median(lib))
    assert s.min_library_size == 0
    assert s.max_library_size == 20
    assert s.sparsity == pytest.approx(9 / 16)
    assert s.mean_genes_detected == pytest.approx((2 + 1 + 0 + 4) / 4)


def test_dense_and_sparse_agree(tiny):
    adata, counts = tiny
    dense = ad.AnnData(X=counts.astype(float), obs=adata.obs.copy(), var=adata.var.copy())
    assert summarize_context(dense) == summarize_context(adata)


def test_explicit_zero_in_sparse_data_not_counted():
    X = sparse.csr_matrix(np.array([[1.0, 0.0], [0.0, 2.0]]))
    X.data[0] = 0.0  # stored explicit zero
    adata = ad.AnnData(X=X, obs=pd.DataFrame(index=["a", "b"]), var=pd.DataFrame(index=["x", "y"]))
    assert sparsity(adata) == pytest.approx(0.75)


def test_summarize_synthetic_contexts(synthetic_contexts):
    df = summarize_contexts(synthetic_contexts)
    assert list(df.index) == ["A", "B", "C"]
    assert (df["n_cells"] == 50).all()
    assert (df["n_genes"] == 30).all()
    assert ((df["sparsity"] >= 0) & (df["sparsity"] <= 1)).all()
    # contexts are scaled by 1.0, 1.2, 1.4, so library sizes should be ordered
    # on average (with independent base rates this is only expected loosely,
    # so we check the mean library size is positive and finite).
    assert np.isfinite(df["mean_library_size"]).all()
    assert (df["mean_library_size"] > 0).all()


def test_unknown_context_label():
    adata = ad.AnnData(
        X=np.ones((2, 2)), obs=pd.DataFrame(index=["a", "b"]), var=pd.DataFrame(index=["x", "y"])
    )
    assert summarize_context(adata).context == "unknown"
