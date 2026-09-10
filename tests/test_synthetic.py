import numpy as np

from virtual_cell.data.synthetic import make_synthetic_contexts, write_synthetic_contexts


def test_synthetic_is_deterministic():
    a = make_synthetic_contexts(("A", "B"), n_cells=20, n_genes=10, seed=3)
    b = make_synthetic_contexts(("A", "B"), n_cells=20, n_genes=10, seed=3)
    for k in a:
        np.testing.assert_array_equal(a[k].X.toarray(), b[k].X.toarray())


def test_synthetic_labels_and_shapes():
    ctx = make_synthetic_contexts(("K562", "rpe1"), n_cells=7, n_genes=5, seed=1)
    assert list(ctx) == ["K562", "rpe1"]
    for label, adata in ctx.items():
        assert adata.shape == (7, 5)
        assert adata.obs["context"].unique().tolist() == [label]
        assert adata.obs_names[0] == f"{label}_cell_0"
        assert not adata.var_names.has_duplicates


def test_write_synthetic_contexts_roundtrip(tmp_path):
    paths = write_synthetic_contexts(tmp_path, ("A",), n_cells=5, n_genes=4, seed=0)
    assert set(paths) == {"A"}
    assert paths["A"].name == "context_A.h5ad"
    assert paths["A"].exists()


def test_shared_base_rate_gives_correlated_contexts():
    from virtual_cell.preprocessing.pseudobulk import basal_mean_correlation, basal_mean_table

    shared = make_synthetic_contexts(("A", "B"), n_cells=100, n_genes=200, seed=0)
    indep = make_synthetic_contexts(
        ("A", "B"), n_cells=100, n_genes=200, seed=0, shared_base_rate=False
    )
    r_shared = basal_mean_correlation(basal_mean_table(shared)).loc["A", "B"]
    r_indep = basal_mean_correlation(basal_mean_table(indep)).loc["A", "B"]
    assert r_shared > 0.7
    assert r_shared > r_indep
