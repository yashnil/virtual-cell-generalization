from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from virtual_cell.data.io import (
    DataIntegrityError,
    check_integrity,
    load_context,
    load_contexts,
)

REPO_SYNTHETIC = Path(__file__).resolve().parents[1] / "data" / "raw" / "synthetic"


def _adata(counts, genes=None, cells=None, context="A"):
    counts = np.asarray(counts)
    n_cells, n_genes = counts.shape
    genes = genes or [f"G{i}" for i in range(n_genes)]
    cells = cells or [f"c{i}" for i in range(n_cells)]
    obs = pd.DataFrame({"context": context}, index=cells)
    var = pd.DataFrame(index=genes)
    return ad.AnnData(X=sparse.csr_matrix(counts), obs=obs, var=var)


def test_load_synthetic_contexts(synthetic_paths):
    for label, path in synthetic_paths.items():
        adata = load_context(path, expected_context=label)
        assert adata.n_obs == 50
        assert adata.n_vars == 30
        assert sparse.issparse(adata.X)


def test_load_contexts_from_mapping_and_iterable(synthetic_paths):
    by_label = load_contexts(synthetic_paths)
    by_path = load_contexts(list(synthetic_paths.values()))
    assert set(by_label) == {"A", "B", "C"}
    assert set(by_path) == {"A", "B", "C"}


def test_counts_are_nonnegative_integers(synthetic_contexts):
    for adata in synthetic_contexts.values():
        values = adata.X.data
        assert values.min() >= 0
        assert np.issubdtype(values.dtype, np.integer)
        check_integrity(adata)  # should not raise


def test_context_label_preserved_exactly(synthetic_paths):
    for label, path in synthetic_paths.items():
        adata = load_context(path)
        assert adata.obs["context"].unique().tolist() == [label]
        assert adata.obs["context"].dtype.name in {"category", "object"}
        assert str(adata.obs["context"].iloc[0]) == label
    with pytest.raises(DataIntegrityError, match="mismatch"):
        load_context(synthetic_paths["A"], expected_context="B")


def test_context_label_case_is_not_normalised(tmp_path):
    adata = _adata(np.ones((3, 2)), context="k562")
    path = tmp_path / "k.h5ad"
    adata.write_h5ad(path)
    with pytest.raises(DataIntegrityError):
        load_context(path, expected_context="K562")
    assert load_context(path, expected_context="k562").obs["context"].iloc[0] == "k562"


def test_duplicate_genes_rejected(tmp_path):
    adata = _adata(np.ones((3, 3)), genes=["G1", "G1", "G2"])
    path = tmp_path / "dup.h5ad"
    adata.write_h5ad(path)
    with pytest.raises(DataIntegrityError, match="Duplicate gene"):
        load_context(path)


def test_duplicate_cells_rejected():
    adata = _adata(np.ones((2, 2)), cells=["c", "c"])
    with pytest.raises(DataIntegrityError, match="Duplicate cell"):
        check_integrity(adata)


def test_negative_counts_rejected():
    adata = _adata([[1, -1], [0, 2]])
    with pytest.raises(DataIntegrityError, match="negative"):
        check_integrity(adata)


def test_non_integer_counts_rejected_unless_allowed():
    adata = _adata([[0.5, 1.0], [0.0, 2.0]])
    with pytest.raises(DataIntegrityError, match="non-integer"):
        check_integrity(adata)
    check_integrity(adata, require_integer=False)


def test_nan_counts_rejected():
    adata = _adata([[np.nan, 1.0], [0.0, 2.0]])
    with pytest.raises(DataIntegrityError, match="NaN"):
        check_integrity(adata)


def test_missing_context_column_rejected():
    adata = _adata(np.ones((2, 2)))
    adata.obs = adata.obs.drop(columns=["context"])
    with pytest.raises(DataIntegrityError, match="missing the context column"):
        check_integrity(adata)
    check_integrity(adata, context_key=None)  # label checks can be skipped


def test_multiple_context_labels_rejected():
    adata = _adata(np.ones((2, 2)))
    adata.obs["context"] = ["A", "B"]
    with pytest.raises(DataIntegrityError, match="exactly one context"):
        check_integrity(adata)


def test_empty_dataset_rejected():
    adata = _adata(np.ones((2, 2)))
    with pytest.raises(DataIntegrityError, match="no cells"):
        check_integrity(adata[:0].copy())
    with pytest.raises(DataIntegrityError, match="no genes"):
        check_integrity(adata[:, :0].copy())


def test_missing_file_and_wrong_suffix(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_context(tmp_path / "missing.h5ad")
    (tmp_path / "x.csv").write_text("a")
    with pytest.raises(ValueError, match=".h5ad"):
        load_context(tmp_path / "x.csv")


def test_duplicate_context_label_across_files(tmp_path):
    a = _adata(np.ones((2, 2)), context="A")
    p1, p2 = tmp_path / "a1.h5ad", tmp_path / "a2.h5ad"
    a.write_h5ad(p1)
    a.write_h5ad(p2)
    with pytest.raises(DataIntegrityError, match="more than one file"):
        load_contexts([p1, p2])


@pytest.mark.skipif(
    not (REPO_SYNTHETIC / "context_A.h5ad").exists(),
    reason="on-disk synthetic data not generated",
)
def test_repo_synthetic_files_load():
    contexts = load_contexts(sorted(REPO_SYNTHETIC.glob("context_*.h5ad")))
    assert set(contexts) == {"A", "B", "C"}
    for adata in contexts.values():
        assert adata.shape == (200, 100)
