"""Invariants of the scPertEval four-context bundle and its pseudobulk pipeline.

Split into two groups:

* Tests on a small synthetic ``.h5ad`` built in-process, which always run and
  pin the pseudobulk mathematics.
* Tests against the real downloaded bundle, skipped when the git-ignored data
  is absent, which pin the audited facts and the frozen design.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from virtual_cell.data import scperteval
from virtual_cell.data.io import DataIntegrityError

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / scperteval.DATA_SUBDIR
SPLITS_DIR = REPO_ROOT / "data" / "splits" / "four_context_v1"

HAVE_DATA = not scperteval.incomplete_files(DATA_DIR)
needs_data = pytest.mark.skipif(
    not HAVE_DATA, reason=f"scPertEval bundle absent or incomplete in {DATA_DIR}"
)


# --- registry (no data needed) --------------------------------------------


def test_registry_covers_four_distinct_contexts():
    assert len(scperteval.DATASETS) == 4
    assert len({d.name for d in scperteval.DATASETS}) == 4
    assert len({d.cell_line for d in scperteval.DATASETS}) == 4
    assert scperteval.CONTEXTS == tuple(d.name for d in scperteval.DATASETS)


def test_registry_urls_and_filenames_are_well_formed():
    for d in scperteval.DATASETS:
        assert d.filename == f"{d.name}_processed_complete.h5ad"
        assert d.url.startswith("https://storage.googleapis.com/scperteval/processed/")
        assert d.expected_bytes > 0
        assert d.md5_base64.endswith("==")


def test_unknown_dataset_is_rejected():
    with pytest.raises(KeyError, match="Unknown scPertEval dataset"):
        scperteval.dataset("not_a_dataset")


# --- intersection semantics (no data needed) -------------------------------


def test_shared_perturbations_excludes_the_control_label():
    labels = {
        "a": ["control", "G1", "G2", "G3"],
        "b": ["control", "G2", "G3", "G4"],
    }
    assert scperteval.shared_perturbations(labels) == ["G2", "G3"]


def test_shared_sets_are_sorted_and_order_independent():
    labels = {"a": ["G3", "G1", "control"], "b": ["G1", "G3"]}
    reversed_labels = {"b": ["G3", "G1"], "a": ["control", "G1", "G3"]}
    assert scperteval.shared_perturbations(labels) == ["G1", "G3"]
    assert scperteval.shared_perturbations(labels) == scperteval.shared_perturbations(
        reversed_labels
    )


def test_empty_intersections_are_rejected():
    with pytest.raises(DataIntegrityError, match="share no perturbations"):
        scperteval.shared_perturbations({"a": ["G1"], "b": ["G2"]})
    with pytest.raises(DataIntegrityError, match="share no genes"):
        scperteval.shared_genes({"a": ["X"], "b": ["Y"]})
    with pytest.raises(DataIntegrityError, match="No contexts"):
        scperteval.shared_perturbations({})


# --- pseudobulk mathematics on a synthetic file ---------------------------


@pytest.fixture
def synthetic_h5ad(tmp_path) -> Path:
    """Tiny dataset with hand-checkable group means."""
    rng = np.random.default_rng(0)
    labels = ["control"] * 6 + ["GENE_A"] * 4 + ["GENE_B"] * 3 + ["GENE_C"] * 2
    genes = ["g0", "g1", "g2", "g3", "g4"]
    X = rng.integers(0, 5, size=(len(labels), len(genes))).astype(np.float32)
    adata = ad.AnnData(
        X=sparse.csr_matrix(X),
        obs=pd.DataFrame({"perturbation": pd.Categorical(labels)}),
        var=pd.DataFrame(index=genes),
    )
    path = tmp_path / "synthetic_processed_complete.h5ad"
    adata.write_h5ad(path)
    return path


def test_pseudobulk_matches_hand_computed_group_means(synthetic_h5ad):
    adata = ad.read_h5ad(synthetic_h5ad)
    X = np.asarray(adata.X.todense())
    labels = adata.obs["perturbation"].astype(str).to_numpy()
    genes = ["g1", "g3"]
    perts = ["GENE_A", "GENE_B"]

    pb = scperteval.pseudobulk(synthetic_h5ad, genes=genes, perturbations=perts)

    cols = [list(adata.var_names).index(g) for g in genes]
    np.testing.assert_allclose(pb.control_mean, X[labels == "control"][:, cols].mean(0), rtol=1e-6)
    for i, p in enumerate(perts):
        np.testing.assert_allclose(
            pb.perturbation_means[i], X[labels == p][:, cols].mean(0), rtol=1e-6
        )
    assert pb.cell_counts.tolist() == [4, 3]
    assert pb.n_control_cells == 6


def test_delta_is_perturbation_mean_minus_control_mean(synthetic_h5ad):
    pb = scperteval.pseudobulk(
        synthetic_h5ad, genes=["g0", "g2"], perturbations=["GENE_A", "GENE_C"]
    )
    np.testing.assert_allclose(
        pb.delta, pb.perturbation_means - pb.control_mean[None, :], rtol=0, atol=0
    )


def test_pseudobulk_respects_the_requested_gene_order(synthetic_h5ad):
    forward = scperteval.pseudobulk(synthetic_h5ad, genes=["g0", "g4"], perturbations=["GENE_A"])
    reverse = scperteval.pseudobulk(synthetic_h5ad, genes=["g4", "g0"], perturbations=["GENE_A"])
    np.testing.assert_allclose(forward.control_mean, reverse.control_mean[::-1])


def test_pseudobulk_ignores_cells_outside_the_frozen_perturbation_set(synthetic_h5ad):
    """GENE_C cells must not leak into GENE_A's mean or the control mean."""
    adata = ad.read_h5ad(synthetic_h5ad)
    X = np.asarray(adata.X.todense())
    labels = adata.obs["perturbation"].astype(str).to_numpy()
    pb = scperteval.pseudobulk(
        synthetic_h5ad, genes=list(adata.var_names), perturbations=["GENE_A"]
    )
    np.testing.assert_allclose(pb.perturbation_means[0], X[labels == "GENE_A"].mean(0), rtol=1e-6)
    np.testing.assert_allclose(pb.control_mean, X[labels == "control"].mean(0), rtol=1e-6)


def test_missing_gene_is_rejected(synthetic_h5ad):
    with pytest.raises(DataIntegrityError, match="genes absent"):
        scperteval.pseudobulk(synthetic_h5ad, genes=["g0", "nope"], perturbations=["GENE_A"])


def test_perturbation_with_no_cells_is_rejected(synthetic_h5ad):
    with pytest.raises(DataIntegrityError, match="no cells"):
        scperteval.pseudobulk(synthetic_h5ad, genes=["g0"], perturbations=["GENE_A", "ABSENT"])


def test_load_perturbation_cells_returns_the_right_cells(synthetic_h5ad):
    adata = ad.read_h5ad(synthetic_h5ad)
    labels = adata.obs["perturbation"].astype(str).to_numpy()
    X, group = scperteval.load_perturbation_cells(
        synthetic_h5ad, genes=list(adata.var_names), perturbations=["GENE_A", "GENE_B"]
    )
    assert X.shape[0] == int((labels == "GENE_A").sum() + (labels == "GENE_B").sum())
    assert sorted(np.bincount(group).tolist()) == [3, 4]
    # group means from the cell-level load must equal the streaming pseudobulk
    pb = scperteval.pseudobulk(
        synthetic_h5ad, genes=list(adata.var_names), perturbations=["GENE_A", "GENE_B"]
    )
    for g in (0, 1):
        np.testing.assert_allclose(
            np.asarray(X[group == g].todense()).mean(0), pb.perturbation_means[g], rtol=1e-5
        )


# --- the real bundle ------------------------------------------------------


@pytest.fixture(scope="module")
def structures() -> dict:
    return {d.name: scperteval.structure_report(d.path(DATA_DIR)) for d in scperteval.DATASETS}


@needs_data
def test_all_four_files_present_with_the_audited_sizes():
    for d in scperteval.DATASETS:
        assert d.path(DATA_DIR).stat().st_size == d.expected_bytes, d.name


@needs_data
def test_x_is_sparse_csr_float32(structures):
    for name, rep in structures.items():
        assert rep["x_encoding"] == "csr_matrix", name
        assert rep["x_dtype"] == "float32", name


@needs_data
def test_no_raw_layers_obsm_or_uns_surprises(structures):
    for name, rep in structures.items():
        assert rep["obs_columns"] == [scperteval.PERTURBATION_KEY], name
        assert rep["var_columns"] == [], name
        for key, value in rep["extras"].items():
            assert value == [], f"{name}: {key} = {value}"


@needs_data
def test_controls_present_and_labelled_exactly(structures):
    for name, rep in structures.items():
        assert rep["has_control"], name
        assert rep["n_control_cells"] > 0, name


@needs_data
def test_no_combination_perturbations(structures):
    for name, rep in structures.items():
        assert rep["n_combination_perturbations"] == 0, name


@needs_data
def test_min_cells_per_perturbation_floor_is_inherited(structures):
    """Verified from the files: the upstream floor is 30 in every context."""
    for name, rep in structures.items():
        assert rep["min_cells_per_pert"] >= 30, name


@needs_data
def test_var_names_unique_in_every_context(structures):
    for name, rep in structures.items():
        assert rep["var_names_unique"], name


@needs_data
def test_x_is_log1p_cp10k_not_raw_counts():
    """expm1(X) row sums must sit at or just below the 1e4 normalisation target.

    Genes were filtered after normalisation, so sums land slightly below 1e4;
    raw counts would be far away and would also be integer-valued.
    """
    for d in scperteval.DATASETS:
        lib = scperteval.library_size_check(d.path(DATA_DIR), n_cells=300)
        assert 0.5 <= lib["median_over_target"] <= 1.0001, (d.name, lib["median"])


@needs_data
def test_frozen_design_matches_a_fresh_recomputation():
    """The committed design files must equal the intersections recomputed now."""
    if not (SPLITS_DIR / "shared_perturbations.txt").exists():
        pytest.skip("Frozen design not generated yet.")
    labels, genes = {}, {}
    for d in scperteval.DATASETS:
        cats, _ = scperteval.read_perturbation_labels(d.path(DATA_DIR))
        labels[d.name] = cats
        genes[d.name] = list(scperteval.read_var_names(d.path(DATA_DIR)))

    frozen_p = (SPLITS_DIR / "shared_perturbations.txt").read_text().split()
    frozen_g = (SPLITS_DIR / "shared_genes.txt").read_text().split()
    frozen_c = (SPLITS_DIR / "contexts.txt").read_text().split()

    assert frozen_p == scperteval.shared_perturbations(labels)
    assert frozen_g == scperteval.shared_genes(genes)
    assert frozen_c == list(scperteval.CONTEXTS)


@needs_data
def test_frozen_design_is_sorted_and_deduplicated():
    if not (SPLITS_DIR / "shared_perturbations.txt").exists():
        pytest.skip("Frozen design not generated yet.")
    for name in ("shared_perturbations.txt", "shared_genes.txt"):
        values = (SPLITS_DIR / name).read_text().split()
        assert values == sorted(values), name
        assert len(values) == len(set(values)), name
