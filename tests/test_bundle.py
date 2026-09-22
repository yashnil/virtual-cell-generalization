"""Response-to-counts translation, the unsupported-gene rule, and the writer.

The unsupported-gene rule is the one that most easily goes wrong silently, so
it is pinned from both sides: the lifted fold change must be zero there, and a
generator driven by it must reproduce the control distribution rather than
emit zeros.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from virtual_cell.arc import bundle, generate
from virtual_cell.arc.panel import Support, build_panel_map


def profile(counts: np.ndarray, target_sum: float = bundle.TARGET_SUM) -> np.ndarray:
    """The mean-``log1p(CP10K)`` profile the public responses are defined in."""
    x = np.asarray(counts, dtype=np.float64)
    lib = x.sum(axis=1, keepdims=True)
    return np.log1p(target_sum * np.divide(x, lib, out=np.zeros_like(x), where=lib > 0)).mean(
        axis=0
    )


# --------------------------------------------------------------------------
# the fold-change convention
# --------------------------------------------------------------------------


def test_zero_delta_is_zero_fold_change():
    basal = np.array([0.0, 0.5, 2.0, 7.0])
    assert np.allclose(bundle.log2_fold_change(basal, np.zeros_like(basal)), 0.0)


def test_fold_change_is_the_ratio_of_profiles_not_two_to_the_delta():
    """``2**delta`` is the wrong answer and this states why the code differs."""
    basal = np.array([4.0])
    delta = np.array([1.0])
    got = bundle.log2_fold_change(basal, delta)
    naive = delta / np.log2(np.e)
    expected = np.log2(
        (np.expm1(5.0) + bundle.FLOOR_FRACTION * bundle.TARGET_SUM)
        / (np.expm1(4.0) + bundle.FLOOR_FRACTION * bundle.TARGET_SUM)
    )
    assert got[0] == pytest.approx(expected)
    assert abs(got[0] - naive[0]) > 0.1


def test_positive_delta_raises_and_negative_lowers():
    basal = np.full(5, 3.0)
    assert np.all(bundle.log2_fold_change(basal, np.full(5, 0.5)) > 0)
    assert np.all(bundle.log2_fold_change(basal, np.full(5, -0.5)) < 0)


def test_a_gene_the_controls_never_detected_cannot_blow_up():
    basal = np.zeros(3)
    out = bundle.log2_fold_change(basal, np.array([0.0, 5.0, -5.0]))
    assert np.all(np.abs(out) <= bundle.LFC_CLIP)
    assert np.all(np.isfinite(out))


def test_fold_change_is_clipped_to_the_stated_range():
    basal = np.array([1e-9, 20.0])
    out = bundle.log2_fold_change(basal, np.array([30.0, -30.0]), clip=2.0)
    assert np.allclose(np.abs(out), 2.0)


def test_fold_change_broadcasts_over_perturbations():
    basal = np.linspace(0.1, 4.0, 6)
    delta = np.random.default_rng(0).normal(scale=0.2, size=(5, 6))
    out = bundle.log2_fold_change(basal, delta)
    assert out.shape == (5, 6)
    assert np.allclose(out[2], bundle.log2_fold_change(basal, delta[2]))


def test_fold_change_rejects_a_mismatched_axis():
    with pytest.raises(ValueError, match="genes"):
        bundle.log2_fold_change(np.zeros(4), np.zeros(5))


def test_predicted_profile_renormalises_to_the_target_sum():
    basal = np.linspace(0.0, 5.0, 40)
    out = bundle.predicted_profile(basal, np.full(40, 0.7))
    assert out.sum() == pytest.approx(bundle.TARGET_SUM)


def test_predicted_profile_of_zero_delta_is_the_control_composition():
    rng = np.random.default_rng(1)
    counts = rng.poisson(3.0, size=(50, 30))
    basal = profile(counts)
    out = bundle.predicted_profile(basal, np.zeros(30))
    assert np.allclose(out / out.sum(), np.expm1(basal) / np.expm1(basal).sum())


# --------------------------------------------------------------------------
# the unsupported-gene rule
# --------------------------------------------------------------------------


def panel_fixture() -> tuple[object, np.ndarray]:
    panel = [f"G{i}" for i in range(8)]
    source = ["G1", "G3", "G5"]
    panel_map = build_panel_map(panel, source)
    basal = np.linspace(1.0, 4.0, 8)
    return panel_map, basal


def test_unsupported_panel_genes_get_exactly_zero_fold_change():
    panel_map, basal = panel_fixture()
    delta = np.array([[0.9, -0.9, 0.4]])
    lfc = bundle.panel_log2_fold_change(panel_map, basal, delta)
    unsupported = ~panel_map.mask(Support.PREDICTED)
    assert np.array_equal(lfc[0][unsupported], np.zeros(int(unsupported.sum())))
    assert np.all(lfc[0][panel_map.mask(Support.PREDICTED)] != 0.0)


def test_supported_genes_land_in_the_right_panel_columns():
    panel_map, basal = panel_fixture()
    delta = np.array([[1.0, 2.0, 3.0]])
    lfc = bundle.panel_log2_fold_change(panel_map, basal, delta)
    for col, value in zip([1, 3, 5], [1.0, 2.0, 3.0], strict=True):
        assert lfc[0, col] == pytest.approx(
            bundle.log2_fold_change(basal[col : col + 1], np.array([value]))[0]
        )


def test_zero_fold_change_reproduces_the_control_distribution_not_zeros():
    """The rule is 'keep the control distribution', which is not zero-filling."""
    rng = np.random.default_rng(2)
    controls = rng.poisson(2.0, size=(400, 60)).astype(np.float64)
    out = generate.transport_controls(controls, np.zeros(60), 400, rng=rng)
    assert out.sum() > 0
    detected = (out > 0).sum(axis=1)
    reference = (controls > 0).sum(axis=1)
    assert np.median(detected) == pytest.approx(np.median(reference), rel=0.1)


def test_panel_fold_change_rejects_a_wrong_sized_control_profile():
    panel_map, _ = panel_fixture()
    with pytest.raises(ValueError, match="control profile"):
        bundle.panel_log2_fold_change(panel_map, np.zeros(3), np.zeros((1, 3)))


def test_panel_fold_change_rejects_a_wrong_sized_response():
    panel_map, basal = panel_fixture()
    with pytest.raises(ValueError, match="source space"):
        bundle.panel_log2_fold_change(panel_map, basal, np.zeros((1, 4)))


# --------------------------------------------------------------------------
# the writer
# --------------------------------------------------------------------------


def write_small(tmp_path, n_cells: int = 6, seed: int = 0):
    import anndata

    rng = np.random.default_rng(seed)
    genes = [f"G{i}" for i in range(10)]
    path = tmp_path / "submission.h5ad"
    blocks = {}
    with bundle.SubmissionWriter(path, genes) as writer:
        for context in ("A", "B"):
            for pert in ("P1", "P2"):
                block = rng.poisson(1.2, size=(n_cells, 10)).astype(np.int64)
                blocks[(context, pert)] = block
                writer.append(block, target_gene=pert, context=context)
    return path, genes, blocks, anndata.read_h5ad(path)


def test_writer_round_trips_through_anndata(tmp_path):
    path, genes, blocks, adata = write_small(tmp_path)
    assert adata.shape == (24, 10)
    assert list(adata.var_names) == genes
    stacked = np.vstack([blocks[k] for k in [("A", "P1"), ("A", "P2"), ("B", "P1"), ("B", "P2")]])
    assert np.array_equal(np.asarray(adata.X.todense()), stacked)


def test_writer_preserves_labels_in_append_order(tmp_path):
    path, _, _, adata = write_small(tmp_path)
    assert list(adata.obs["target_gene"][:6]) == ["P1"] * 6
    assert list(adata.obs["context"][:12]) == ["A"] * 12
    assert list(adata.obs["context"][12:]) == ["B"] * 12


def test_writer_uses_int64_offsets_and_int32_columns(tmp_path):
    """``indptr`` is cumulative and can exceed 2**31; ``indices`` cannot.

    Writing both as int64 would add eight gigabytes to a full submission for
    no information, and writing both as int32 would overflow the offsets.
    """
    import h5py

    path, _, _, _ = write_small(tmp_path)
    with h5py.File(path, "r") as f:
        assert f["X/indptr"].dtype == np.int64
        assert f["X/indices"].dtype == np.int32


def test_writer_rejects_negative_counts(tmp_path):
    with bundle.SubmissionWriter(tmp_path / "x.h5ad", ["A", "B"]) as writer:
        with pytest.raises(ValueError, match="non-negative"):
            writer.append(np.array([[1, -1]]), target_gene="P", context="A")


def test_writer_rejects_a_wrong_gene_count(tmp_path):
    with bundle.SubmissionWriter(tmp_path / "x.h5ad", ["A", "B"]) as writer:
        with pytest.raises(ValueError, match="block must be"):
            writer.append(np.zeros((2, 3)), target_gene="P", context="A")


def test_writer_rejects_duplicate_gene_names(tmp_path):
    with pytest.raises(ValueError, match="unique"):
        bundle.SubmissionWriter(tmp_path / "x.h5ad", ["A", "A"])


def test_writer_handles_an_all_zero_block(tmp_path):
    import anndata

    path = tmp_path / "z.h5ad"
    with bundle.SubmissionWriter(path, ["A", "B", "C"]) as writer:
        writer.append(np.zeros((3, 3)), target_gene="P", context="A")
        writer.append(np.ones((2, 3)), target_gene="Q", context="A")
    adata = anndata.read_h5ad(path)
    assert adata.shape == (5, 3)
    assert adata.X.sum() == 6


# --------------------------------------------------------------------------
# reading back
# --------------------------------------------------------------------------


def test_inspect_reports_the_written_shape_and_labels(tmp_path):
    path, genes, _, _ = write_small(tmp_path)
    report = bundle.inspect_bundle(path, genes)
    assert report.n_cells == 24
    assert report.n_genes == 10
    assert report.gene_order_matches
    assert report.integer_valued and report.non_negative and report.finite
    assert report.contexts == ("A", "B")
    assert report.n_perturbations_per_context == {"A": 2, "B": 2}


def test_inspect_detects_a_wrong_gene_order(tmp_path):
    path, genes, _, _ = write_small(tmp_path)
    assert not bundle.inspect_bundle(path, list(reversed(genes))).gene_order_matches


def test_inspect_counts_cells_per_group(tmp_path):
    path, genes, _, _ = write_small(tmp_path)
    report = bundle.inspect_bundle(path, genes)
    assert set(report.cells_per_group.unique()) == {6}


def test_inspect_flags_a_non_integer_matrix(tmp_path):
    import anndata

    path = tmp_path / "frac.h5ad"
    adata = anndata.AnnData(
        X=sparse.csr_matrix(np.array([[0.5, 1.0], [2.0, 0.0]], dtype=np.float32)),
        obs=pd.DataFrame(
            {"target_gene": pd.Categorical(["P", "P"]), "context": pd.Categorical(["A", "A"])},
            index=["c0", "c1"],
        ),
        var=pd.DataFrame(index=pd.Index(["G0", "G1"])),
    )
    adata.write_h5ad(path)
    assert not bundle.inspect_bundle(path, ["G0", "G1"]).integer_valued


def test_inspect_reports_the_true_minimum_cell_total(tmp_path):
    """``min(initial=0)`` would report 0 here and hide the real floor."""
    path = tmp_path / "m.h5ad"
    with bundle.SubmissionWriter(path, ["G0", "G1", "G2"]) as writer:
        writer.append(np.array([[5, 0, 2], [1, 1, 1]]), target_gene="P", context="A")
    report = bundle.inspect_bundle(path, ["G0", "G1", "G2"])
    assert report.min_counts_per_cell == 3
    assert report.max_counts_per_cell == 7


def test_inspect_still_reports_a_genuinely_empty_cell(tmp_path):
    path = tmp_path / "e.h5ad"
    with bundle.SubmissionWriter(path, ["G0", "G1"]) as writer:
        writer.append(np.array([[0, 0], [4, 1]]), target_gene="P", context="A")
    assert bundle.inspect_bundle(path, ["G0", "G1"]).min_counts_per_cell == 0
