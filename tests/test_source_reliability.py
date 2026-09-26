"""Block split-half machinery for the Kaden source-reliability diagnostic."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from virtual_cell.analysis import source_reliability as sr


def _groups(sizes, n_control):
    g = np.concatenate([np.full(n, i) for i, n in enumerate(sizes)] + [np.full(n_control, -2)])
    return np.random.default_rng(1).permutation(g)


def test_blocks_are_balanced_and_cover_every_cell():
    group = _groups([23, 40, 57], 101)
    layout = sr.assign_blocks(
        group,
        n_groups=3,
        control_code=-2,
        n_blocks=10,
        n_control_blocks=7,
        rng=np.random.default_rng(0),
    )
    assert (layout.cell_block >= 0).all()
    for g, n in enumerate([23, 40, 57]):
        sizes = np.bincount(layout.cell_block[group == g] - g * 10, minlength=10)
        assert sizes.sum() == n and sizes.max() - sizes.min() <= 1
    ctrl = np.bincount(layout.cell_block[group == -2] - layout.control_offset, minlength=7)
    assert ctrl.sum() == 101 and ctrl.max() - ctrl.min() <= 1


def test_ignored_cells_get_no_block():
    group = np.array([0, 0, -1, -2, -2, 0, -1])
    layout = sr.assign_blocks(
        group,
        n_groups=1,
        control_code=-2,
        n_blocks=2,
        n_control_blocks=2,
        rng=np.random.default_rng(0),
    )
    assert (layout.cell_block[group == -1] == -1).all()


def test_half_masks_split_exactly_in_half():
    m = sr.draw_half_masks(50, 20, np.random.default_rng(3))
    assert m.shape == (50, 20) and (m.sum(axis=1) == 10).all()
    with pytest.raises(ValueError):
        sr.draw_half_masks(2, 5, np.random.default_rng(0))


def test_half_means_equal_direct_cell_means():
    rng = np.random.default_rng(4)
    x = rng.poisson(2.0, size=(60, 8)).astype(float)
    group = np.repeat([0, 1], 30)
    layout = sr.assign_blocks(
        group, n_groups=2, control_code=-2, n_blocks=6, n_control_blocks=2, rng=rng
    )
    sums = np.zeros((layout.n_rows, 8))
    np.add.at(sums, layout.cell_block, x)
    counts = np.bincount(layout.cell_block, minlength=layout.n_rows)
    blocks = sums[: layout.control_offset].reshape(2, 6, 8)
    bcounts = counts[: layout.control_offset].reshape(2, 6)
    masks = sr.draw_half_masks(2, 6, rng)
    a, b, na, nb = sr.half_means(blocks, bcounts, masks)
    for g in range(2):
        chosen = np.flatnonzero(masks[g])
        in_a = np.isin(layout.cell_block - g * 6, chosen) & (group == g)
        in_b = ~np.isin(layout.cell_block - g * 6, chosen) & (group == g)
        np.testing.assert_allclose(a[g], x[in_a].mean(axis=0))
        np.testing.assert_allclose(b[g], x[in_b].mean(axis=0))
        assert na[g] + nb[g] == 30


def test_accumulate_block_sums_matches_dense(tmp_path):
    rng = np.random.default_rng(5)
    n_cells, genes = 90, [f"G{i}" for i in range(12)]
    x = rng.poisson(1.0, size=(n_cells, len(genes))).astype(np.float32)
    labels = np.array(["control"] * 30 + ["P1"] * 30 + ["P2"] * 30)
    adata = ad.AnnData(
        X=sparse.csr_matrix(x),
        obs=pd.DataFrame(
            {"perturbation": pd.Categorical(labels)}, index=[f"c{i}" for i in range(n_cells)]
        ),
        var=pd.DataFrame(index=genes),
    )
    path = tmp_path / "toy.h5ad"
    adata.write_h5ad(path)
    group = np.array([{"P1": 0, "P2": 1}.get(v, -2) for v in labels])
    layout = sr.assign_blocks(
        group, n_groups=2, control_code=-2, n_blocks=4, n_control_blocks=3, rng=rng
    )
    subset = ["G3", "G0", "G7"]
    sums, counts = sr.accumulate_block_sums(path, genes=subset, layout=layout, chunk_size=17)
    want = np.zeros_like(sums)
    np.add.at(want, layout.cell_block, x[:, [3, 0, 7]])
    np.testing.assert_allclose(sums, want, rtol=1e-6)
    assert counts.sum() == n_cells


def test_split_half_reliability_recovers_planted_signal():
    """A latent response plus independent cell noise: block split-halves give the
    analytic half reliability, and pure noise gives ~0."""
    rng = np.random.default_rng(6)
    n_genes, n_cells, sigma = 400, 200, 1.0
    latent = rng.normal(0, 0.2, n_genes)
    cells = latent + rng.normal(0, sigma, (n_cells, n_genes))
    group = np.zeros(n_cells, dtype=int)
    layout = sr.assign_blocks(
        group, n_groups=1, control_code=-2, n_blocks=20, n_control_blocks=2, rng=rng
    )
    sums = np.zeros((layout.n_rows, n_genes))
    np.add.at(sums, layout.cell_block, cells)
    counts = np.bincount(layout.cell_block, minlength=layout.n_rows)
    blocks = sums[:20].reshape(1, 20, n_genes)
    bc = counts[:20].reshape(1, 20)
    rhos = []
    for _ in range(50):
        a, b, _, _ = sr.half_means(blocks, bc, sr.draw_half_masks(1, 20, rng))
        rhos.append(sr.row_pearson(a, b)[0])
    var_l = latent.var()
    expected = var_l / (var_l + sigma**2 / (n_cells / 2))
    assert abs(np.mean(rhos) - expected) < 0.05


def test_main_effect_is_more_reliable_than_its_perturbations():
    """Averaging many noisy perturbations sharing a common shift yields a
    reliable main effect even when each perturbation alone is unreliable."""
    rng = np.random.default_rng(7)
    shared = rng.normal(0, 0.1, 300)
    halves = [shared + rng.normal(0, 1.0, (200, 300)) for _ in range(2)]
    per = sr.row_pearson(halves[0], halves[1])
    main = sr.pearson(halves[0].mean(axis=0), halves[1].mean(axis=0))
    assert np.median(per) < 0.05 and main > 0.5


def test_agreement_statistics():
    a = np.array([[1.0, 2.0, 3.0, -4.0]])
    assert sr.row_pearson(a, 2 * a)[0] == pytest.approx(1.0)
    assert sr.row_cosine(a, -a)[0] == pytest.approx(-1.0)
    assert np.isnan(sr.row_pearson(np.ones((1, 4)), a)[0])
    assert sr.norm_ratio(2 * a[0], a[0]) == pytest.approx(2.0)
    assert sr.sign_agreement_top(a, a, 2)[0] == 1.0
    assert sr.sign_agreement_top(a, -a, 2)[0] == 0.0


def test_noise_ceiling_and_disattenuation():
    assert sr.noise_ceiling(0.64, 0.25) == pytest.approx(0.4)
    assert sr.noise_ceiling(-0.2, 0.9) == 0.0
    d = sr.disattenuated(
        np.array([0.2, 0.2]), np.array([0.64, 0.05]), np.array([0.25, 0.9]), min_reliability=0.1
    )
    assert d[0] == pytest.approx(0.5) and np.isnan(d[1])


def test_quality_band_thresholds():
    kw = {"high": 0.5, "moderate": 0.2}
    assert sr.quality_band(0.5, **kw) == "high"
    assert sr.quality_band(0.2, **kw) == "moderate"
    assert sr.quality_band(0.19, **kw) == "low / uninformative"
    assert sr.quality_band(float("nan"), **kw) == "not estimable"
