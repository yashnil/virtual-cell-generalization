"""C4 KOLF streaming statistics: equal to a direct computation and to the C3 protocol."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from virtual_cell.competition_v2 import fusion_c3, licensing, sources, sources_c4


def _write(tmp_path, rng):
    n_genes = 30
    genes = [f"G{i}" for i in range(n_genes - 2)] + ["G0", "G1"]  # duplicated symbols
    targets = ["T0", "T1", "T2", "T3"]
    labels = np.array(["NTC"] * 300 + ["T0"] * 90 + ["T1"] * 60 + ["T2"] * 30 + ["T3"] * 50)
    rng.shuffle(labels)
    batch = rng.choice(["ALPHA", "BETA", "GAMMA"], size=len(labels))
    mean = rng.gamma(1.0, 5.0, size=n_genes)
    lam = mean[None, :] * rng.gamma(2.0, 0.5, size=(len(labels), n_genes))
    lam[labels == "T0", 3] *= 0.2
    x = rng.poisson(lam).astype(np.float32)
    x[:, 5] += 1  # no zero-depth cells
    obs = pd.DataFrame(
        {
            "gene_target": pd.Categorical(labels),
            "batch": pd.Categorical(batch),
            "total_counts": x.sum(1).astype(np.float32),
        },
        index=[f"c{i}" for i in range(len(labels))],
    )
    a = ad.AnnData(sparse.csc_matrix(x), obs=obs, var=pd.DataFrame(index=pd.Index(genes)))
    path = tmp_path / "kolf.h5ad"
    a.write_h5ad(path)
    return path, x, labels, batch, genes, targets


def test_streaming_statistics_match_direct(tmp_path):
    rng = np.random.default_rng(0)
    path, x, labels, batch, genes, targets = _write(tmp_path, rng)
    stats, halves = sources_c4.kolf_statistics(
        path, targets, seed=1, half_seed=7, null_n=0, block_genes=7, log=lambda m: None
    )
    ugenes, proj = sources.symbol_projection(np.array(genes))
    xp = sparse.csr_matrix(x.astype(np.float64)) @ proj
    np.testing.assert_array_equal(stats["genes"], ugenes)
    for i, t in enumerate(targets):
        rows = labels == t
        np.testing.assert_allclose(
            stats["target_count_sums"][i], np.asarray(xp[rows].sum(0)).ravel(), rtol=1e-6
        )
        cpm = xp[rows].toarray() / x[rows].sum(1, keepdims=True) * 1e6
        np.testing.assert_allclose(stats["target_mean_cpm"][i], cpm.mean(0), rtol=1e-5)
        assert stats["n_cells"][i] == rows.sum()
    # reliability equals the C3 split-half protocol on the same cells and seed
    kept = sorted(t for t in targets if (labels == t).sum() >= 40)
    want = fusion_c3.split_half(xp, labels, kept, "NTC", seed=7)
    got = sources_c4.reliability_from_halves(
        halves["mean_a"], halves["mean_b"], halves["control_mean_cpm"], list(halves["kept"])
    )
    assert list(got["rel"]) == list(want["rel"]) == ["T0", "T1", "T3"]
    for t in want["rel"]:
        assert got["rel"][t] == pytest.approx(want["rel"][t], rel=1e-9, abs=1e-12)
        assert got["noise"][t] == pytest.approx(want["noise"][t], rel=1e-9)
    np.testing.assert_array_equal(got["genes"], want["genes"])


def test_loads_as_a_source(tmp_path):
    rng = np.random.default_rng(1)
    path, *_, targets = _write(tmp_path, rng)
    stats, _ = sources_c4.kolf_statistics(
        path, targets, seed=1, half_seed=7, null_n=2, null_size=40, log=lambda m: None
    )
    np.savez(tmp_path / "s.npz", **stats)
    src = sources.load_source(tmp_path / "s.npz")
    assert src.name == sources_c4.KOLF_NAME
    assert set(src.usable()) == {"T0", "T1", "T2", "T3"}
    assert licensing.STATUS[sources_c4.KOLF_NAME] == licensing.GREEN
    assert licensing.STATUS["JURKAT_GSE249595"] == licensing.UNKNOWN
