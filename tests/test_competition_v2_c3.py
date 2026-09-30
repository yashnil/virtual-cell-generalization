"""C3 fusion: C1a is reproduced bit for bit; weights redistribute, never shrink."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from virtual_cell.competition_v2 import fusion, fusion_c3, sources

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "outputs" / "competition_v2" / "c1_license_clean" / "sources"


@pytest.fixture(scope="module")
def real():
    need = [
        SRC / "K562_GWPS_CPM_full_statistics.npz",
        SRC / "H1_2025_full_statistics.npz",
        SRC / "CD4_DE_statistics.npz",
    ]
    if not all(p.exists() for p in need):
        pytest.skip("prepared C1 source statistics are not available")
    k562 = sources.load_source(need[0])
    h1 = sources.load_source(need[1])
    with np.load(need[2]) as d:
        cd4 = {k: d[k] for k in d.files}
    genes = np.intersect1d(k562.genes, h1.genes)[:600]
    targets = np.array(sorted(set(k562.targets[:40]) | set(h1.targets[:40])))
    rng = np.random.default_rng(0)
    controls = {s: rng.dirichlet(np.ones(len(genes))) for s in ("log2fc", "bulk_delta")}
    return {"K562": k562, "H1": h1}, cd4, targets, genes, controls


def test_equal_weights_reproduce_c1a_exactly(real):
    srcs, cd4, targets, genes, controls = real
    want = fusion.fused_effects(list(srcs.values()), [1.0, 1.0], cd4, 1.0, targets, genes, controls)
    comp = fusion_c3.components(srcs, cd4, targets, genes)
    got = fusion_c3.fused(comp, controls)
    for space in want:
        assert got[space].dtype == want[space].dtype
        np.testing.assert_array_equal(got[space], want[space])


def test_weights_do_not_change_single_source_cells(real):
    srcs, cd4, targets, genes, controls = real
    comp = fusion_c3.components(srcs, cd4, targets, genes)
    base = fusion_c3.combine(comp, "log2fc", None)
    moved = fusion_c3.combine(comp, "log2fc", None, weights={"K562": 3.0, "H1": 0.5, "CD4": 0.2})
    classes = fusion_c3.consensus_classes(comp)
    single = classes == 1
    np.testing.assert_allclose(moved[single], base[single], rtol=1e-6, atol=1e-7)
    assert not np.allclose(moved[classes >= 2], base[classes >= 2])


def test_per_target_weights_match_scalar_weights(real):
    srcs, cd4, targets, genes, controls = real
    comp = fusion_c3.components(srcs, cd4, targets, genes)
    n = len(targets)
    w = {"K562": 2.0, "H1": 1.0, "CD4": 0.5}
    a = fusion_c3.combine(comp, "bulk_delta", controls["bulk_delta"], weights=w)
    b = fusion_c3.combine(
        comp, "bulk_delta", controls["bulk_delta"], weights={k: np.full(n, v) for k, v in w.items()}
    )
    np.testing.assert_allclose(a, b, rtol=1e-5, atol=1e-6)


def test_scale_multiplies_single_source_cells(real):
    srcs, cd4, targets, genes, controls = real
    comp = fusion_c3.components(srcs, cd4, targets, genes)
    base = fusion_c3.combine(comp, "log2fc", None)
    scaled = fusion_c3.combine(comp, "log2fc", None, scales={"K562": 2.0, "H1": 2.0, "CD4": 2.0})
    np.testing.assert_allclose(scaled, 2 * base, rtol=1e-5, atol=1e-6)


def test_consensus_classes_synthetic():
    e1 = np.array([[1.0, 1.0, 1.0, 0.0]], dtype=np.float32)
    e2 = np.array([[2.0, -1.0, 0.0, 0.0]], dtype=np.float32)
    m = np.array([[True, True, True, False]])
    comp = {"log2fc": {"a": (e1, m), "b": (e2, m & (e2 != 0))}, "bulk_delta": {}}
    np.testing.assert_array_equal(fusion_c3.consensus_classes(comp)[0], [2, 3, 1, 0])
    np.testing.assert_array_equal(
        fusion_c3.consensus_factor(fusion_c3.consensus_classes(comp), 0.5)[0], [1, 0.5, 1, 1]
    )


def test_row_metrics_perfect_and_sign():
    t = np.array([[1.0, -2.0, 3.0, 0.5]])
    out, sse, sst = fusion_c3.row_metrics(2 * t, t, np.ones_like(t, dtype=bool), top_k=3)
    assert out["cosine"][0] == pytest.approx(1)
    assert out["norm_ratio"][0] == pytest.approx(2)
    assert out["sign_acc"][0] == 1
    assert sse == pytest.approx(sst)


def test_zero_weight_keeps_single_source_responses(real):
    srcs, cd4, targets, genes, controls = real
    extra = [t for t in cd4["targets"].astype(str) if t not in set(targets)][:20]
    targets = np.array(sorted(set(targets) | set(extra)))
    comp = fusion_c3.components(srcs, cd4, targets, genes)
    base = fusion_c3.combine(comp, "log2fc", None)
    zero = fusion_c3.combine(comp, "log2fc", None, weights={"CD4": 0.0})
    only_cd4 = comp["CD4"][1] & ~comp["log2fc"]["K562"][1] & ~comp["log2fc"]["H1"][1]
    assert only_cd4.any()
    np.testing.assert_allclose(zero[only_cd4], base[only_cd4], rtol=1e-6, atol=1e-7)
    per_target = fusion_c3.combine(comp, "log2fc", None, weights={"CD4": np.zeros(len(targets))})
    np.testing.assert_allclose(per_target[only_cd4], base[only_cd4], rtol=1e-6, atol=1e-7)
