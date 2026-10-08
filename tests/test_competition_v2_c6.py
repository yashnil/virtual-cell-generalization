"""C6 uncertainty module: moments, delta-method SE, tau², EB shrinkage, fusion identity."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from scipy import sparse

from virtual_cell.competition_v2 import fusion, fusion_c3, sources
from virtual_cell.competition_v2 import uncertainty as U

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "outputs" / "competition_v2" / "c1_license_clean" / "sources"


def test_moments_match_direct_computation():
    rng = np.random.default_rng(0)
    x = rng.poisson(rng.gamma(1.0, 3.0, size=(300, 12))).astype(float) + (np.arange(12) == 0)
    labels = np.array(["a"] * 120 + ["b"] * 100 + ["ctrl"] * 80)
    groups = ["a", "b", "ctrl"]
    ids = np.array([groups.index(v) for v in labels])
    half = U.half_assignment(labels, ["a", "b"], seed=3)
    m = U.Moments(groups, 12)
    for lo in range(0, 300, 37):
        m.update(sparse.csr_matrix(x[lo : lo + 37]), ids[lo : lo + 37], half[lo : lo + 37])
    r = m.result()
    cpm = x / x.sum(1, keepdims=True) * 1e6
    for gi, g in enumerate(groups):
        sel = labels == g
        np.testing.assert_allclose(r["mean"][gi], cpm[sel].mean(0), rtol=1e-10)
        np.testing.assert_allclose(r["var"][gi], cpm[sel].var(0, ddof=1), rtol=1e-8)
    for gi, g in enumerate(["a", "b"]):
        for h in (0, 1):
            sel = (labels == g) & (half == h)
            assert r["half_n"][gi, h] == sel.sum()
            np.testing.assert_allclose(r["half_mean"][gi, h], cpm[sel].mean(0), rtol=1e-10)
    assert abs(r["half_n"][0, 0] - r["half_n"][0, 1]) <= 1


def test_delta_sigma2_matches_monte_carlo():
    rng = np.random.default_rng(1)
    n_t, n_c, mu_t, mu_c, sd = 200, 4000, 80.0, 50.0, 40.0
    f = 0.7
    sims = []
    for _ in range(4000):
        mt = rng.normal(mu_t, sd, n_t).mean()
        mc = rng.normal(mu_c, sd, n_c).mean()
        m = f * mt + (1 - f) * mc
        sims.append(np.log2((m + 1) / (mc + 1)))
    m0 = f * mu_t + (1 - f) * mu_c
    got = U.delta_sigma2(m0, mu_c, f, sd**2, n_t, sd**2, n_c)
    assert got == pytest.approx(np.var(sims), rel=0.08)


def test_pooled_tau2_recovers_heterogeneity_and_zero():
    rng = np.random.default_rng(2)
    P, G, tau2, s2 = 4000, 3, np.array([0.0, 0.25, 1.0]), 0.5
    truth = rng.normal(0, 1, (P, G))
    eff, sig, msk = {}, {}, {}
    for n in ("a", "b"):
        eff[n] = truth + rng.normal(0, np.sqrt(tau2), (P, G)) + rng.normal(0, np.sqrt(s2), (P, G))
        sig[n] = np.full((P, G), s2)
        msk[n] = np.ones((P, G), dtype=bool)
    est = U.pooled_tau2(eff, sig, msk)
    # between-source variance of two draws is 2 * tau2_source; DL estimates it per source
    np.testing.assert_allclose(est, tau2, atol=0.08)


def test_eb_multiplier_shrinks_noisy_cells_and_keeps_missing():
    e = np.array([[1.0, 0.1], [-1.0, -0.1], [0.5, 0.05]])
    s = np.array([[0.1, 1.0], [0.1, 1.0], [np.nan, 1.0]])
    m, tau2 = U.eb_multiplier(e, s, np.array([True, True, False]))
    assert tau2[0] == pytest.approx(1.0 - 0.1) and tau2[1] == 0.0
    assert m[0, 0] == pytest.approx(0.9 / 1.0)
    assert m[2, 0] == 1.0  # non-finite sigma: C1 treatment
    assert (m[:, 1] == 0).all()  # no signal above noise: shrunk to zero


def test_cd4_sigma2_is_mean_variance_over_usable_conditions():
    st = {
        "available": np.array([[True], [True], [False]]),
        "quality_pass": np.array([[True], [True], [True]]),
        "n_cells": np.array([[50], [50], [50]]),
        "lfcSE": np.array([[[0.2, 0.4]], [[0.4, np.nan]], [[9.0, 9.0]]]),
        "log2fc": np.zeros((3, 1, 2)),
    }
    got = U.cd4_sigma2(st)
    assert got[0, 0] == pytest.approx((0.04 + 0.16) / 4)
    assert got[0, 1] == pytest.approx(0.16 / 1)


@pytest.mark.skipif(
    not (SRC / "K562_GWPS_CPM_full_statistics.npz").exists(), reason="frozen C1 statistics absent"
)
def test_unit_weight_fusion_equals_c1_equal_fusion():
    k562 = sources.load_source(SRC / "K562_GWPS_CPM_full_statistics.npz")
    h1 = sources.load_source(SRC / "H1_2025_full_statistics.npz")
    with np.load(SRC / "CD4_DE_statistics.npz") as d:
        cd4 = {k: d[k] for k in d.files}
    genes = np.intersect1d(k562.genes, h1.genes)[:500]
    targets = np.array(sorted(set(k562.targets[:30]) | set(h1.targets[:30])))
    rng = np.random.default_rng(0)
    controls = {s: rng.dirichlet(np.ones(len(genes))) for s in ("log2fc", "bulk_delta")}
    want = fusion.fused_effects([k562, h1], [1.0, 1.0], cd4, 1.0, targets, genes, controls)
    comp = fusion_c3.components({"K562": k562, "H1": h1}, cd4, targets, genes)
    ones = {"K562": 1.0, "H1": 1.0, "CD4": 1.0}
    for space in want:
        got = U.fuse(comp, space, controls[space], ones)
        np.testing.assert_allclose(got, want[space], rtol=1e-5, atol=1e-6)
