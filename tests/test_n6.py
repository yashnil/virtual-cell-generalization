"""N6 outcome-rule checks on synthetic bootstrap tables (reports/n6_protocol.md §5)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "research_v3"))
import analyse_n6 as a  # noqa: E402

GATE_OK = {
    "viperturb_pooled_split_half_reliability": 0.3,
    "viperturb_reliable_energy_boot_lo": 5.0,
    "n_panel": 600,
}


def _boot(vip, gwps=0.82, jurkat=0.56, rpe1=0.38, sd=0.01, ctrl=None, n=2000, seed=0):
    rng = np.random.default_rng(seed)
    cols = {
        "VIPerturb_K562": vip,
        "K562_GWPS": gwps,
        "K562_GWPS_vipdepth": ctrl or gwps,
        "RPE1": rpe1,
        "Jurkat": jurkat,
    }
    b = pd.DataFrame({k: v + sd * rng.normal(size=n) for k, v in cols.items()})
    return {k: float(b[k].mean()) for k in b}, b


def test_outcome_a_when_vip_approaches_gwps():
    C, b = _boot(0.80)
    assert a.classify(C, b, GATE_OK, adj2_vip_lo=0.05)["outcome"] == "A"


def test_outcome_b_when_vip_resembles_non_k562():
    C, b = _boot(0.57)
    assert a.classify(C, b, GATE_OK, adj2_vip_lo=0.05)["outcome"] == "B"


def test_outcome_c_when_intermediate():
    C, b = _boot(0.69)
    assert a.classify(C, b, GATE_OK, adj2_vip_lo=0.05)["outcome"] == "C"


def test_outcome_d_on_any_gate_failure():
    C, b = _boot(0.80)
    low = dict(GATE_OK, viperturb_pooled_split_half_reliability=0.05)
    assert a.classify(C, b, low, adj2_vip_lo=0.05)["outcome"] == "D"
    C, b = _boot(0.80, ctrl=0.65)  # positive control shifts at VIPerturb depth
    assert a.classify(C, b, GATE_OK, adj2_vip_lo=0.05)["outcome"] == "D"
    small = dict(GATE_OK, n_panel=200)
    C, b = _boot(0.80)
    assert a.classify(C, b, small, adj2_vip_lo=0.05)["outcome"] == "D"


def test_not_robust_to_adjustment_is_inconclusive():
    C, b = _boot(0.80)
    assert a.classify(C, b, GATE_OK, adj2_vip_lo=-0.01)["outcome"].startswith("inconclusive")


def test_build_viperturb_matches_manual_delta(tmp_path, monkeypatch):
    import anndata as ad
    import build_n6 as b
    from scipy import sparse

    rng = np.random.default_rng(3)
    labels = ["NO-TARGET"] * 80 + ["WRB"] * 40 + ["PRMT7"] * 36
    n, G = len(labels), 300
    dense0 = rng.poisson(2.0, size=(n, G)).astype(np.int32)
    lab0 = np.array(labels)
    dense0[lab0 == "WRB", G - 2] = 0  # on-target knockdown
    dense0[lab0 == "PRMT7", G - 1] = 0
    X = sparse.csr_matrix(dense0)
    var = [f"g{i}" for i in range(G - 2)] + ["WRB", "PRMT7"]
    obs = pd.DataFrame(
        {"gene": pd.Categorical(labels), "guide": pd.Categorical([x + "g0" for x in labels])},
        index=[f"c{i}" for i in range(n)],
    )
    path = tmp_path / "vip.h5ad"
    ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index(var, name="gene"))).write_h5ad(path)
    monkeypatch.setattr(b, "VIP", path)
    genes_hgnc = ["g1", "g5", "GET1", "PRMT7"]
    vip_genes = ["g1", "g5", "WRB", "PRMT7"]
    fm, pm, n_full, n_parts, g = b.build_viperturb(
        ["GET1", "PRMT7"], ["WRB", "PRMT7"], genes_hgnc, vip_genes
    )
    dense = X.toarray().astype(float)
    val = np.log1p(1e4 * dense / dense.sum(axis=1, keepdims=True))
    cols = [var.index(v) for v in vip_genes]
    lab = np.array(labels)
    manual = val[lab == "WRB"][:, cols].mean(0) - val[lab == "NO-TARGET"][:, cols].mean(0)
    np.testing.assert_allclose(fm[0] - fm[2], manual, rtol=1e-10)
    assert n_full.tolist() == [40, 36, 80]
