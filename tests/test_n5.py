"""N5 checks on synthetic data only (reports/n5_protocol.md §6)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

from virtual_cell.analysis import source_compatibility as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "research_v3"))


def _fields(rng, P=300, G=400, cos=0.6):
    base = rng.normal(size=(P, G))
    other = rng.normal(size=(P, G))
    L_T = base
    L_S = cos * base + np.sqrt(1 - cos**2) * other
    return L_S, L_T


def _noisy(L, rng, sd, R=5):
    a = [L + sd * rng.normal(size=L.shape) for _ in range(R)]
    b = [L + sd * rng.normal(size=L.shape) for _ in range(R)]
    full = L + sd / np.sqrt(2) * rng.normal(size=L.shape)
    return full, a, b


def test_latent_cosine_recovers_planted_value_under_unequal_noise():
    rng = np.random.default_rng(0)
    L_S, L_T = _fields(rng, cos=0.6)
    S_full, S_a, S_b = _noisy(L_S, rng, sd=3.0)  # very noisy source
    T_full, T_a, T_b = _noisy(L_T, rng, sd=1.0)
    t = sc.latent_cosine_terms(S_full, T_full, S_a, S_b, T_a, T_b)
    C = sc.pooled_latent_cosine(t["num"], t["den_S"], t["den_T"])
    truth = float(
        np.sum(sc.centre(L_S) * sc.centre(L_T))
        / np.sqrt(np.sum(sc.centre(L_S) ** 2) * np.sum(sc.centre(L_T) ** 2))
    )
    assert abs(C - truth) < 0.03
    raw = sc.pooled_raw_cosine(S_full, T_full)
    raw_c = raw["num"].sum() / np.sqrt(raw["ss"].sum() * raw["tt"].sum())
    assert raw_c < truth - 0.2  # raw cosine is attenuated by source noise


def test_fixed_effect_regression_recovers_contrast():
    rng = np.random.default_rng(1)
    P, sources = 500, ["A", "B", "C"]
    pert = np.tile(np.arange(P), 3)
    src = np.repeat(np.array(sources), P)
    rel = rng.uniform(0, 1, size=3 * P)
    alpha = rng.normal(size=P)[pert]
    beta = {"A": 0.2, "B": 0.0, "C": -0.1}
    y = alpha + np.array([beta[s] for s in src]) + 0.5 * rel + 0.05 * rng.normal(size=3 * P)
    est = sc.fe_regression(y, pert, src, np.column_stack([rel, rel**2]), sources, "B")
    assert abs(est["A"] - 0.2) < 0.02 and abs(est["C"] + 0.1) < 0.02


def test_matched_comparator_is_source_only_and_deterministic():
    rel = {"GWPS": np.full(10, 0.30), "X": np.full(10, 0.31), "Y": np.full(10, 0.9)}
    cells = {"GWPS": np.full(10, 200.0), "X": np.full(10, 50.0), "Y": np.full(10, 200.0)}
    assert sc.matched_comparator(rel, cells, "GWPS")["matched"] == "X"


def test_part_assignment_sizes_and_disjointness():
    import build_n5_gwps as b

    groups = np.repeat(np.arange(4), [30, 31, 57, 100])
    rows = np.arange(len(groups))
    out = b.assign_parts(groups, 4, (1, 2), rows)
    for r in range(out.shape[0]):
        codes = out[r][out[r] >= 0]
        part, grp = codes // 4, codes % 4
        assert np.array_equal(grp, groups[out[r] >= 0])
        for g, n in enumerate([30, 31, 57, 100]):
            counts = np.bincount(part[grp == g], minlength=3)
            assert counts.tolist() == [n // 2, n // 4, n // 4]
