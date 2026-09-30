"""Final-panel readiness: the panel-agnostic C1 pipeline on a small mock D/E/F bundle.

Data-gated: the mock bundle is subsampled from the real validation controls, and the
sources are the frozen C1 statistics. Covered:

* panel loading (D/E/F, arbitrary count, shuffled order, strict gene axis);
* the registry license gate (blocked sources are never fused);
* coverage (unsupported / single / multi-source targets);
* the source-preparation policy;
* end-to-end emission invariants and determinism;
* the frozen unsupported-target fallback;
* a regression against the committed C1 builder's output.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest
import yaml

from virtual_cell.competition_v2 import final as fin
from virtual_cell.competition_v2 import fusion, generator, licensing, panel, sources

ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT / "data" / "raw" / "arc2026" / "controls"
REG = ROOT / "configs" / "source_registry.yaml"
SUPPORT = ROOT / "data" / "splits" / "arc_target_support_license_clean_v1.csv"
FROZEN = ROOT / "outputs" / "competition_v2" / "c1_license_clean" / "sources"
NEEDED = [
    VAL / "context_A.h5ad",
    FROZEN / "K562_GWPS_CPM_full_statistics.npz",
    ROOT / "data/raw/competition_v2/gencode.v47.annotation.gtf.gz",
]
pytestmark = pytest.mark.skipif(
    not all(p.exists() for p in NEEDED), reason="validation controls / frozen C1 statistics absent"
)
CELLS = 5
N_CTRL = 120


def _pick_targets() -> list[str]:
    sup = pd.read_csv(SUPPORT).set_index("arc_target")
    n = sup.n_usable_contexts_c1.clip(upper=3)
    pick = []
    for k in (3, 2, 1, 0):
        pick += list(n.index[n == k][:2])
    return pick


@pytest.fixture(scope="module")
def mock(tmp_path_factory):
    d = tmp_path_factory.mktemp("mock_def")
    targets = list(np.random.default_rng(0).permutation(_pick_targets()))
    for old, new in zip("ABC", "DEF", strict=True):
        a = ad.read_h5ad(VAL / f"context_{old}.h5ad")[:N_CTRL].copy()
        a.obs["context"] = pd.Categorical([new] * a.n_obs)
        a.obs["target_gene"] = pd.Categorical(a.obs["target_gene"].astype(str))
        a.obs = a.obs[["target_gene", "context"]]
        a.obs_names = pd.Index([f"{new}{i}" for i in range(a.n_obs)], dtype=object)
        a.var_names = pd.Index(a.var_names.astype(str), dtype=object)
        a.write_h5ad(d / f"context_{new}.h5ad")
    shutil.copyfile(VAL / "gene_names.csv", d / "gene_names.csv")
    pd.DataFrame({"target_gene": targets}).to_csv(d / "pert_counts.csv", index=False)
    (d / "manifest.json").write_text(
        json.dumps(
            {
                "partition": "MOCK-TEST",
                "panel_id": "pytest-DEF",
                "contexts": ["D", "E", "F"],
                "pert_col": "target_gene",
                "context_col": "context",
                "control_label": "non-targeting",
                "n_genes": 18533,
                "n_constructs": len(targets),
                "cells_per_pert": CELLS,
            }
        )
    )
    reg = yaml.safe_load(REG.read_text())
    for s in reg["sources"]:  # mechanics tests: reuse frozen statistics as given
        if s["role"] == "c1_source":
            s["panel_dependent"] = False
    (d / "registry.yaml").write_text(yaml.safe_dump(reg))
    return d, targets


@pytest.fixture(scope="module")
def emitted(mock, tmp_path_factory):
    d, _ = mock
    p = panel.load_panel(d)
    reg = panel.load_registry(d / "registry.yaml")
    out = tmp_path_factory.mktemp("emit")
    prepared = fin.prepare_sources(p, reg, out, log=lambda m: None)
    fin.emit(p, reg, prepared, out / "prediction.h5ad", jobs=2, log=lambda m: None)
    return p, reg, prepared, out


# ----------------------------------------------------------------------------- panel


def test_panel_loads_def_contexts_count_and_order(mock):
    d, targets = mock
    p = panel.load_panel(d)
    assert p.contexts == ("D", "E", "F")
    assert list(p.targets) == targets and len(targets) == 8
    assert p.cells_per_pert == CELLS and p.n_cells == 3 * 8 * CELLS
    assert np.array_equal(p.genes, pd.read_csv(VAL / "gene_names.csv").gene_name.astype(str))
    assert set(p.checksums) >= {"gene_names.csv", "pert_counts.csv", "context_D.h5ad"}
    assert p.summary()["target_list_sha256"] == panel.list_sha(targets)


def test_panel_rejects_gene_order_and_inconsistencies(mock, tmp_path):
    d, _ = mock
    bad = tmp_path / "bad"
    shutil.copytree(d, bad)
    g = pd.read_csv(bad / "gene_names.csv")
    g.iloc[[0, 1]] = g.iloc[[1, 0]].to_numpy()
    g.to_csv(bad / "gene_names.csv", index=False)
    with pytest.raises(panel.PanelError, match="gene axis"):
        panel.load_panel(bad, checksums=False)
    with pytest.raises(panel.PanelError, match="disagrees"):
        panel.load_panel(d, cells_per_pert=CELLS + 1, checksums=False)
    (bad / "context_E.h5ad").unlink()
    with pytest.raises(panel.PanelError, match="cannot identify control files"):
        panel.load_panel(bad, checksums=False)


# ----------------------------------------------------------------------------- registry


def test_registry_never_fuses_blocked_or_audit_sources(tmp_path):
    reg = panel.load_registry(REG)
    names = [e.name for e in reg.c1_sources()]
    assert names == ["K562", "H1", "CD4"]
    assert "HCT116" not in names and "KOLF2.1J_iPSC" not in names
    raw = yaml.safe_load(REG.read_text())
    for s in raw["sources"]:
        if s["name"] == "HCT116":
            s["role"], s["qualification"], s["c1_order"], s["kind"] = (
                "c1_source", "QUALIFIED", 4, "cell_statistics",
            )  # fmt: skip
    (tmp_path / "r.yaml").write_text(yaml.safe_dump(raw))
    with pytest.raises(licensing.LicenseError):
        panel.load_registry(tmp_path / "r.yaml").c1_sources()
    for s in raw["sources"]:
        if s["name"] == "HCT116":
            s["license_status"] = "GREEN"  # a registry cannot override licensing.STATUS
    (tmp_path / "r2.yaml").write_text(yaml.safe_dump(raw))
    with pytest.raises(licensing.LicenseError):
        panel.load_registry(tmp_path / "r2.yaml")


def test_coverage_counts_unsupported_single_and_multi(emitted, mock):
    p, reg, prepared, _ = emitted
    _, _, _, _, stats = fin.load_c1_inputs(reg, prepared)
    table, summary = panel.coverage_audit(p.targets, reg, stats)
    sup = pd.read_csv(SUPPORT).set_index("arc_target").n_usable_contexts_c1
    assert (table.n_c1_sources == sup.loc[p.targets].to_numpy()).all()
    assert summary["c1_source_count_distribution"] == {"0": 2, "1": 2, "2": 2, "3+": 2}
    assert set(summary["unsupported_targets"]) == set(sup.index[sup == 0]) & set(p.targets)
    assert "HCT116" in summary["audit_sources"]
    assert summary["audit_sources"]["HCT116"]["qualification"] == "BLOCKED"


def test_prepare_policy_reuses_for_validation_and_rebuilds_for_new_panel(
    mock, tmp_path, monkeypatch
):
    reg = panel.load_registry(REG)
    val = panel.load_panel(VAL, checksums=False)
    got = fin.prepare_sources(val, reg, tmp_path / "v", log=lambda m: None)
    assert all(v["action"].startswith("reused") for v in got.values())
    calls = {}

    def subset(path, retained):
        with np.load(path) as z:
            s = {k: z[k] for k in z.files}
        rows = pd.Index(s["targets"].astype(str)).get_indexer(retained)
        assert (rows >= 0).all()
        for k, v in s.items():
            if isinstance(v, np.ndarray) and v.ndim >= 1 and len(v) == len(s["targets"]):
                s[k] = v[rows]
        return s

    def fake_k562(raw, retained, log=print):
        calls["k562"] = list(retained)
        return subset(FROZEN / "K562_GWPS_CPM_full_statistics.npz", retained), None

    def fake_cd4(raw, retained):
        calls["cd4"] = list(retained)
        with np.load(FROZEN / "CD4_DE_statistics.npz") as z:
            s = {k: z[k] for k in z.files}
        idx = pd.Index(s["targets"].astype(str)).get_indexer(retained)
        per_target = ("log2fc", "adjusted_p", "lfcSE", "available", "quality_pass", "n_cells",
                      "rows_per_result")  # fmt: skip
        return {**s, "targets": s["targets"][idx], **{k: s[k][:, idx] for k in per_target}}

    monkeypatch.setattr(sources, "prepare_k562", fake_k562)
    monkeypatch.setattr(sources, "prepare_cd4", fake_cd4)
    d, targets = mock
    new = panel.load_panel(d, checksums=False)
    got = fin.prepare_sources(new, reg, tmp_path / "n", log=lambda m: None)
    expected = fin.retained_for_panel(new, reg)
    assert expected[: len(targets)] == targets
    assert calls["k562"] == expected == calls["cd4"]
    assert got["K562"]["action"].startswith("prepared") and got["CD4"]["action"].startswith(
        "prepared"
    )
    assert got["H1"]["action"].startswith("reused (panel-independent)")


# ----------------------------------------------------------------------------- emission


def _obs(path, col):
    with h5py.File(path, "r") as f:
        node = f["obs"][col]
        cats = np.array(
            [c.decode() if isinstance(c, bytes) else str(c) for c in node["categories"][()]]
        )
        return cats[node["codes"][()]]


def test_emission_invariants(emitted):
    p, _, _, out = emitted
    pred = out / "prediction.h5ad"
    val = fin.validate(p, pred)
    assert val["all_pass"], {k: v for k, v in val["checks"].items() if not v}
    ctx = _obs(pred, "context")
    tgt = _obs(pred, "target_gene")
    assert list(dict.fromkeys(ctx)) == ["D", "E", "F"]
    assert np.array_equal(tgt[: 8 * CELLS : CELLS], p.targets)  # shuffled order kept
    assert "non-targeting" not in set(tgt)
    counts = pd.Series(list(zip(ctx, tgt, strict=True))).value_counts()
    assert (counts == CELLS).all() and len(counts) == 24
    with h5py.File(pred, "r") as f:
        data = f["X/data"][()]
    assert (data > 0).all() and (data == np.floor(data)).all()


def test_emission_is_deterministic(emitted, tmp_path):
    p, reg, prepared, out = emitted
    fin.emit(p, reg, prepared, tmp_path / "again.h5ad", jobs=1, log=lambda m: None)
    with h5py.File(out / "prediction.h5ad", "r") as a, h5py.File(tmp_path / "again.h5ad", "r") as b:
        for k in ("data", "indices", "indptr"):
            assert np.array_equal(a["X"][k][()], b["X"][k][()])


def test_unsupported_targets_use_the_frozen_zero_effect_fallback(emitted):
    p, reg, prepared, _ = emitted
    cell, w, cd4, cd4_w, stats = fin.load_c1_inputs(reg, prepared)
    pairs = generator.promoter_pairs(
        generator.gencode_tss(ROOT / "data/raw/competition_v2/gencode.v47.annotation.gtf.gz"),
        p.targets,
        p.genes,
    )
    _, _, p_cpm, p_bulk, _, effects = fin.context_moments(p, 0, "D", cell, w, cd4, cd4_w, pairs)
    _, summary = panel.coverage_audit(p.targets, reg, stats)
    rows = [list(p.targets).index(t) for t in summary["unsupported_targets"]]
    assert rows and all(not np.any(effects[s][rows]) for s in effects)
    ctrl = ad.read_h5ad(p.control_files["D"])
    _, _, mean, bulk = generator.control_template(ctrl.X, CELLS, generator.POOL_K, generator.SEED)
    zero = {s: np.zeros_like(v) for s, v in effects.items()}
    z_cpm, z_bulk = generator.expected_moments(
        zero, {"log2fc": mean, "bulk_delta": bulk}, p.targets, p.genes, pairs=pairs
    )
    np.testing.assert_array_equal(p_cpm[rows], z_cpm[rows])
    np.testing.assert_array_equal(p_bulk[rows], z_bulk[rows])


def test_single_source_target_keeps_its_source_response(emitted):
    p, reg, prepared, _ = emitted
    cell, w, cd4, cd4_w, stats = fin.load_c1_inputs(reg, prepared)
    table, _ = panel.coverage_audit(p.targets, reg, stats)
    one = table.index[table.n_c1_sources == 1][0]
    src = [c.split(":")[1] for c in table.columns if c.startswith("c1:") and table.loc[one, c]][0]
    controls = {
        "log2fc": np.full(len(p.genes), 1 / len(p.genes)),
        "bulk_delta": np.full(len(p.genes), 1 / len(p.genes)),
    }
    full = fusion.fused_effects(cell, w, cd4, cd4_w, [one], p.genes, controls)
    solo = (
        fusion.fused_effects([stats[src]], [1.0], None, 0.0, [one], p.genes, controls)
        if src != "CD4"
        else fusion.fused_effects([], [], cd4, 1.0, [one], p.genes, controls)
    )
    for s in full:
        np.testing.assert_allclose(full[s], solo[s], rtol=1e-6, atol=1e-7)


# ----------------------------------------------------------------------------- regression


COMMITTED = ROOT / "outputs" / "final" / "c1_committed_builder_abc" / "prediction.h5ad"


@pytest.mark.skipif(not COMMITTED.exists(), reason="committed-builder regression output absent")
def test_generic_blocks_equal_the_committed_c1_builder_on_abc():
    val = panel.load_panel(VAL, checksums=False)
    reg = panel.load_registry(REG)
    prepared = fin.prepare_sources(
        val, reg, ROOT / "outputs" / "final" / "_unused", log=lambda m: None
    )
    assert all(v["action"].startswith("reused") for v in prepared.values())
    cell, w, cd4, cd4_w, _ = fin.load_c1_inputs(reg, prepared)
    pairs = generator.promoter_pairs(
        generator.gencode_tss(ROOT / "data/raw/competition_v2/gencode.v47.annotation.gtf.gz"),
        val.targets,
        val.genes,
    )
    ci, ctx = 2, "C"
    template, depths, p_cpm, p_bulk, _, _ = fin.context_moments(
        val, ci, ctx, cell, w, cd4, cd4_w, pairs
    )
    with h5py.File(COMMITTED, "r") as f:
        ip = f["X/indptr"]
        for ti in (0, 1, 150, 299):
            t = val.targets[ti]
            got = generator.dual_moment_counts(
                template,
                p_cpm[ti],
                p_bulk[ti],
                depths=depths,
                seed=generator.seed_for(f"{ctx}:{t}"),
            )
            lo = (ci * len(val.targets) + ti) * 400
            a, b = int(ip[lo]), int(ip[lo + 400])
            want = np.zeros((400, len(val.genes)), dtype=np.float32)
            rows = np.repeat(np.arange(400), np.diff(ip[lo : lo + 401]))
            want[rows, f["X/indices"][a:b]] = f["X/data"][a:b]
            np.testing.assert_array_equal(got.astype(np.float32), want)
