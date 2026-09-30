"""C1 reimplementation vs the vendored AtlasShift (C0) code: numerical equality.

The vendored files are loaded here, in the test-suite only, so our package never imports
them. Synthetic inputs exercise every function; the real-data tests (skipped when the
prepared statistics are absent) repeat the check on a deterministic subset of the C0
inputs, including X-Atlas, which is *read* here only to prove the implementation.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from virtual_cell.competition_v2 import atlas, fusion, generator, sources

ROOT = Path(__file__).resolve().parents[1]
VENDORED = ROOT / "third_party" / "atlasshift"
C0_DATA = ROOT / "outputs" / "competition_v2" / "atlasshift_c0" / "data"


def _load(name: str, file: str):
    spec = importlib.util.spec_from_file_location(name, VENDORED / file)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def up():
    saved = sys.modules.get("model")
    model = _load("model", "model.py")  # predict.py does ``from model import ...``
    predict = _load("atlasshift_predict", "predict.py")
    yield model, predict
    if saved is None:
        sys.modules.pop("model", None)
    else:
        sys.modules["model"] = saved


def _raw_stats(rng, targets, genes, n_batches_ctrl=True, name="SRC"):
    n_t, n_g = len(targets), len(genes)
    counts = rng.poisson(rng.gamma(0.5, 20, size=(n_t, n_g))).astype(np.float32)
    counts[1] = 0  # a target with no counts
    ctrl = rng.dirichlet(np.ones(n_g), size=n_t).astype(np.float32)
    cpm = (counts / np.maximum(counts.sum(1, keepdims=True), 1) * 1e6).astype(np.float32)
    ccpm = (ctrl * 1e6).astype(np.float32)
    measured = np.ones(n_g, dtype=bool)
    measured[-2:] = False
    n_cells = rng.integers(5, 200, size=n_t)
    n_cells[0] = 10  # below minimum_cells
    return {
        "source": np.asarray(name),
        "targets": np.asarray(targets),
        "genes": np.asarray(genes),
        "n_cells": n_cells,
        "target_count_sums": counts,
        "target_mean_cpm": cpm,
        "matched_control_probability": ctrl,
        "matched_control_mean_cpm": ccpm,
        "global_control_probability": ctrl[0],
        "global_control_mean_cpm": ccpm[0],
        "measured_genes": measured,
    }


@pytest.fixture(scope="module")
def synthetic(tmp_path_factory):
    rng = np.random.default_rng(7)
    genes = np.asarray([f"G{i}" for i in range(40)])
    tmp = tmp_path_factory.mktemp("c1eq")
    paths = []
    for k, name in enumerate(["A", "B"]):
        targets = np.asarray([f"G{i}" for i in range(k, 12 + k)] + ["NOTAGENE"])
        raw = _raw_stats(rng, targets, genes[k:], name=name)
        path = tmp / f"{name}.npz"
        np.savez(path, **raw)
        paths.append(path)
    cd4_targets = np.asarray([f"G{i}" for i in range(3, 15)])
    shape = (3, len(cd4_targets), len(genes) - 5)
    log2fc = rng.normal(0, 0.5, size=shape).astype(np.float32)
    log2fc[1, 2, 4] = np.nan
    cd4 = {
        "targets": cd4_targets,
        "genes": genes[5:],
        "conditions": np.asarray(["Rest", "Stim48hr", "Stim8hr"]),
        "log2fc": log2fc,
        "available": rng.random(shape[:2]) > 0.2,
        "quality_pass": rng.random(shape[:2]) > 0.2,
        "n_cells": rng.integers(10, 60, size=shape[:2]),
    }
    cd4_path = tmp / "cd4.npz"
    np.savez(cd4_path, **cd4)
    panel_targets = np.asarray([f"G{i}" for i in range(0, 16, 2)])
    return {
        "paths": paths,
        "cd4": cd4,
        "cd4_path": cd4_path,
        "genes": genes,
        "targets": panel_targets,
        "tmp": tmp,
    }


def _same(a, b):
    assert a.dtype == b.dtype, (a.dtype, b.dtype)
    np.testing.assert_array_equal(a, b)


def test_source_loading_is_identical(up, synthetic):
    model, _ = up
    for path in synthetic["paths"]:
        u = model.load_xatlas(path, prior_counts=100000)
        o = sources.load_source(path)
        for field in [
            "probability",
            "control_probability",
            "mean_cpm",
            "control_mean_cpm",
            "n_cells",
            "measured",
        ]:
            _same(getattr(o, field), getattr(u, field))


@pytest.mark.parametrize("space", ["log2fc", "bulk_delta"])
@pytest.mark.parametrize("center", [True, False])
def test_centred_source_response_is_identical(up, synthetic, space, center):
    model, _ = up
    for path in synthetic["paths"]:
        u_src = model.load_xatlas(path)
        o_src = sources.load_source(path)
        ue, um = model.aligned_effect(
            u_src, u_src.targets, synthetic["genes"], space=space, common_subtract=int(center)
        )
        oe, om = atlas.source_effect(
            o_src, o_src.targets, synthetic["genes"], space=space, center=center
        )
        _same(oe, ue)
        _same(om, um)


@pytest.mark.parametrize("space", ["log2fc", "bulk_delta"])
def test_fused_effect_and_cd4_family_are_identical(up, synthetic, space):
    model, _ = up
    genes, targets = synthetic["genes"], synthetic["targets"]
    u_src = [model.load_xatlas(p) for p in synthetic["paths"]]
    o_src = [sources.load_source(p) for p in synthetic["paths"]]
    ue, ud = model.fuse_source_centered(
        u_src, [2, 1], targets, genes, space=space, common_subtract=1
    )
    oe, od = fusion.fuse_centered(o_src, [2, 1], targets, genes, space=space)
    _same(oe, ue)
    _same(od, ud)
    uc, ua = model.aligned_cd4(
        synthetic["cd4_path"], targets, genes, common_subtract=1, center_scope="source"
    )
    oc, oa = atlas.cd4_effect(synthetic["cd4"], targets, genes)
    _same(oc, uc)
    _same(oa, ua)
    control = np.random.default_rng(3).dirichlet(np.ones(len(genes)))
    uf, uw = model.add_cd4_family(ue, ud, uc, ua, 0.5, space=space, control_probability=control)
    of, ow = fusion.add_cd4(oe, od, oc, oa, 0.5, space=space, control_probability=control)
    _same(of, uf)
    _same(ow, uw)


@pytest.mark.parametrize("space", ["log2fc", "bulk_delta"])
def test_mapped_fold_change_and_promoter_cap_are_identical(up, synthetic, space):
    model, _ = up
    rng = np.random.default_rng(11)
    genes, targets = synthetic["genes"], synthetic["targets"]
    control = rng.dirichlet(np.ones(len(genes)))
    effect = rng.normal(0, 3, size=(len(targets), len(genes))).astype(np.float32)
    amp = generator.AMPLITUDE[space]
    u = model.desired_mean(control, effect, space=space, amplitude=amp, clip=3)
    o = generator.expected_composition(control, effect, space=space, amplitude=amp)
    _same(o, u)
    pairs = pd.DataFrame(
        {
            "target": ["G0", "G2", "G4", "G6"],
            "neighbor": ["G1", "G3", "G30", "G2"],
            "distance": [10, 800, 4000, 5000],
        }
    )
    path = synthetic["tmp"] / "pairs.csv"
    pairs.to_csv(path, index=False)
    up_p, up_changed = model.apply_promoter_prior(u, control, targets, genes, path, 0.15)
    our_p, our_changed = generator.promoter_cap(o, control, targets, genes, pairs)
    _same(our_p, up_p)
    assert len(our_changed) == len(up_changed) > 0


def test_promoter_pairs_are_identical():
    prepare = _load("atlasshift_prepare", "prepare.py")
    table = pd.DataFrame(
        {
            "gene": ["A", "B", "C", "D", "D", "E"],
            "chromosome": ["1", "1", "1", "2", "2", "1"],
            "strand": ["+", "-", "+", "+", "-", "+"],
            "tss": [1000, 1400, 7000, 50, 60, 5900],
        }
    )
    targets, genes = np.asarray(["A", "D", "E"]), np.asarray(["A", "B", "C", "D", "E"])
    up_pairs = prepare.pairs(table, targets, genes)[["target", "neighbor", "distance"]]
    ours = generator.promoter_pairs(table, targets, genes)
    pd.testing.assert_frame_equal(
        ours.reset_index(drop=True), up_pairs.reset_index(drop=True), check_dtype=False
    )


def test_control_template_and_integer_emission_are_identical(up, tmp_path):
    model, predict = up
    rng = np.random.default_rng(5)
    n_ctrl, n_g, cells = 64, 30, 12
    raw = rng.poisson(rng.gamma(0.8, 3, size=(n_ctrl, n_g))).astype(np.float32)
    raw[:, 0] += 1
    genes = np.asarray([f"G{i}" for i in range(n_g)])
    path = tmp_path / "ctrl.h5ad"
    ad.AnnData(sparse.csr_matrix(raw), var=pd.DataFrame(index=genes)).write_h5ad(path)
    ut, udep, umean, ubulk = predict.control_template(path, genes, cells, 4, 99)
    ot, odep, omean, obulk = generator.control_template(sparse.csr_matrix(raw), cells, 4, 99)
    for a, b in [(ot, ut), (odep, udep), (omean, umean), (obulk, ubulk)]:
        _same(a, b)
    effect = rng.normal(0, 1, size=(1, n_g)).astype(np.float32)
    p = generator.expected_composition(omean, effect, space="log2fc", amplitude=0.6)[0]
    b = generator.expected_composition(obulk, effect, space="bulk_delta", amplitude=0.3)[0]
    seed = generator.seed_for("A:G3")
    assert seed == model.seed_for("A:G3", 20260910)
    uc = model.dual_moment_counts(ut, p, b, depths=udep, seed=seed)
    oc = generator.dual_moment_counts(ot, p, b, depths=odep, seed=seed)
    _same(oc, uc)
    flat = np.full_like(odep, int(odep.mean()))
    _same(
        generator.dual_moment_counts(ot, p, p, depths=flat, seed=1),
        model.dual_moment_counts(ut, p, p, depths=flat, seed=1),
    )
    for fn, tpl in [(generator.dual_moment_counts, ot), (model.dual_moment_counts, ut)]:
        with pytest.raises(ValueError, match="projection too large"):
            fn(tpl, p, b, depths=flat, seed=1)  # equal depths cannot carry a bulk shift


def test_agreement_is_mean_pairwise_cosine_without_panel_targets():
    rng = np.random.default_rng(0)
    n_t, n_g = 4, 6
    a = rng.normal(size=(n_t, n_g)).astype(np.float32)
    b = rng.normal(size=(n_t, n_g)).astype(np.float32)
    c = rng.normal(size=(n_t, n_g)).astype(np.float32)
    ones = np.ones((n_t, n_g), dtype=bool)
    only_one = ones.copy()
    only_one[3] = False
    vectors = {"a": (a, ones), "b": (b, ones), "c": (c, only_one)}
    exclude = np.zeros(n_g, dtype=bool)
    exclude[0] = True
    agree, n_src = fusion.source_agreement(vectors, exclude)

    def cos(x, y):
        return x @ y / np.linalg.norm(x) / np.linalg.norm(y)

    k = ~exclude
    want0 = np.mean([cos(a[0, k], b[0, k]), cos(a[0, k], c[0, k]), cos(b[0, k], c[0, k])])
    assert agree[0] == pytest.approx(want0)
    assert agree[3] == pytest.approx(cos(a[3, k], b[3, k]))
    assert list(n_src) == [3, 3, 3, 2]
    single = {"a": (a, ones), "c": (c, only_one)}
    assert np.isnan(fusion.source_agreement(single, exclude)[0][3])


def test_shrinkage_family():
    agree = np.array([0.1, np.nan, 0.5, 0.3])
    np.testing.assert_array_equal(fusion.shrinkage_factor("S0", 4), np.ones(4))
    np.testing.assert_array_equal(fusion.shrinkage_factor("S1", 4, scalar=0.5), np.full(4, 0.5))
    lam = fusion.shrinkage_factor("S2", 4, floor=0.5, agreement=agree)
    np.testing.assert_allclose(lam, 0.5 + 0.5 * np.array([1 / 6, 0.5, 5 / 6, 0.5]))
    assert np.all(np.diff(lam[[0, 3, 2]]) > 0)  # monotone in agreement
    with pytest.raises(ValueError):
        fusion.shrinkage_factor("S1", 4, scalar=1.5)


# ---------------------------------------------------------------------------
# real C0 inputs, deterministic subset
# ---------------------------------------------------------------------------

NEEDED = [
    "K562_GWPS_CPM_full_statistics.npz",
    "HCT116_full_statistics.npz",
    "HEK293T_full_statistics.npz",
    "H1_2025_full_statistics.npz",
    "CD4_DE_statistics.npz",
    "official_pairs.csv",
    "pert_counts.csv",
    "gene_names.csv",
]


@pytest.mark.skipif(
    not all((C0_DATA / n).exists() for n in NEEDED), reason="C0 prepared statistics not present"
)
def test_real_c0_subset_is_reproduced(up):
    model, _ = up
    targets = pd.read_csv(C0_DATA / "pert_counts.csv").target_gene.astype(str).to_numpy()[::25]
    genes = pd.read_csv(C0_DATA / "gene_names.csv").gene_name.astype(str).to_numpy()
    genes = np.unique(np.concatenate([genes[::20], targets[np.isin(targets, genes)]]))
    names, weights = NEEDED[:4], [2, 1, 1, 2]
    u_src = [model.load_xatlas(C0_DATA / n, prior_counts=100000) for n in names]
    o_src = [sources.load_source(C0_DATA / n) for n in names]
    control = np.random.default_rng(1).dirichlet(np.ones(len(genes)))
    uc, ua = model.aligned_cd4(
        C0_DATA / "CD4_DE_statistics.npz", targets, genes, common_subtract=1, center_scope="source"
    )
    with np.load(C0_DATA / "CD4_DE_statistics.npz") as d:
        oc, oa = atlas.cd4_effect({k: d[k] for k in d.files}, targets, genes)
    _same(oc, uc)
    pairs = pd.read_csv(C0_DATA / "official_pairs.csv")
    for space in ["log2fc", "bulk_delta"]:
        for us, os_ in zip(u_src, o_src, strict=True):  # centred source response
            _same(
                atlas.source_effect(os_, os_.targets, genes, space=space)[0],
                model.aligned_effect(us, us.targets, genes, space=space, common_subtract=1)[0],
            )
        ue, ud = model.fuse_source_centered(
            u_src, weights, targets, genes, space=space, common_subtract=1
        )
        oe, od = fusion.fuse_centered(o_src, weights, targets, genes, space=space)
        _same(oe, ue)  # fused effect
        ue, _ = model.add_cd4_family(ue, ud, uc, ua, 0.5, space=space, control_probability=control)
        oe, _ = fusion.add_cd4(oe, od, oc, oa, 0.5, space=space, control_probability=control)
        _same(oe, ue)
        amp = generator.AMPLITUDE[space]
        um = model.desired_mean(control, ue, space=space, amplitude=amp, clip=3)
        om = generator.expected_composition(control, oe, space=space, amplitude=amp)
        _same(om, um)  # mapped fold change -> composition
        up_p, _ = model.apply_promoter_prior(
            um, control, targets, genes, C0_DATA / "official_pairs.csv", 0.15
        )
        our_p, _ = generator.promoter_cap(om, control, targets, genes, pairs)
        _same(our_p, up_p)  # expected perturbed mean


# ---------------------------------------------------------------------------
# licensing gate and license-clean coverage
# ---------------------------------------------------------------------------


def test_blocked_and_unknown_sources_are_refused():
    from virtual_cell.competition_v2 import licensing

    licensing.assert_sources_allowed(["K562", "H1", "CD4", "GENCODE"])
    for bad in (["HCT116"], ["HEK293T"], ["K562", "HEK293T"], ["SOMETHING_NEW"]):
        with pytest.raises(licensing.LicenseError):
            licensing.assert_sources_allowed(bad)
    assert licensing.STATUS["HCT116"] == licensing.STATUS["HEK293T"] == licensing.BLOCKED


COVERAGE = ROOT / "data" / "splits" / "arc_target_support_license_clean_v1.csv"


@pytest.mark.skipif(not COVERAGE.exists(), reason="license-clean coverage not built")
def test_license_clean_coverage_has_no_xatlas_and_all_targets():
    table = pd.read_csv(COVERAGE)
    assert len(table) == 300 and table.arc_target.is_unique
    assert not any("XAtlas" in c or "HCT116" in c or "HEK293T" in c for c in table.columns)
    usable = table[["K562_GWPS__usable", "H1_2025_full__usable", "CD4_DE__usable"]].sum(axis=1)
    assert (usable == table.n_usable_contexts_c1).all()
    assert not table.Kaden_RPE1__usable.any()
    assert (table.n_usable_contexts_c1 <= table.n_usable_contexts_c0).all()


def test_c1_candidate_code_paths_never_name_blocked_sources():
    """The C1 candidate builder may not reference X-Atlas files at all."""
    text = (ROOT / "scripts" / "competition_v2" / "build_c1_candidate.py").read_text()
    for token in [
        "HCT116_full_statistics",
        "HEK293T_full_statistics",
        "xatlas_orion",
        "atlasshift_c0/data",
    ]:
        assert token not in text
