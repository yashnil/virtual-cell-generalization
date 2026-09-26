"""Figure-source reproducibility and the invariants this phase must not break.

Covers: deterministic extraction from frozen artifacts, finite source values,
exact Arc tier counts, provenance records, frozen-manifest integrity, the frozen
Arc model parameters, and that the dry-run bundle was not regenerated.
Data-gated tests skip when the git-ignored artifacts are absent.
"""

from __future__ import annotations

import io
import json
import re
from pathlib import Path

import pandas as pd
import pytest

from virtual_cell.visualization import common, sources

ROOT = Path(__file__).resolve().parents[1]
FROZEN_EXTRACTORS = [n for n in sources.EXTRACTORS if n != "fig7_source_reliability"]


def _have(ex_name: str) -> bool:
    try:
        ex = sources.EXTRACTORS[ex_name]
        return all((ROOT / s).exists() for s in ex().sources)
    except FileNotFoundError:
        return False


def _csv_text(table: pd.DataFrame) -> str:
    buf = io.StringIO()
    table.to_csv(buf, index=False, float_format="%.10g")
    return buf.getvalue()


@pytest.mark.parametrize("name", list(sources.EXTRACTORS))
def test_figure_source_is_reproducible_from_artifacts(name):
    if not _have(name):
        pytest.skip("source artifacts absent")
    path = common.SOURCES_DIR / f"{name}.csv"
    if not path.exists():
        pytest.skip("figure source not generated yet")
    assert _csv_text(sources.EXTRACTORS[name]().table) == path.read_text()


@pytest.mark.parametrize("name", list(sources.EXTRACTORS))
def test_figure_source_has_no_unintended_nan_or_inf(name):
    path = common.SOURCES_DIR / f"{name}.csv"
    if not path.exists():
        pytest.skip("figure source not generated yet")
    allow = sources.EXTRACTORS[name]().allow_nan if _have(name) else []
    assert common.numeric_is_finite(pd.read_csv(path), allow)


@pytest.mark.parametrize("name", list(sources.EXTRACTORS))
def test_provenance_records_current_source_hashes(name):
    prov = common.SOURCES_DIR / f"{name}.provenance.json"
    if not prov.exists():
        pytest.skip("figure source not generated yet")
    record = json.loads(prov.read_text())
    for key in (
        "source_files",
        "report",
        "extraction_script",
        "figure_script",
        "generated",
        "git_head",
    ):
        assert key in record
    assert (ROOT / record["figure_script"]).exists()
    assert (ROOT / record["report"]).exists()
    for src in record["source_files"]:
        p = ROOT / src["path"]
        if p.exists():
            assert common.sha256(p) == src["sha256"], src["path"]


def test_frozen_figure_sources_come_from_frozen_artifacts():
    for name in FROZEN_EXTRACTORS:
        prov = common.SOURCES_DIR / f"{name}.provenance.json"
        if not prov.exists():
            continue
        for src in json.loads(prov.read_text())["source_files"]:
            assert src["freezes"], f"{name}: {src['path']} is not pinned by any freeze"


def test_numeric_is_finite_flags_inf_and_unexpected_nan():
    t = pd.DataFrame({"a": [1.0, float("nan")], "b": [1.0, 2.0]})
    assert not common.numeric_is_finite(t)
    assert common.numeric_is_finite(t, allow_nan=["a"])
    assert not common.numeric_is_finite(pd.DataFrame({"a": [float("inf")]}), allow_nan=["a"])


# ----------------------------------------------------------------- Arc tiers


def test_arc_tier_counts_are_exact():
    s = pd.read_csv(ROOT / "data/splits/arc_target_support_v1.csv")
    assert len(s) == 300 and s.arc_target.is_unique
    assert s.support_tier.value_counts().to_dict() == {0: 214, 1: 79, 2: 7}


def test_fig6_source_sums_to_300_with_frozen_tiers():
    t = sources.arc_target_support().table
    assert int(t.n_targets.sum()) == 300
    assert t.groupby("tier").n_targets.sum().to_dict() == {0: 214, 1: 79, 2: 7}
    assert t.set_index("provenance").n_targets.to_dict() == {
        "both arch1 and Kaden": 7,
        "arch1 only": 6,
        "Kaden only": 73,
        "no public perturbation": 214,
    }


# --------------------------------------------------- frozen model invariants


def test_frozen_arc_model_constants_unchanged():
    src = (ROOT / "scripts/run_arc_dry_run.py").read_text()
    assert 'MAIN_EFFECT = "M3b_basal_shrunk"' in src
    assert "TIER_WEIGHTS = {2: 0.50, 1: 0.25, 0: 0.0}" in src
    assert "SEED = 20260921" in src
    assert "CELLS_PER_PERT = 400" in src
    assert re.search(r"smoothing=0\.5\b", src)
    assert 'BETA_SOURCES = ("arch1", "kaden25rpe1")' in src


def test_frozen_dry_run_summary_values_unchanged():
    path = ROOT / "outputs/arc_dry_run_v1/summary.json"
    if not path.exists():
        pytest.skip("dry-run outputs absent")
    s = json.loads(path.read_text())
    assert s["tier_weights"] == {"2": 0.5, "1": 0.25, "0": 0.0}
    assert s["main_effect_estimator"] == "M3b_basal_shrunk"
    assert s["fitted_main_effect_scale"] == pytest.approx(0.11890236823138407, abs=0)
    assert s["submitted"] is False


def test_dry_run_bundle_was_not_regenerated():
    bundle = ROOT / "outputs/arc_dry_run_v1/arc_dry_run_v1.h5ad"
    summary = ROOT / "outputs/arc_dry_run_v1/summary.json"
    if not bundle.exists():
        pytest.skip("dry-run bundle absent")
    s = json.loads(summary.read_text())
    assert bundle.stat().st_size == s["bundle_bytes"] == 16_489_293_744
    # the bundle is written before summary.json, which the arc_count_space freeze pins
    assert bundle.stat().st_mtime <= summary.stat().st_mtime + 1


def test_every_freeze_manifest_still_verifies():
    manifests = sorted((ROOT / "data/provenance/scperteval").glob("*_freeze.txt"))
    assert len(manifests) == 15
    checked = 0
    for m in manifests:
        for line in m.read_text().splitlines():
            parts = line.split()
            if len(parts) != 2 or len(parts[0]) != 64:
                continue
            p = ROOT / parts[1]
            if not p.exists():
                assert parts[1].startswith(("data/raw", "data/processed", "outputs/")), parts[1]
                continue
            assert common.sha256(p) == parts[0], f"{m.name}: {parts[1]} changed"
            checked += 1
    assert checked > 0
