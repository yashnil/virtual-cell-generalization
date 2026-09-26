"""The Kaden source-reliability diagnostic: predeclaration, leakage and
construction checks against its outputs (skipped when they are absent)."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/run_kaden_source_reliability.py"
MODULE = ROOT / "src/virtual_cell/analysis/source_reliability.py"
OUT = ROOT / "outputs/kaden_source_reliability_v1"
RAW = ROOT / "data/raw/scperteval"

needs_outputs = pytest.mark.skipif(
    not (OUT / "summary.json").exists(), reason="diagnostic outputs absent"
)
needs_raw = pytest.mark.skipif(
    not (RAW / "kaden25rpe1_processed_complete.h5ad").exists(), reason="raw data absent"
)


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def test_script_reads_no_arc_outcome():
    """Only the Arc gene list and target list are read -- never a context file,
    and nothing that could hold a hidden outcome or leaderboard value."""
    tree = ast.parse(SCRIPT.read_text())
    arc_calls = {
        n.func.attr
        for n in ast.walk(tree)
        if isinstance(n, ast.Call)
        and isinstance(n.func, ast.Attribute)
        and isinstance(n.func.value, ast.Name)
        and n.func.value.id == "arc2026"
    }
    assert arc_calls == {"load_gene_names", "load_pert_counts"}
    code = SCRIPT.read_text()
    assert "context_A" not in code and "context_B" not in code and "context_C" not in code
    assert "vcc" not in code.replace("vcc2026", "")


def test_script_does_not_touch_the_frozen_arc_model():
    code = SCRIPT.read_text()
    for forbidden in (
        "mean_response",
        "context_main_effect",
        "arc.generate",
        "arc.bundle",
        "SubmissionWriter",
        "TIER_WEIGHTS",
    ):
        assert forbidden not in code


@needs_outputs
def test_predeclaration_matches_the_code_that_ran():
    pre = json.loads((OUT / "predeclaration.json").read_text())
    assert pre["script_sha256"] == _sha(SCRIPT)
    assert pre["module_sha256"] == _sha(MODULE)
    assert pre["seed"] == 20260925 and pre["n_repeats"] == 50
    # written before any result table
    t0 = (OUT / "predeclaration.json").stat().st_mtime
    for f in ("per_perturbation_reliability.csv", "main_effect_reliability.csv", "summary.json"):
        assert (OUT / f).stat().st_mtime >= t0


@needs_outputs
def test_frozen_main_effect_cosines_are_reproduced():
    s = json.loads((OUT / "summary.json").read_text())
    for rec in s["frozen_reproduction"].values():
        assert abs(rec["recomputed"] - rec["frozen"]) < 1e-9


@needs_outputs
@needs_raw
def test_matched_perturbation_intersections_recompute():
    from virtual_cell.data import scperteval

    def labels(n):
        return set(scperteval.cell_labels(RAW / f"{n}_processed_complete.h5ad").tolist()) - {
            "control"
        }

    a, k, r = labels("arch1"), labels("kaden25rpe1"), labels("replogle22rpe1")
    arc = set(pd.read_csv(ROOT / "data/raw/arc2026/controls/pert_counts.csv").target_gene)
    s = json.loads((OUT / "summary.json").read_text())
    assert s["n_shared_arch1_kaden"] == len(a & k)
    assert s["n_shared_kaden_rpe1"] == len(k & r)
    assert sorted(a & k & arc) == s["arc_matched_targets"]
    tiers = pd.read_csv(ROOT / "data/splits/arc_target_support_v1.csv")
    assert set(tiers[tiers.support_tier == 2].arc_target) == set(s["arc_matched_targets"])


@needs_outputs
@needs_raw
def test_gene_axes_recompute():
    from virtual_cell.data import scperteval

    def var(n):
        return set(scperteval.read_var_names(RAW / f"{n}_processed_complete.h5ad"))

    a, k, r = var("arch1"), var("kaden25rpe1"), var("replogle22rpe1")
    panel = set(pd.read_csv(ROOT / "data/raw/arc2026/controls/gene_names.csv").iloc[:, 0])
    frozen = (ROOT / "outputs/arc_dry_run_v1/response_genes.txt").read_text().split()
    assert sorted(panel & a & k) == frozen
    axes = json.loads((OUT / "summary.json").read_text())["axes"]
    assert axes == {
        "G_ARC": len(frozen),
        "G_RPE": len(k & r),
        "G_3": len(a & k & r),
        "G_4CTX": 6640,
    }


@needs_outputs
def test_verdict_follows_from_the_saved_tables():
    """Re-apply the predeclared rules to the written tables."""
    s = json.loads((OUT / "summary.json").read_text())
    arc = pd.read_csv(OUT / "arc_supported_target_reliability.csv")
    me = pd.read_csv(OUT / "main_effect_reliability.csv")
    agree = pd.read_csv(OUT / "main_effect_agreement.csv")
    k = arc[arc.dataset == "kaden25rpe1"]
    med, det = k.reliability_centred.median(), k.signal_detected_centred.mean()
    beta = (
        "unsuitable"
        if med < 0.2 and det < 0.5
        else ("suitable" if med >= 0.5 and det >= 0.75 else "partial")
    )
    m = float(
        me[
            (me.dataset == "kaden25rpe1") & (me.axis == "G_ARC") & (me.panel == "all")
        ].spearman_brown.iloc[0]
    )
    mv = "reliable" if m >= 0.8 else ("unreliable" if m < 0.5 else "moderate")
    v = s["verdict"]
    assert v["beta_source"]["verdict"] == beta
    assert v["m_hat_source"]["verdict"] == mv
    assert v["arch1_kaden_main_effect_noise_explanation"] == agree.noise_explanation.iloc[0]
    assert v["case"] in {"A", "B", "C", "D", "E"}


@needs_outputs
def test_arc_target_table_keeps_every_supported_target():
    arc = pd.read_csv(OUT / "arc_supported_target_reliability.csv")
    tiers = pd.read_csv(ROOT / "data/splits/arc_target_support_v1.csv")
    supported = tiers[tiers.support_tier > 0]
    assert set(arc.target) == set(supported.arc_target)
    assert len(arc) == int(supported.n_contexts_perturbed.sum())  # one row per (target, source)
    assert (arc[arc.dataset == "kaden25rpe1"].shape[0]) == 80
