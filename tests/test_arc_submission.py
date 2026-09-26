"""Arc submission candidate v1: frozen model identity, manifest consistency,
predeclared expectations, and a scorecard that cannot plot placeholders."""

from __future__ import annotations

import json
import runpy
from pathlib import Path

import pandas as pd
import pytest

from virtual_cell.visualization import common, scorecard

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/arc_submission_v1"
MANIFEST = OUT / "submission_manifest.json"
SHA_FILE = ROOT / "data/provenance/arc_submission_v1_sha256.txt"

needs_manifest = pytest.mark.skipif(not MANIFEST.exists(), reason="submission outputs absent")


def _valid() -> dict:
    rec = {k: 0.1 for k in scorecard.SCORE_KEYS}
    rec.update({"partition": "val", "panel": "vcc2026-val-1", "anchor_set": "official"})
    return rec


def test_scorecard_accepts_a_complete_record():
    assert scorecard.validate_scorecard(_valid())["overall"] == 0.1


@pytest.mark.parametrize("key", scorecard.SCORE_KEYS + scorecard.TEXT_KEYS)
def test_scorecard_rejects_missing_fields(key):
    rec = _valid()
    del rec[key]
    with pytest.raises(scorecard.ScorecardError):
        scorecard.validate_scorecard(rec)


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), None, "0.3", True])
def test_scorecard_rejects_non_finite_or_non_numeric_scores(bad):
    rec = _valid()
    rec["pds"] = bad
    with pytest.raises(scorecard.ScorecardError):
        scorecard.validate_scorecard(rec)


def test_scorecard_members_map_to_frozen_accounting():
    acct = ROOT / "outputs/arc_bridge_v1/score_accounting.csv"
    if not acct.exists():
        pytest.skip("frozen accounting absent")
    members = set(pd.read_csv(acct).member)
    assert {m for _, _, m in scorecard.MEMBERS} == members


def test_no_scorecard_is_drawn_without_a_real_score(tmp_path, monkeypatch):
    if scorecard.SCORECARD_PATH.exists():
        pytest.skip("a real score exists")
    monkeypatch.setattr(common, "FIGURES_DIR", tmp_path)
    runpy.run_path(
        str(ROOT / "scripts/figures/plot_arc_submission_scorecard.py"), run_name="not_main"
    )
    with pytest.raises(scorecard.ScorecardError):
        scorecard.load_scorecard()
    assert not list((ROOT / "reports/figures").glob("10_arc_submission_scorecard.*"))


def test_expectations_predeclare_every_member_and_pattern():
    text = (ROOT / "reports/arc_submission_v1_expectations.md").read_text()
    assert "before any hidden score exists" in text
    for member in (
        "PDS",
        "Expression accuracy",
        "DE direction fidelity",
        "DE direction reach",
        "DE significance overlap",
        "DE log-FC accuracy",
    ):
        assert member in text
    for n in range(1, 7):
        assert f"| **{n}** |" in text


@needs_manifest
def test_manifest_records_the_exact_frozen_model():
    m = json.loads(MANIFEST.read_text())
    assert m["submitted"] is False
    assert m["model_id"] == "arc_count_space_baseline_v1"
    model = m["model"]
    assert model["m_hat"]["estimator"] == "M3b_basal_shrunk"
    assert model["m_hat"]["fitted_scalar"] == pytest.approx(0.11890236823138407, abs=0)
    assert model["tier_weights"] == {"2": 0.5, "1": 0.25, "0": 0.0}
    assert model["tier_counts"] == {"tier_2": 7, "tier_1": 79, "tier_0": 214}
    assert model["generator"] == {
        "name": "G1 control transport",
        "smoothing": 0.5,
        "seed": 20260921,
        "cells_per_perturbation": 400,
    }
    assert (m["cells"], m["genes"], m["perturbations"]) == (360_000, 18_533, 300)
    assert m["contexts"] == ["A", "B", "C"]


@needs_manifest
def test_manifest_matches_the_sha_file_and_files_on_disk():
    m = json.loads(MANIFEST.read_text())
    recorded = {}
    for line in SHA_FILE.read_text().splitlines():
        parts = line.split()
        if len(parts) == 2 and len(parts[0]) == 64:
            recorded[parts[1]] = parts[0]
    assert recorded[m["source_prediction"]["path"]] == m["source_prediction"]["sha256"]
    assert recorded[m["vcc_package"]["path"]] == m["vcc_package"]["sha256"]
    vcc = ROOT / m["vcc_package"]["path"]
    assert vcc.stat().st_size == m["vcc_package"]["bytes"]
    assert common.sha256(vcc) == m["vcc_package"]["sha256"]
    assert m["source_prediction"]["sha256"] == (
        "e2aa9acd58eb0060586140849a3389e38b09ebaa0bc30c510e78a87914dfcadb"
    )


@needs_manifest
def test_candidate_audit_passed_and_matched_frozen_dry_run():
    a = json.loads((OUT / "candidate_audit.json").read_text())
    assert a["failures"] == []
    assert all(a["hard_checks"].values()) and len(a["hard_checks"]) == 12
    assert all(a["local_checks_13"].values()) and len(a["local_checks_13"]) == 13
    assert a["local_checks_equal_frozen"] is True
    assert a["vcc_prep_dry_run"]["identical_to_frozen"] is True
    assert a["statistics"]["nnz"] < a["statistics"]["nnz_cap"]
    pkg = json.loads((OUT / "vcc_prep_package.json").read_text())
    assert pkg["verified_targets"] is True and pkg["dropped"] == [] and pkg["dry_run"] is False
