"""Competition-v2 invariants: the frozen official v1 result, the postmortem arithmetic,
and the separation of competition artifacts from frozen research files."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "outputs" / "arc_submission_v1" / "official_score_raw.json"
PARSED = ROOT / "outputs" / "arc_submission_v1" / "official_score_parsed.json"
POSTMORTEM = ROOT / "outputs" / "competition_v2" / "v1_postmortem.json"
V1_TIERS = ROOT / "data" / "splits" / "arc_target_support_v1.csv"
V1_TIERS_SHA = "2280f40b2ce874eef335dab169216d7aad476a887dd49b04845a5a5c4b3ed6f4"
MEMBERS = ["score_pds", "score_mse", "score_nmae", "score_fid", "score_reach", "score_jac"]


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_v1_tier_file_is_untouched():
    assert _sha(V1_TIERS) == V1_TIERS_SHA


def test_vendored_atlasshift_is_not_imported_by_research_package():
    for py in (ROOT / "src" / "virtual_cell").rglob("*.py"):
        for line in py.read_text().splitlines():
            code = line.strip()
            if code.startswith(("import ", "from ")) or "sys.path" in code:
                assert "third_party" not in code and "atlasshift" not in code.lower(), py


@pytest.mark.skipif(not RAW.exists(), reason="official score not captured locally")
def test_official_overall_is_mean_of_six_members():
    raw = json.loads(RAW.read_text())
    assert raw["entry_id"] == "zYdT8klGWw8UXg3r8KJx"
    assert raw["partition"] == "val" and raw["panel_id"] == "vcc2026-val-1"
    mean = sum(raw[m] for m in MEMBERS) / 6
    assert raw["score_avg"] == pytest.approx(mean, abs=1e-12)


@pytest.mark.skipif(not PARSED.exists(), reason="official score not captured locally")
def test_parsed_score_matches_raw_and_manifest():
    raw = json.loads(RAW.read_text())
    parsed = json.loads(PARSED.read_text())
    assert parsed["overall"] == raw["score_avg"]
    assert parsed["members"]["de_direction_fidelity"]["official_score"] == raw["score_fid"]
    manifest = json.loads(
        (ROOT / "outputs" / "arc_submission_v1" / "submission_manifest.json").read_text()
    )
    assert parsed["submitted_package"]["sha256"] == manifest["vcc_package"]["sha256"]


@pytest.mark.skipif(not POSTMORTEM.exists(), reason="postmortem not run locally")
def test_postmortem_arithmetic():
    pm = json.loads(POSTMORTEM.read_text())
    overall = pm["official"]["overall"]
    fid = pm["official"]["members"]["de_direction_fidelity"]
    assert pm["A_fid_over_6"] == pytest.approx(fid / 6)
    assert pm["A_overall_if_fid_were_0"] == pytest.approx(overall - fid / 6)
    pds = pm["3_pds_structure"]
    assert pds["fraction_panel_target_specific"] == pytest.approx(86 / 300)
    assert pds["realistic_ceiling_raw_supported_perfect_tier0_chance"] == pytest.approx(
        (86 + 214 * 0.5) / 300
    )
