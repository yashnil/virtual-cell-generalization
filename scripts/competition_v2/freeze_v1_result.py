"""Freeze the first official VCC 2026 hidden-validation result for submission v1.

Reads the machine-readable status from ``vcc --json status <entry>`` (or, with
``--from-raw``, the already-captured copy), then writes:

* ``outputs/arc_submission_v1/official_score_raw.json``  — the CLI output, verbatim
* ``outputs/arc_submission_v1/official_score_parsed.json`` — scores + provenance
* ``data/figure_sources/arc_submission_v1_scorecard.json`` — input of Figure 10

The ``.vcc`` checksum is recomputed from disk and must equal the frozen manifest.
Nothing here re-derives, edits, or reinterprets the predeclared expectations.

Reproduce: ``uv run python scripts/competition_v2/freeze_v1_result.py --from-raw``
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "arc_submission_v1"
ENTRY = "zYdT8klGWw8UXg3r8KJx"
SUBMITTED_HEAD = "04f463a"  # HEAD when the package was submitted (expectations commit)

#: (parsed key, official status field, raw-metric status field, frozen accounting member)
MEMBERS = [
    ("pds", "score_pds", "pds_cosine", "pds_cosine"),
    (
        "expression_accuracy",
        "score_mse",
        "expr_mse_unbiased_capped_norm",
        "expr_mse_unbiased_capped_norm",
    ),
    ("de_lfc_accuracy", "score_nmae", "de_wilcoxon_lfc_nmae", "de_wilcoxon_lfc_nmae"),
    (
        "de_direction_fidelity",
        "score_fid",
        "de_wilcoxon_direction_fidelity_yield_raw",
        "de_wilcoxon_direction_fidelity_yield_raw",
    ),
    (
        "de_direction_reach",
        "score_reach",
        "de_wilcoxon_direction_reach_raw",
        "de_wilcoxon_direction_reach_raw",
    ),
    ("de_significance_overlap", "score_jac", "de_wilcoxon_sig_jaccard", "de_wilcoxon_sig_jaccard"),
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(16 * 1024**2), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--from-raw", action="store_true")
    args = parser.parse_args()
    raw_path = OUT / "official_score_raw.json"
    if args.from_raw:
        raw = json.loads(raw_path.read_text())
    else:
        text = subprocess.run(
            ["vcc", "--json", "status", ENTRY], check=True, capture_output=True, text=True
        ).stdout
        raw = json.loads(text)
        raw_path.write_text(json.dumps(raw, indent=2) + "\n")
    if raw["status"] != "published" or not raw["is_terminal"]:
        raise SystemExit(f"entry {ENTRY} is not a published, terminal score")

    manifest = json.loads((OUT / "submission_manifest.json").read_text())
    vcc_path = ROOT / manifest["vcc_package"]["path"]
    vcc_sha = sha256(vcc_path)
    if vcc_sha != manifest["vcc_package"]["sha256"]:
        raise SystemExit("submitted .vcc on disk no longer matches the frozen manifest")

    members = {
        key: {"official_score": raw[score], "raw_metric": raw[metric], "member": member}
        for key, score, metric, member in MEMBERS
    }
    overall = raw["score_avg"]
    mean_of_members = sum(m["official_score"] for m in members.values()) / len(members)
    parsed = {
        "entry_id": raw["entry_id"],
        "model_name": raw["model_name"],
        "submission_timestamp_utc": raw["submission_date"],
        "partition": raw["partition"],
        "panel": raw["panel_id"],
        "anchor_set": raw["anchor_version"],
        "leaderboard_rank_at_capture": raw["rank"],
        "overall": overall,
        "overall_equals_mean_of_six": abs(overall - mean_of_members) < 1e-12,
        "members": members,
        "submitted_package": {
            "path": manifest["vcc_package"]["path"],
            "sha256": vcc_sha,
            "sha256_matches_manifest": True,
            "bytes": manifest["vcc_package"]["bytes"],
        },
        "source_prediction": manifest["source_prediction"],
        "package_built_at_git_sha": manifest["git_sha"],
        "submitted_from_git_head": SUBMITTED_HEAD,
        "model_id": manifest["model_id"],
        "expectations": "reports/arc_submission_v1_expectations.md (unmodified)",
    }
    (OUT / "official_score_parsed.json").write_text(json.dumps(parsed, indent=2) + "\n")

    scorecard = {
        "overall": overall,
        **{k: v["official_score"] for k, v in members.items()},
        "partition": raw["partition"],
        "panel": raw["panel_id"],
        "anchor_set": raw["anchor_version"],
        "entry_id": raw["entry_id"],
        "model_name": raw["model_name"],
        "scored_utc": raw["submission_date"],
        "vcc_package_sha256": vcc_sha,
        "raw_status_json": raw,
    }
    (ROOT / "data" / "figure_sources" / "arc_submission_v1_scorecard.json").write_text(
        json.dumps(scorecard, indent=2) + "\n"
    )
    print(json.dumps({k: parsed[k] for k in ("entry_id", "overall", "partition", "panel")}))


if __name__ == "__main__":
    main()
