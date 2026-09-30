"""Freeze the official VCC hidden-validation result for C1 and draw the V1 -> C1 scorecard.

Two sources, never mixed:

* ``--entry <id>``: ``vcc --json status <id>`` is captured verbatim to
  ``official_score_raw.json``. This is the authoritative path, and it fills raw metrics,
  partition, panel, anchors and timestamp.
* ``--user-reported``: the scaled scores and rank the human read off the leaderboard
  (2026-09-28). The installed ``vcc`` 0.2.0 has no submission-list command, so without an
  entry id the CLI cannot retrieve them. Every field the CLI would supply is recorded as
  ``null``, and the record is marked ``provenance: user-reported``.

The submitted ``.vcc`` checksum is recomputed from disk and must equal the C1 manifest.
Outputs:

* ``outputs/competition_v2/c1_license_clean/official_score_{raw,parsed}.json``
* ``reports/competition_v2/figures/c2_A_official_scorecard.{png,svg}`` and its source CSV

    uv run python scripts/competition_v2/freeze_c1_result.py --entry <id>
    uv run python scripts/competition_v2/freeze_c1_result.py --user-reported
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from virtual_cell.visualization import style  # noqa: E402

BASE = ROOT / "outputs" / "competition_v2" / "c1_license_clean"
V1 = ROOT / "outputs" / "arc_submission_v1" / "official_score_parsed.json"
FIG = ROOT / "reports" / "competition_v2" / "figures"
EXPECTED_VCC_SHA = "fcc4e2508798805d93b9f1bb3a6957ca7cc8fbaa9f8318296161e4184bc43ea7"
CODE_HEAD = "5a2831419ba45547046f0fec84dac86dbc21ec52"  # commit that contains the C1 code
BUILD_HEAD = "d375f935"  # HEAD when the package was built (C1 code then uncommitted)

#: (key, label, official score field, official raw field)
MEMBERS = [
    ("pds", "PDS", "score_pds", "pds_cosine"),
    ("expression_accuracy", "Expression (MSE)", "score_mse", "expr_mse_unbiased_capped_norm"),
    ("de_lfc_accuracy", "DE log-FC (NMAE)", "score_nmae", "de_wilcoxon_lfc_nmae"),
    (
        "de_direction_fidelity",
        "DE fidelity (FID)",
        "score_fid",
        "de_wilcoxon_direction_fidelity_yield_raw",
    ),
    ("de_direction_reach", "DE reach (REACH)", "score_reach", "de_wilcoxon_direction_reach_raw"),
    ("de_significance_overlap", "DE overlap (JAC)", "score_jac", "de_wilcoxon_sig_jaccard"),
]

#: As read off the leaderboard by the human on 2026-09-28 (rounded as displayed).
USER_REPORTED = {
    "rank": "370 / 1207",
    "overall": 0.1394,
    "pds": 0.602,
    "expression_accuracy": 0.057,
    "de_lfc_accuracy": 0.120,
    "de_direction_fidelity": -0.019,
    "de_direction_reach": 0.076,
    "de_significance_overlap": 0.001,
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(16 * 1024**2), b""):
            h.update(block)
    return h.hexdigest()


def parsed_from_raw(raw: dict) -> dict:
    if raw["status"] != "published" or not raw["is_terminal"]:
        raise SystemExit("entry is not a published, terminal score")
    return {
        "provenance": "vcc status --json",
        "entry_id": raw["entry_id"],
        "model_name": raw["model_name"],
        "submission_timestamp_utc": raw["submission_date"],
        "partition": raw["partition"],
        "panel": raw["panel_id"],
        "anchor_set": raw["anchor_version"],
        "rank": raw["rank"],
        "overall": raw["score_avg"],
        "members": {k: {"official_score": raw[s], "raw_metric": raw[r]} for k, _, s, r in MEMBERS},
    }


def parsed_user() -> dict:
    return {
        "provenance": "user-reported (leaderboard, 2026-09-28); NOT retrieved by the CLI. "
        "vcc 0.2.0 has no submission-list command and no C1 entry id is recorded locally. "
        "Re-run with --entry <id> to replace this record.",
        "entry_id": None,
        "model_name": None,
        "submission_timestamp_utc": None,
        "partition": None,
        "panel": None,
        "anchor_set": None,
        "rank": USER_REPORTED["rank"],
        "overall": USER_REPORTED["overall"],
        "members": {
            k: {"official_score": USER_REPORTED[k], "raw_metric": None} for k, *_ in MEMBERS
        },
    }


def scorecard(c1: dict) -> None:
    v1 = json.loads(V1.read_text())
    rows = [{"member": "Overall", "V1": v1["overall"], "C1": c1["overall"]}]
    for key, label, *_ in MEMBERS:
        rows.append(
            {
                "member": label,
                "V1": v1["members"][key]["official_score"],
                "C1": c1["members"][key]["official_score"],
            }
        )
    table = pd.DataFrame(rows)
    table["delta"] = table.C1 - table.V1
    (FIG / "sources").mkdir(parents=True, exist_ok=True)
    table.to_csv(FIG / "sources" / "c2_A_official_scorecard.csv", index=False)

    style.apply()
    fig, ax = plt.subplots(figsize=(style.FULL_WIDTH, 3.6))
    y = np.arange(len(table))
    ax.barh(y - 0.2, table.V1, 0.38, color=style.GREY, label="V1 (official)")
    ax.barh(y + 0.2, table.C1, 0.38, color=style.BLUE, label="C1 license-clean (official)")
    for i, r in table.iterrows():
        ax.text(max(r.C1, 0) + 0.01, i + 0.2, f"{r.C1:+.3f}", va="center", fontsize=7)
    style.reference_line(ax, 0, axis="x")
    ax.set_yticks(y, table.member)
    ax.invert_yaxis()
    ax.set_xlabel("official normalised score (0 = mean-response baseline, 1 = replicate)")
    src = "user-reported" if c1["entry_id"] is None else f"entry {c1['entry_id']}"
    ax.set_title(f"VCC hidden validation: V1 rank 883 → C1 rank {c1['rank']} ({src})", fontsize=9)
    ax.legend(loc="lower right", fontsize=7)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(FIG / f"c2_A_official_scorecard.{ext}", dpi=200)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--entry")
    group.add_argument("--user-reported", action="store_true")
    args = parser.parse_args()
    vcc_path = BASE / "c1_license_clean_val.vcc"
    vcc_sha = sha256(vcc_path)
    if vcc_sha != EXPECTED_VCC_SHA:
        raise SystemExit("C1 .vcc on disk no longer matches the frozen package")
    if args.entry:
        text = subprocess.run(
            ["vcc", "--json", "status", args.entry], check=True, capture_output=True, text=True
        ).stdout
        raw = json.loads(text)
        (BASE / "official_score_raw.json").write_text(json.dumps(raw, indent=2) + "\n")
        parsed = parsed_from_raw(raw)
    else:
        parsed = parsed_user()
    mean6 = np.mean([m["official_score"] for m in parsed["members"].values()])
    parsed.update(
        overall_minus_mean_of_six=float(parsed["overall"] - mean6),
        submitted_package={"path": str(vcc_path.relative_to(ROOT)), "sha256": vcc_sha},
        package_built_at_git_head=BUILD_HEAD,
        c1_code_git_sha=CODE_HEAD,
    )
    (BASE / "official_score_parsed.json").write_text(json.dumps(parsed, indent=2) + "\n")
    scorecard(parsed)
    print(json.dumps({k: parsed[k] for k in ("provenance", "entry_id", "rank", "overall")}))


if __name__ == "__main__":
    main()
