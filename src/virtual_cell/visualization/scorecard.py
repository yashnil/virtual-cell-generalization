"""Schema and validation for the (future) Arc submission scorecard.

The scorecard figure is drawn **only** from a real score file copied from
``vcc status <entry> --json`` after scoring. No placeholder or expected value
is ever plotted: :func:`load_scorecard` refuses a file that is missing,
incomplete, or carries a non-finite score.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

from virtual_cell.visualization.common import ROOT

SCORECARD_PATH = ROOT / "data" / "figure_sources" / "arc_submission_v1_scorecard.json"
SCHEMA_PATH = ROOT / "data" / "figure_sources" / "arc_submission_v1_scorecard.schema.json"
ACCOUNTING_PATH = ROOT / "outputs" / "arc_bridge_v1" / "score_accounting.csv"

#: (json key, label, frozen score_accounting member)
MEMBERS = [
    ("pds", "PDS (perturbation discrimination)", "pds_cosine"),
    ("expression_accuracy", "Expression accuracy", "expr_mse_unbiased_capped_norm"),
    ("de_lfc_accuracy", "DE log-FC accuracy", "de_wilcoxon_lfc_nmae"),
    ("de_direction_fidelity", "DE direction fidelity", "de_wilcoxon_direction_fidelity_yield_raw"),
    ("de_direction_reach", "DE direction reach", "de_wilcoxon_direction_reach_raw"),
    ("de_significance_overlap", "DE significance overlap", "de_wilcoxon_sig_jaccard"),
]
SCORE_KEYS = ["overall"] + [k for k, _, _ in MEMBERS]
TEXT_KEYS = ["partition", "panel", "anchor_set"]

SCHEMA = {
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "title": "Arc submission v1 scorecard (normalized official scores)",
    "description": "Fill ONLY from `vcc status <entry> --json` after real scoring. Scores are "
    "the official normalized values (0 = mean-response baseline, 1 = replicate anchor).",
    "type": "object",
    "required": SCORE_KEYS + TEXT_KEYS,
    "properties": {
        **{k: {"type": "number"} for k in SCORE_KEYS},
        **{k: {"type": "string"} for k in TEXT_KEYS},
        "entry_id": {"type": "string"},
        "model_name": {"type": "string"},
        "scored_utc": {"type": "string"},
        "vcc_package_sha256": {"type": "string"},
        "raw_status_json": {"type": "object"},
    },
}


class ScorecardError(ValueError):
    """The score file is absent or does not describe a real, complete score."""


def validate_scorecard(record: dict) -> dict:
    missing = [k for k in SCORE_KEYS + TEXT_KEYS if k not in record]
    if missing:
        raise ScorecardError(f"missing fields: {missing}")
    for k in SCORE_KEYS:
        v = record[k]
        if isinstance(v, bool) or not isinstance(v, int | float) or not math.isfinite(v):
            raise ScorecardError(f"{k} must be a finite number, got {v!r}")
    for k in TEXT_KEYS:
        if not isinstance(record[k], str) or not record[k].strip():
            raise ScorecardError(f"{k} must be a non-empty string")
    return record


def load_scorecard(path: Path = SCORECARD_PATH) -> dict:
    if not path.exists():
        raise ScorecardError(
            f"{path} does not exist: no real score yet. The scorecard figure is generated "
            "only after `vcc status <entry> --json` returns published scores."
        )
    return validate_scorecard(json.loads(path.read_text()))


def write_schema(path: Path = SCHEMA_PATH) -> Path:
    path.write_text(json.dumps(SCHEMA, indent=2) + "\n")
    return path
