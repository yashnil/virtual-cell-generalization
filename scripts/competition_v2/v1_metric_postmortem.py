"""V1 official-score postmortem (sections 2 and 3 of competitive baseline expansion v2).

Pure arithmetic on the frozen official result and the published cell-eval2 anchor
ranges (``vcc2026-metrics.md`` section 8, cell-eval2 @ 5e64833). No model is run
and no expectation is edited.

Writes ``outputs/competition_v2/v1_postmortem.json``.

Reproduce: ``uv run python scripts/competition_v2/v1_metric_postmortem.py``
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PARSED = ROOT / "outputs" / "arc_submission_v1" / "official_score_parsed.json"
OUT = ROOT / "outputs" / "competition_v2"

N_PANEL = 300
TIERS = {"2": 7, "1": 79, "0": 214}
N_SUPPORTED = TIERS["2"] + TIERS["1"]

#: published anchor ranges over contexts A/B/C (lo, hi)
PDS_BASELINE = 0.5
PDS_REPLICATE = (0.927, 0.984)
FID_BASELINE = (0.505, 0.522)
FID_REPLICATE = (0.795, 0.832)
FID_HALF_SPLIT_RAW = (0.41, 0.46)  # "a genuine half of the reference", spec section 4
FID_HALF_SPLIT_SCALED = (-0.34, -0.15)
#: "scaled value of an exact reproduction of the reference" (spec section 8)
EXACT_REPRODUCTION = {
    "pds": (1.03, 1.17),
    "expression_accuracy": (1.0, 1.0),
    "de_direction_fidelity": (1.51, 1.72),
    "de_direction_reach": (1.02, 1.05),
    "de_significance_overlap": (2.48, 2.85),
    "de_lfc_accuracy": (1.58, 1.75),
}


def main() -> None:
    parsed = json.loads(PARSED.read_text())
    members = {k: v["official_score"] for k, v in parsed["members"].items()}
    raw = {k: v["raw_metric"] for k, v in parsed["members"].items()}
    overall = parsed["overall"]

    # ---- A. contribution of each member to Overall -------------------------
    contribution = {k: v / 6 for k, v in members.items()}
    fid_zero = overall - contribution["de_direction_fidelity"]
    others = [v for k, v in members.items() if k != "de_direction_fidelity"]
    without_fid_mean_of_five = sum(others) / 5

    # ---- B. direction precision vs yield (bounds from the single raw value) --
    # F = precision * yield with both in [0, 1], so each factor is at least F.
    fid_raw = raw["de_direction_fidelity"]
    fid_bounds = {
        "raw_fid": fid_raw,
        "precision_lower_bound": fid_raw,
        "yield_lower_bound": fid_raw,
        "if_precision_is_chance_0.5_then_yield": fid_raw / 0.5,
        "if_yield_is_1_then_precision": fid_raw,
        "half_split_reference_raw_range": FID_HALF_SPLIT_RAW,
        "half_split_reference_scaled_range": FID_HALF_SPLIT_SCALED,
        "within_half_split_raw_range": FID_HALF_SPLIT_RAW[0] - 0.01
        <= fid_raw
        <= FID_HALF_SPLIT_RAW[1],
        "scaled_anchor_consistency": [
            (fid_raw - b) / (r - b) for b, r in zip(FID_BASELINE, FID_REPLICATE, strict=True)
        ],
    }

    # ---- 3. structural PDS limit -------------------------------------------
    frac_supported = N_SUPPORTED / N_PANEL
    # Tier-0 targets share one expected profile; sampling noise alone orders them,
    # so their expected normalized rank of the true match is 0.5.
    realistic_ceiling_raw = (N_SUPPORTED * 1.0 + TIERS["0"] * 0.5) / N_PANEL
    # Absolute bound: one shared Tier-0 vector ranks all 214 Tier-0 truths first.
    tier0_best = 1 - ((TIERS["0"] - 1) / 2) / (N_PANEL - 1)
    absolute_ceiling_raw = (N_SUPPORTED * 1.0 + TIERS["0"] * tier0_best) / N_PANEL

    def scale(u: float) -> list[float]:
        return [(u - PDS_BASELINE) / (r - PDS_BASELINE) for r in PDS_REPLICATE]

    implied_supported_pds = 0.5 + (raw["pds"] - 0.5) * N_PANEL / N_SUPPORTED
    pds = {
        "tier_counts": TIERS,
        "fraction_panel_target_specific": frac_supported,
        "fraction_panel_shared_context_response": TIERS["0"] / N_PANEL,
        "realistic_ceiling_raw_supported_perfect_tier0_chance": realistic_ceiling_raw,
        "realistic_ceiling_scaled_range": sorted(scale(realistic_ceiling_raw)),
        "absolute_ceiling_raw_tier0_best_ordering": absolute_ceiling_raw,
        "absolute_ceiling_scaled_range": sorted(scale(absolute_ceiling_raw)),
        "observed_raw": raw["pds"],
        "observed_scaled": members["pds"],
        "implied_mean_pds_of_86_supported_targets_if_tier0_is_0.5": implied_supported_pds,
        "global_amplitude_invariance": "cosine distance is invariant to a positive global "
        "scalar on the predicted effect, so scalar calibration changes PDS only through "
        "count-sampling noise",
    }

    # Oracle-on-supported ceiling for every member (per-perturbation means are linear
    # in the per-target scaled value; Tier-0 held at the mean-response baseline = 0).
    oracle = {
        k: [frac_supported * lo, frac_supported * hi] for k, (lo, hi) in EXACT_REPRODUCTION.items()
    }
    oracle["pds"] = sorted(scale(realistic_ceiling_raw))
    oracle_overall = [sum(v[i] for v in oracle.values()) / 6 for i in (0, 1)]

    result = {
        "official": {"overall": overall, "members": members, "raw": raw},
        "A_contribution_to_overall": contribution,
        "A_fid_over_6": contribution["de_direction_fidelity"],
        "A_overall_if_fid_were_0": fid_zero,
        "A_mean_of_other_five": without_fid_mean_of_five,
        "A_fid_share_of_negative_overall": contribution["de_direction_fidelity"] / overall,
        "B_fid_anatomy_bounds": fid_bounds,
        "3_pds_structure": pds,
        "3_oracle_supported_ceiling_per_member_scaled": oracle,
        "3_oracle_supported_ceiling_overall_range": oracle_overall,
        "anchors_source": "cell-eval2 5e64833 docs/vcc2026_metrics/vcc2026-metrics.md s8",
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "v1_postmortem.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
