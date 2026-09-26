# Arc submission v1: official hidden-validation result (frozen)

Captured 2026-09-25 (local) from `vcc --json status zYdT8klGWw8UXg3r8KJx`.
Verbatim CLI output: `outputs/arc_submission_v1/official_score_raw.json`; parsed
record with provenance: `outputs/arc_submission_v1/official_score_parsed.json`;
figure: `reports/figures/10_arc_submission_scorecard.{png,svg}`.
Reproduce the freeze: `uv run python scripts/competition_v2/freeze_v1_result.py`.

**This file records the result. It does not edit the predeclared expectations**
(`reports/arc_submission_v1_expectations.md`, committed before submission).

## Record

| field | value |
|---|---|
| submission ID | `zYdT8klGWw8UXg3r8KJx` |
| model name | virtual-cell-generalization baseline v1 (`arc_count_space_baseline_v1`) |
| submission timestamp | 2026-09-26T03:40:47.080101Z |
| status | published, terminal |
| partition / panel | `val` / `vcc2026-val-1` |
| anchor set | `vcc2026-valA-r4+vcc2026-valB-r4+vcc2026-valC-r4` |
| leaderboard rank at capture | 882 |
| package built at git SHA | `a4e9a27ed5effc01ececb282a7ff04ef6f20487d` |
| submitted from git HEAD | `04f463a` (candidate freeze + expectations) |
| submitted `.vcc` SHA-256 | `16a5b17c429582d2c93dc94356cf350d815eb62f8881097aad6441f5a0390910` (recomputed from disk; equals the manifest) |
| source `.h5ad` SHA-256 | `e2aa9acd58eb0060586140849a3389e38b09ebaa0bc30c510e78a87914dfcadb` |

## Scores

| member | official scaled | raw metric | control-emitting reference (frozen) |
|---|---|---|---|
| **Overall** | **−0.0620** | — | −0.311 |
| PDS (`pds_cosine`) | +0.0216 | 0.5096 | 0.000 |
| Expression accuracy (`expr_mse_unbiased_capped_norm`) | 0.0000 (clamped) | 1.0831 | 0.000 |
| DE log-FC (`de_wilcoxon_lfc_nmae`) | −0.0062 | 1.0047 | +0.002 |
| DE direction fidelity (`..._direction_fidelity_yield_raw`) | **−0.3488** | 0.4074 | −1.712 |
| DE direction reach (`..._direction_reach_raw`) | −0.0029 | 0.0765 | −0.080 |
| DE significance overlap (`de_wilcoxon_sig_jaccard`) | −0.0358 | 0.0171 | −0.078 |

Overall equals the mean of the six members to 1e-12.

## Against the predeclared expectations

The bands (`clearly positive ≥ +0.05`, `near baseline`, `below baseline ≤ −0.05`) and
the pattern table were fixed before submission and are applied here unchanged.

| member | expectation (verbatim gist) | band observed | met? |
|---|---|---|---|
| PDS | strongest relative member, small absolute gain | near baseline (+0.022), and the strongest member | **yes** |
| Expression | among the weakest; near or below 0 | near baseline (0, clamped; raw 1.083 > b) | **yes** |
| DE direction fidelity | below baseline, possibly well below | below baseline (−0.349) | **yes** |
| DE direction reach | limited | near baseline (−0.003) | **yes** |
| DE significance overlap | uncertain | near baseline (−0.036) | (no directional claim) |
| DE log-FC | limited, likely near 0 | near baseline (−0.006) | **yes** |
| Overall | risk: toward or below −0.311 | −0.062: well above the control reference, below 0 | partially — the feared floor was not reached |

**Pattern.** No member is clearly positive. PDS and direction fidelity are both near
or below baseline, and every other member is near baseline. This is **pattern 2**
(supported effects barely transfer to Arc's contexts: PDS +0.022) combined with the
predeclared reading of fidelity. Fidelity landed at −0.35, not near −1.7, so the
near-control Tier-0 majority does **not** fully dominate it. Pattern 5 (all six near or
below baseline) also technically holds. Its predeclared licence is "do not increase
model complexity; re-examine the public-to-Arc evidence chain first". The competitive
baseline expansion (`reports/competition_v2/`) does that re-examination.
It adds no capacity to V1.

The result is **frozen**. Nothing in V1 is re-tuned on it.
