# C1 official hidden-validation result (frozen)

Frozen 2026-09-28 by `scripts/competition_v2/freeze_c1_result.py`. Outputs:

* `outputs/competition_v2/c1_license_clean/official_score_parsed.json`
* figure `reports/competition_v2/figures/c2_A_official_scorecard.{png,svg}`, with source
  `figures/sources/c2_A_official_scorecard.csv`

## Provenance: read this first

**The scores below were reported by the human from the leaderboard. They were not
retrieved by the CLI.**

* `vcc` 0.2.0 has no submission-list or history command. `vcc status` needs an entry id.
* No C1 entry id is recorded anywhere locally: not in the repo, the vcc state file,
  pending uploads or shell history.
* The only programmatic route, a direct `GET /api/cli/submissions` with the stored
  token, was declined by the session's permission policy.
* Separately, `vcc status zYdT8klGWw8UXg3r8KJx` (V1) now returns `not_found`, although
  the same call succeeded on 2026-09-26. V1's official record is the frozen copy in
  `outputs/arc_submission_v1/official_score_raw.json`.

Fields the CLI would supply (entry id, raw metrics, partition, panel, anchor set,
timestamp) are therefore recorded as **pending**. To complete them, run

    uv run python scripts/competition_v2/freeze_c1_result.py --entry <C1 entry id>

This replaces the record with the verbatim `vcc --json status` output and redraws the
figure.

Consistency check on the reported numbers: the mean of the six members is 0.1395, against
the reported Overall 0.1394. They agree to display rounding.

## Record

| field | value |
|---|---|
| submitted | yes (per the human; the leaderboard shows it) |
| entry id | **pending** (not retrievable without it) |
| rank | **370 / 1207** |
| Overall | **0.1394** |
| partition / panel / anchor set | pending; expected `val` / `vcc2026-val-1` / `vcc2026-valA-r4+valB-r4+valC-r4` as for V1, **unverified** |
| timestamp | pending |
| submitted package | `outputs/competition_v2/c1_license_clean/c1_license_clean_val.vcc`, SHA-256 `fcc4e2508798805d93b9f1bb3a6957ca7cc8fbaa9f8318296161e4184bc43ea7` (recomputed from disk at freeze; matches) |
| git | package built at `d375f935` with the C1 code uncommitted; that code is committed as `5a2831419ba45547046f0fec84dac86dbc21ec52`. The package bytes are unchanged since |

## Scaled members, V1 → C1 (official)

| member | V1 | C1 | Δ |
|---|---|---|---|
| **Overall** | −0.0620 | **0.1394** | **+0.2014** |
| PDS | 0.0216 | 0.602 | +0.580 (**27.9×**) |
| Expression (MSE) | 0 | 0.057 | +0.057 |
| DE log-FC (NMAE) | −0.0062 | 0.120 | +0.126 |
| DE fidelity (FID) | −0.3488 | −0.019 | +0.330 |
| DE reach (REACH) | −0.0029 | 0.076 | +0.079 |
| DE overlap (JAC) | −0.0358 | 0.001 | +0.037 |
| rank | 883 | 370 / 1207 | −513 places |

Raw metrics: pending, since only the scaled values were reported. V1's raw metrics are
in the V1 freeze.

## Hidden vs public (C1a, same pipeline)

| member | hidden | H1 fold (local) | K562 fold (local) |
|---|---|---|---|
| PDS | 0.602 | 0.626 | 0.437 |
| MSE | 0.057 | 0 (raw worse than baseline, clamped) | 0 (clamped) |
| NMAE | 0.120 | −0.120 | +0.117 |
| FID | −0.019 | −0.122 | −0.158 |
| REACH | 0.076 | −0.012 | +0.136 |
| JAC | 0.001 | +0.048 | +0.001 |
| Overall | **0.139** | +0.070 | +0.089 |

The hidden Overall (0.139) is **above** both public folds (+0.070 / +0.089). The public
estimate was conservative, which is the opposite of the discount the C1 report warned
about. PDS transferred almost exactly: 0.602 against 0.626 on H1. FID is the member
furthest from its public value (−0.02 hidden vs −0.12 to −0.16 local). That fits the
local anchors having been emitted through the same over-calling G0 generator, so the
local FID scale is miscalibrated.

## Binding interpretation (section B)

* The atlas-transfer backbone is **validated**. PDS rose from 0.0216 to 0.602, so
  perturbation identity is no longer the primary bottleneck.
* The weak members are now expression accuracy (0.057), DE reach (0.076), DE fidelity
  (−0.019) and significance overlap (0.001).
* The source-atlas architecture is **not** altered in C2.

This file is frozen. Later analysis goes in `c2_expression_de_calibration.md`.
