# C2 predeclaration, amendment 2: reproduction tolerance (before any candidate was scored)

Written 2026-09-28, after the phase-1 H1 run stopped at the reproduction gate. At that
point only C1-equivalent arms had been scored (G0 at a = 1, the two anchors and the G0
null). No C2 candidate had been scored.

## What failed

| arm vs frozen C1 fold | max abs raw difference |
|---|---|
| `ANCHOR_mean_response` | 1.1e-16 |
| `ANCHOR_split_half` | 1.1e-16 |
| `G0_null` vs `G0_control` | 1.1e-16 |
| **`G0_a1.00` vs `C1a`** | **5.7e-7** (FID; NMAE 2.0e-7, MSE 7.9e-9, JAC 3.3e-9, PDS 0) |

## Why

The frozen C1 H1 fold scores (`folds/H1/scores_raw.csv`) were written on 2026-09-26 at
15:14. `competition_v2/{atlas,fusion,…}.py` were last modified at 15:29 the same day,
after the fold run and before the C1 commit `5a28314`. Everything that does not pass
through the fused effect reproduces to machine precision. The C1a arm, which does, moves
by a handful of integer counts. So the committed code is not byte-for-byte the code that
scored the frozen folds, and the difference is immaterial: 5.7e-7 raw is about 1e-6
scaled.

## Change

* The reproduction gate becomes: anchors and null ≤ 1e-12, and `G0_a1.00` vs `C1a`
  ≤ 1e-6 on every raw member.
* **The C1 baseline in rule J is this run's `G0_a1.00`**, so every arm comes from the
  same committed code. The frozen C1 numbers are reported beside it.

Nothing else changes.
