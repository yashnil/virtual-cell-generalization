# X-Atlas/Orion permission status

| field | value |
|---|---|
| dataset | X-Atlas/Orion, HCT116 and HEK293T genome-wide CRISPRi FiCS Perturb-seq |
| rights holder | Xaira Therapeutics (https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion) |
| license | CC BY-NC-SA 4.0 (Hugging Face tag `license:cc-by-nc-sa-4.0`, re-checked 2026-09-26; local `data/raw/competition_v2/xatlas_orion/LICENSE.md`) |
| local copy | `slaf-project/X-Atlas-Orion` @ `598aa5442a4ee9a8f1a383a843ab5ce142904d04`, 8,689 files verified against LFS etags |
| **status** | **PENDING** |
| last updated | 2026-09-26 |

Allowed statuses: **PENDING**, **APPROVED**, **DENIED**.

## Rules

* The status becomes **APPROVED** only on explicit written permission from Xaira
  Therapeutics that covers use in a prize-bearing Virtual Cell Challenge 2026 entry,
  including the ShareAlike question for submitted predictions. The permission text,
  sender, date and channel must be recorded below.
* Public accessibility, the Challenge FAQ, or an informal reply are **not** permission.
* Until APPROVED, no X-Atlas response or statistic may enter a submission candidate.
  `virtual_cell.competition_v2.licensing` enforces this in code (`HCT116`, `HEK293T` =
  `BLOCKED_PENDING_PERMISSION`).
* If permission arrives during a running experiment, it is recorded here with its date.
  The running license-clean experiment is **not** altered. Adding X-Atlas back is a
  separate, later change.

## Log

| date | event | evidence |
|---|---|---|
| 2026-09-26 | Status set to PENDING at the start of C1. No request or reply is recorded in this repository. | — |
| 2026-09-26 | No permission received during C1. C1 was built without X-Atlas. | — |

## What permission would be worth (public folds; details in `license_clean_c1_v1.md` §13)

Source: `outputs/competition_v2/c1_license_clean/xatlas_ablation.csv`. Both arms use the same C0
backbone and weights, and the same folds.

| fold | coverage without → with X-Atlas | sources per target | PDS (raw) | local Overall |
|---|---|---|---|---|
| H1 (150 targets) | 141 → 149 | 1.61 → 3.53 | 0.8277 → 0.8270 | +0.073 → +0.081 |
| K562 (390 targets) | 355 → 389 | 1.20 → 3.13 | 0.6659 → 0.6915 | +0.091 → +0.119 |
| CD4, effect-level (346 targets) | 328 → 346 | 1.27 → 3.21 | 0.6989 → 0.6726 | — |

On the Arc panel, 13 of the 300 targets have no license-clean evidence (all 13 are
covered by X-Atlas).

Summary: +0.018 mean local Overall. PDS is mixed across the three folds (flat, up,
down). The largest single effect is coverage of the 13 unsupported Arc targets.
