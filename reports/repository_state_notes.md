# Repository state notes — terminology and supersession map

Written 2026-09-25, during the visualization + Kaden source-reliability phase.
Purpose: stop the same word meaning different things in different documents.
Frozen reports are **not edited** to make them consistent; where one carries an
outdated or ambiguous statement it is preserved and clarified here.

## 1. Terminology / supersession map

| term | what it is | where it came from | status |
|---|---|---|---|
| **source agreement** | agreement among the *source-context responses* of a perturbation measured in several sources | zero-shot + transferability-confidence phases | **valid** four-context confidence result (Spearman vs −D 0.55–0.79). Never tested externally. Not used in the Arc model: Tier 2 has only two sources and Tier 1 one |
| **neighbour agreement** | agreement among the STRING/DepMap neighbours standing in for a perturbation measured **nowhere** | unseen-perturbation phase (Tier-0 prior diagnostic) | **refuted externally** (arch1 +0.314 → Feng −0.057 per line). Diagnostic only. `unseen_perturbation_v1_freeze.txt` #5 is superseded by `feng_multicontext_v1_freeze.txt` #5 |
| **scale calibration** (four-context) | one leave-one-source-out least-squares scalar `s` shrinking the centred source mean in LOCO transfer | foundations; pathway residual models | valid. **Same estimator, different values by representation**: gene level (6,640 genes) 0.428–0.469 (`transferability_foundations_v1.md`); Hallmark pathway space 0.52–0.59 (`pathway_residual_model_v1.md`, `research_conclusions_freeze.txt` #1); Reactome 0.53–0.58; external arch1 0.517 and Feng 0.564 (their own gene axes) |
| **Arc tier weights** `w[2]=0.50, w[1]=0.25, w[0]=0` | shrinkage on `beta_hat` in the Arc model, chosen by nested validation on the four public contexts **with each tier's source-count structure** | Arc count-space phase | **frozen Arc parameters**. They play the role scale calibration played in the four-context work but are separate numbers, selected separately |
| **`m_hat` scalar** | one leave-one-source-out scalar rescaling the basal-weighted source main effect (`M3b_basal_shrunk`) | Arc count-space phase | separate again: 0.49–0.87 on the public folds, **0.119** when fitted on Arc's two actual sources (arch1, Kaden) |
| **Tier 0 / 1 / 2** | Arc targets perturbed in 0 / exactly 1 / ≥ 2 public datasets (identifier presence only) | unseen-perturbation phase | frozen: 214 / 79 / 7 = 300 |
| **direct transfer** | transferring a measured perturbation response from the sources | throughout | the only strategy that survived external validation |
| **unseen-perturbation prior** | predicting a perturbation measured nowhere from STRING/DepMap/pathway features | unseen-perturbation phase | failed externally on arch1 and all 19 Feng lines; **do not reopen** |

## 2. Documentation disagreements found by the 2026-09-25 audit

| # | disagreement | resolution |
|---|---|---|
| 1 | `README.md` (research-question list) said raw source agreement "did not replicate externally and is now diagnostic only", while the same README (Transferability section) says only STRING **neighbour** agreement was refuted | The evidence supports the second statement: the only external agreement test (`feng_multicontext_external_validation_v1.md`, research log 2026-09-21) was of neighbour agreement. README question C corrected. |
| 2 | Two "scale calibration" ranges (0.43–0.47 vs 0.52–0.59) | Same estimator in two representations (§1). Both correct. README now says which is which. |
| 3 | README said the scale calibration "is the shrinkage scalar the Arc model uses" | Not literally true: the Arc model uses separately selected tier weights (0.50/0.25) and a separate `m_hat` scalar (0.119). README clarified. |
| 4 | `plans.MD` §62 still calls the Molina & Zhang reproduction a hard gate, and its phase plan (§35–37: representations, interaction model, deep distributional generator) predates the terminations | `plans.MD` is the original specification and is not rewritten. A status banner at its top points here. Binding current state: the M&Z reproduction is **closed** (their data is not public) and replaced by the independent four-context gate (passed); interaction/pathway modelling is **terminated**; no deep generator is justified (count-space report §L12). |
| 5 | Freeze-manifest counts in the research log: "18 manifests / 199 files", "19", "20", then "14 / 220", "15 / 258" | The early counts included the separate data-checksum files (arc2026, scPertEval, acquired, DepMap, MSigDB, STRING, later Feng); from the count-space phase on only `*_freeze.txt` manifests are counted. The Feng-phase "20" fits neither convention (it would be 21). **Current:** 15 freeze manifests, 258 digest entries, 225 unique files; plus 7 data-checksum files over 19 raw files. |
| 6 | Kaden was judged uninformative (`external_validation_v1_freeze.txt` #2; split-half Spearman–Brown 0.115–0.128) yet supplies `beta_hat` for 80 of 86 supported Arc targets and half of `m_hat`; the count-space report never mentions its reliability | Investigated by the new diagnostic, `reports/kaden_source_reliability_diagnostic_v1.md`; see §3. |
| 7 | `arc_count_space_baseline_v1.md` "Validation" section ends with "`uv run pytest` — see below." with nothing below | The report is frozen and stays as written. The number it omitted is in the research log entry for that phase: **491 passed**. |

## 3. Kaden reliability tension

The tension is real and is now measured (`reports/kaden_source_reliability_diagnostic_v1.md`):

- **Per perturbation Kaden is weak**: median split-half reliability 0.17 on the
  Arc response axis, against 0.91 for arch1 and 0.68 for Replogle RPE1. The
  frozen 0.115–0.128 replicates on the 6,640-gene axis (0.123).
- **Averaged, its main effect is moderately reliable** (0.76 over 1,836
  perturbations), but only 0.40 over its 80 Arc targets.
- **The arch1–Kaden main-effect disagreement is not noise**: noise ceiling 0.86,
  observed 0.095. Kaden also fails to agree with the same-cell-line Replogle
  RPE1 screen beyond what noise allows.
- Predeclared verdict: **CASE E (mixed / inconclusive)**. The frozen Arc model
  is unchanged, and no reliability-qualified v2 follows from this diagnostic.
- The historical statements stay as written: `external_validation_v1_freeze.txt`
  #2 ("Kaden is uninformative" as an *external benchmark*) and the count-space
  report's use of Kaden as a *source*. They answer different questions. Both
  now read against this diagnostic.

## 4. Current counts, for reference

- Freeze manifests: 15 (`data/provenance/scperteval/*_freeze.txt`), 258 digests, all verifying.
- Data-checksum files: 7, over 19 raw files, all verifying.
- Arc panel: 300 targets = 7 Tier 2 + 79 Tier 1 + 214 Tier 0; 18,533 genes; 400 cells per perturbation.
- Arc submissions made: **none**. One validated dry-run bundle exists (`outputs/arc_dry_run_v1/`).
