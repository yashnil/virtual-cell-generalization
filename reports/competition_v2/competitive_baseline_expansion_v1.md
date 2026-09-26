# Arc competitive baseline expansion (v2 phase), report 1

**Competition track.** Date: 2026-09-26. No submission was made in this phase. The
second submission slot is unspent.

The frozen research conclusions are unchanged, and nothing here contradicts them. The
failed hypothesis was that *annotation-derived priors* (STRING / DepMap / pathways)
predict an **unseen** perturbation's effect. This phase expands **direct**
perturbation evidence, which tests a different hypothesis (§15).

All 16 frozen manifests verify with 0 mismatches at the end of the phase. The frozen v1
tier file (`data/splits/arc_target_support_v1.csv`) is unchanged; a test pins its hash.

| artifact | path |
|---|---|
| official v1 result (frozen) | `reports/arc_submission_v1_result.md`, `outputs/arc_submission_v1/official_score_{raw,parsed}.json`, `reports/figures/10_arc_submission_scorecard.{png,svg}` |
| metric postmortem | `outputs/competition_v2/v1_postmortem.json` (`scripts/competition_v2/v1_metric_postmortem.py`) |
| public FID anatomy + amplitude ladder | `outputs/competition_v2/public_fid_anatomy/` (`run_public_fid_anatomy.py`, `analyse_public_fid_anatomy.py`) |
| AtlasShift method audit | `reports/competition_v2/atlas_shift_method_audit.md` |
| AtlasShift source audit | `reports/competition_v2/atlas_shift_source_audit.md` |
| coverage matrix | `data/splits/arc_target_support_competition_v2.csv`, `outputs/competition_v2/target_coverage_summary.json` |
| C0 bundle (packaged, **not submitted**) | `outputs/competition_v2/atlasshift_c0/` + `c0_manifest.json` |
| C0 public benchmark (H1 held out) + promoter ablation | `outputs/competition_v2/atlasshift_public_h1/` |
| V1 vs C0 mechanistic comparison | `outputs/competition_v2/mechanistic_comparison/` |
| STATE feasibility | `reports/competition_v2/state_feasibility_audit.md` |
| vendored external code | `third_party/atlasshift/` (MIT, commit `d24ce4f`, byte-identical), `third_party/README.md` |
| competition namespace | `src/virtual_cell/competition_v2/`, `scripts/competition_v2/` |

---

## 1. Why was V1's score negative?

**Direction fidelity alone.** FID/6 = **−0.0581**, which is 94 % of the −0.0620 Overall. The
other five members sum to −0.023, so with FID at 0 and everything else unchanged the
Overall would have been **−0.0039**.

| member | scaled | contribution (/6) |
|---|---|---|
| PDS | +0.0216 | +0.0036 |
| Expression | 0.0000 (clamped; raw 1.083) | 0 |
| NMAE | −0.0062 | −0.0010 |
| **FID** | **−0.3488** | **−0.0581** |
| REACH | −0.0029 | −0.0005 |
| JAC | −0.0358 | −0.0060 |

FID is `k / max(n_pred, n_real)`, which equals directional precision × min(1, yield). A
negative score therefore does not by itself mean wrong signs. The official specification
(cell-eval2 `5e64833`, `vcc2026-metrics.md` §4) states: *"a genuine half of the reference
scores −0.15 to −0.34 here, because a half calls fewer genes and its raw value, 0.41–0.46,
falls below the mean-response arm's 0.505."*

V1's raw FID was **0.407**, the level of a half-depth replicate of the truth. From one raw
number, precision and yield are each bounded below by 0.407. If precision were at chance
(0.5), yield would be 0.81.

Read jointly with the other members: REACH ≈ 0, JAC slightly negative, NMAE ≈ 0 and PDS
slightly positive. That is a submission whose predicted effects are **real but weak**. It
calls too few genes (low FID, JAC and REACH), its fold-change sizes are near the null
(NMAE ≈ 1), and it carries a small amount of perturbation identity (PDS). §8 settles
sign versus yield on public data.

## 2. Why is V1 structurally noncompetitive even if FID were fixed?

* **With FID at 0 the Overall is −0.004.** Fixing FID completely gets V1 to the baseline,
  not above it.
* **Only 86 / 300 (28.7 %) targets carry any target-specific information.** The other
  214 share one expected profile, the shrunk `m_hat`.
* **PDS ceiling.** Suppose the supported targets were perfect and Tier-0 stayed at
  chance. Raw PDS would be (86 + 214 · 0.5)/300 = **0.643**, which scales to
  **0.30–0.34** with the official replicate range 0.927–0.984. Even if one shared vector
  somehow ranked all 214 Tier-0 truths first, the bound is raw 0.746 (scaled 0.51–0.58).
  Observed: raw **0.5096**. That implies a mean PDS of only **0.534** for the 86
  supported targets, so direct transfer from other contexts barely discriminates on A/B/C.
* **Scalar calibration cannot move PDS.** Cosine distance is invariant to a positive
  global scalar. On public data the ladder *lowers* PDS as amplitude grows (§9),
  because the shared `m_hat` direction grows too.
* **Oracle ceiling of the architecture.** Supported targets perfect, Tier-0 at the
  baseline: Overall 0.41–0.46. The architecture is not the binding limit. The binding
  limits are that 71 % of targets have no evidence, and that the evidence V1 had for
  the other 29 transferred poorly.

**Answer: no.** Scalar calibration cannot make V1 competitive.

## 3. How many Arc targets gain direct evidence from the expanded atlas? (primary result)

Targets counted by the number of independent **cell contexts** with ≥ 20 directly
perturbed cells (CD4 additionally requires the publisher quality flags). The v1 column
restates the frozen tiers.

| contexts with direct evidence | V1 (frozen) | competition v2 |
|---|---|---|
| 0 | **214** | **0** |
| 1 | 79 | 1 |
| 2 | 7 | 18 |
| 3+ | 0 | **281** |

**All 300 targets now have direct evidence, and 281 have it in three or more contexts.**
All 214 former Tier-0 targets are rescued, by X-Atlas HEK293T (208), X-Atlas HCT116 (204),
K562 GWPS (183), CD4 (179) and H1 2025 (6).

| source (context) | Arc targets measured | usable |
|---|---|---|
| X-Atlas/Orion HEK293T | 300 | 294 |
| X-Atlas/Orion HCT116 | 300 | 290 |
| K562 GWPS | 272 | 269 |
| CD4 DE (Rest / Stim8hr / Stim48hr) | 291 | 251 |
| Kaden 2025 RPE1 (v1) | 80 | 80 |
| H1 2025 full (train + validation + test) | 25 | 25 |
| `arch1` = H1 2025 train (v1) | 13 | 13 |

The old figure of "86/300 public coverage" held only for the seven scPertEval datasets v1
audited. **Strategic consequence.** The VCC FAQ (verbatim, site bundle, 2026-09-25) says
*"the validation and final test rounds use different panels."* The D/E/F targets are
unknown until 2026-10-22. Only **genome-wide** sources guarantee coverage of an unknown
panel: K562 GWPS (9,866 targets), X-Atlas HCT116 (18,294), HEK293T (18,312) and CD4
(11,526 contrasts). The raw data must be kept, because re-preparing for a new panel
reads it.

## 4. How much did the full VCC 2025 data change coverage?

A little. `arch1` is exactly the 150 H1 training perturbations. The complete public 2025
release (train 150 + validation 50 + test 100 = 300 targets, all raw counts, 18,080
genes) contributes **25** Arc 2026 targets: train 13, validation 4 (ACLY, HDAC8, NFE2L1,
SLIRP) and test 8 (ADNP, ANKZF1, BRPF1, MED15, MED25, PLAGL2, ZFP62, ZNF32). That is
**+12** over `arch1`, 6 of them former Tier-0 targets.

The 2025 data matters more as a **public held-out fold** than as coverage. It is Arc's
own protocol, which makes it the most Arc-like public benchmark available (§6, §7).

## 5. Was AtlasShift reproduced exactly enough to package?

**Yes. It is packaged as `C0_ATLASSHIFT_REPRODUCTION` and not submitted.**

* **Upstream code ran unmodified.** Commit `d24ce4f`, vendored byte-identically, in its
  own Python 3.13 environment with `requirements.txt` pins exactly (numpy 2.5.2, pandas
  3.0.5, anndata 0.13.3.post0, …) plus `slafdb` 0.3.2 and `pylance` 12.0.0. The
  pipeline was `prepare.py` (all 6 jobs) → `predict.py` → `compact.py` → `pack.py`.
* **Every URL source matched `sources.json` byte count and SHA-256**, verified by our
  downloader and again by upstream `prepare.py`.
* **X-Atlas was pinned.** Upstream scans the *live* Hugging Face dataset. We pinned
  revision `598aa544…` locally (LFS hashes verified) and passed it via `--atlas-root`.
  This is the one deliberate deviation from upstream, made for reproducibility. The
  results could differ from the authors' own run only if they scanned a different
  revision.
* **Validation:** all 13 frozen v1 local checks, an exact target-set match, no stored
  zeros, and `vcc prep --dry-run` exit 0 (`verified_targets: true`,
  `normalization: counts-preserved`, nothing dropped, no genes reordered). `pack.py`
  packaged through `vcc.prep.run_prep` (vcc-cli 0.2.0) with every official check.

| artifact | SHA-256 |
|---|---|
| `prediction.h5ad` (360,000 × 18,533; 2,437,685,571 nnz, over 2³¹ and under the 4.75e9 cap) | `155e2440aed8add083e2259f288c19e66c77106177d3305235f93c577bc9d51f` |
| `prediction_compact.h5ad` | `d47a15187e299a4561a9a75eed4876594d5495fb06461d473313801b11ea7dea` |
| `c0_atlasshift_reproduction_val.vcc` (3.35 GB) | `bd33b1adb7e128518206bf22a7495cbcc2e5dd8799f38caceacc89e79f3ea0ae` |

Code, source and prepared-statistics checksums: `outputs/competition_v2/atlasshift_c0/c0_manifest.json`.

**Not established.** Whether C0 scores the reported 0.1546 on the hidden panel today. That
figure remains an external reference point. Upstream records no seed-level or
environment-level provenance of the authors' own bundle to compare against.

## 6. What does AtlasShift do that V1 does not?

Measured on the two A/B/C bundles' emitted counts (`mechanistic_comparison/`):

| | V1 | C0 |
|---|---|---|
| A. targets with reproducible target-specific signal (cross-context reproducibility above the V1 Tier-0 noise null, p99 = 0.015) | **89** (all 86 supported + 3 at the null rate) | **300** |
| median same-target cross-context cosine of the deviation | 0.003 | **0.91** |
| … for V1 Tier 2 / 1 / 0 targets | 0.16 / 0.07 / −0.00 | 0.91 / 0.91 / 0.91 |
| B. pseudobulk effect norm, median (5–95 %) | 3.90 (3.74–4.17), almost all sampling noise | 1.44 (1.03–2.16), almost all signal |
| C. deviation effective rank (A/B/C) | 288–292 (noise is full-rank) | 235–237 (structured) |
| F. context adaptation | `m_hat` per context; beta is context-free | **none**: the same fused effect is multiplied onto each context's control |
| G. per-cell structure vs real controls (context A) | library 20,084 vs 20,109; genes detected 6,101 vs 6,147; heterogeneity 33.3 | library 20,037; genes detected **7,240**; heterogeneity **21.9** (4-cell averaged template) |

Four differences follow.

1. **Evidence.** C0 has a direct, multi-context measurement of every target. V1 has
   one for 86.
2. **No context mean shift.** C0 adds no `m_hat`. Per-source centering removes each
   atlas's generic response, and the target context's control is the base. V1 spent a
   parameter on the context main effect and shrank it to 0.119, so Tier-0 was
   near-control.
3. **Emission.** C0's pseudobulk is pinned to its expectation. V1's G1 transport emits the
   full sampling noise of single resampled controls. That noise cancels in PDS (PDS was
   unchanged between the two generators on public data, §8). It does change the DE yield
   members, and it made V1's own supported signal almost invisible in the emitted counts.
4. **Unrealistic cells.** C0's cells detect ~18 % more genes than real cells and are less
   heterogeneous. The metric does not penalise this directly, but it produces spurious
   Wilcoxon calls (§8).

**E. DE yield on public held-out data (H1 2025 train, 150 targets, C0 rebuilt without the
H1 source).** C0 calls a median of 5,450 genes against 1,296 real ones (yield 0.97). Its
precision is **0.540**, lower than the oracle mean-response arm's 0.587. So its FID of
0.525 sits *below* that arm's 0.572. The DE members are at or just below baseline: FID
−0.14, REACH +0.01, JAC +0.03, NMAE −0.07 (local anchors). **C0's public advantage is
almost entirely PDS:** raw 0.827, local scaled +0.655.

**Why a broader atlas method can beat our conservative model.** PDS is the only member
that rewards knowing *which* perturbation is which, and it can only be earned with
target-specific evidence. C0 has that evidence for all 300 targets and carries it
cleanly into counts. V1 has it for 86 and buried it in emission noise. The DE members sit
near baseline for both methods.

## 7. Does the promoter-neighbour correction matter publicly?

The rule (GENCODE v47): neighbours with TSS within 5 kb of the target TSS are capped at
`control × (0.15 + 0.85 · clip(log10(max(d, 500)/500), 0, 1))`. On the Arc panel it touches
**74 targets and 75 neighbour genes**, 40 of them within 500 bp. The rationale is that
dCas9-KRAB silences the shared or adjacent promoter chromatin, so this is an **assay**
effect. There is no leakage risk.

Ablation on H1 (distance threshold untouched): 47 pairs on that panel, 44 entries capped.

| | PDS raw | FID | local Overall |
|---|---|---|---|
| C0 | **0.827** | 0.5246 | **+0.081** |
| C0 without the promoter prior | 0.809 | 0.5248 | +0.076 |

**Small, positive and free.** It gives +0.018 raw PDS and +0.005 local Overall, with no
effect on DE members. That is one fold, so treat it as directional.

## 8. Is V1's FID failure sign error or low DE yield?

**Low yield (B), not wrong signs (A).** This is shown on four public leave-one-context-out
folds (K562, RPE1, HepG2, Jurkat; 300 perturbations × 400 cells each; frozen G1 generator),
where the mixed Arc-prevalence panel reproduces the hidden failure level: raw FID
0.30–0.44 against the hidden 0.407.

| arm (mean over 4 folds) | FID | precision where called | yield | n_pred / n_real (pooled) |
|---|---|---|---|---|
| Tier 2 | 0.522 | 0.681 | 0.786 | 0.63 |
| Tier 1 | 0.420 | 0.665 | 0.662 | 0.34 |
| **Tier 0** | **0.340** | 0.618 | 0.574 | **0.18** |
| oracle mean response (local b) | 0.422 | 0.686 | 0.659 | 0.39 |
| split half (local r) | 0.787 | 0.938 | 0.841 | 0.99 |

* **Precision is 0.62–0.68 in every tier**, well above chance, so the signs are mostly
  right.
* **Yield collapses with effect magnitude.** On the mixed panel, for the strongest-effect
  tercile, Tier-0 yield is **0.13** (FID 0.085) against 0.95 for the weakest tercile.
  V1 misses big responses because it predicts almost nothing for Tier-0.
* **Most of the remaining yield gap is the generator.** Same V1 effects, K562 fold, a = 1:
  G1 → dual-moment raises yield from 0.58 to 0.92 and FID from 0.353 to 0.559. PDS is
  unchanged (0.582 vs 0.580).
* **Caveat, and a correction to my own first reading.** The dual-moment generator also
  "calls" a median of 3,837 genes when emitting *unperturbed* controls (H1), at chance
  precision. That comes from its averaged-template cell structure. So generator-driven
  yield raises FID to about baseline, never materially above it.

Frozen-harness defect, recorded rather than fixed: `score_generator` in the frozen
`scripts/run_count_generator_benchmark.py` passes the **full-axis** target-gene index to
the DE members. Their tables live on the control-gated tested axis, so `_drop_target`
removes an unrelated gene and keeps the target. The impact is one gene per perturbation
out of ~370 reference calls, immaterial to every conclusion here. Both the frozen v1
count-space numbers and the §8/§9 anatomy carry it. The new H1 harness corrects it.

## 9. Does amplitude calibration help V1 publicly?

**Only through DE yield, and at PDS's expense. It cannot repair Tier 0.** The ladder was
predeclared, applied to the frozen V1 direction on the mixed panel, and averaged over the
four folds with local anchors. Nothing was selected.

| a | PDS | FID | REACH | JAC | NMAE | local Overall |
|---|---|---|---|---|---|---|
| 0.5 | 0.147 | −0.475 | −0.269 | −0.139 | −0.251 | −0.165 |
| 0.75 | 0.137 | −0.360 | −0.235 | −0.109 | −0.214 | −0.130 |
| **1.0** | 0.131 | −0.213 | −0.207 | −0.072 | −0.181 | −0.090 |
| 1.25 | 0.126 | −0.096 | −0.181 | −0.034 | −0.149 | −0.056 |
| 1.5 | 0.113 | +0.006 | −0.161 | −0.005 | −0.130 | −0.029 |
| 2.0 | 0.104 | +0.151 | −0.124 | +0.044 | −0.105 | +0.012 |
| 3.0 | 0.099 | +0.311 | −0.091 | +0.103 | −0.125 | +0.049 |

Expression is clamped at 0 throughout. The raw expression error grows with a, from 1.01 to
2.24 on K562. Precision stays flat at 0.60–0.64; yield rises from 0.48 to 0.89.

Local anchors were emitted through G1, so the local b is itself yield-limited: 0.42
against the official 0.505–0.522. The local scaled DE values are therefore **optimistic**
relative to the official scale. **Calibration is not the v2 route.**

## 10. Is STATE inference feasible?

**Pretrained STATE: feasible to run, useless for this panel.** From the audit,
spot-checked here against `_infer.py` L798–806:

* Every public ST checkpoint encodes perturbations one-hot, and **0 of the 300 targets**
  are in the Replogle map.
* Unknown targets silently fall back to the control one-hot vector. The output is
  roughly the context-conditioned control, scaled score ≈ 0.
* Checkpoints emit 2,000 or 6,546 genes of clipped log-normalised expression, not raw
  counts.
* License: the model is non-commercial (the Arc State Model License), and the VCC FAQ
  allows checkpoints for non-commercial entrants.

The only route to all 300 targets is **training a new ESM2-featurised ST** on a CUDA GPU:
≈ 9 T4-hours per 40k-step run (2025 Colab log), $5–15 per run. Inference is minutes on an
M1 CPU. That would be a new model trained on the same kind of direct data C0 already
uses directly. Plan and commands: `reports/competition_v2/state_feasibility_audit.md`.

## 11. Which route has the highest expected value?

| route | evidence | expected value |
|---|---|---|
| V1 calibration | best public ladder point gives local Overall +0.05, with PDS falling and Tier-0 unrepairable | **lowest** |
| AtlasShift reproduction (C0) | 300/300 targets with reproducible signal; public H1 PDS 0.827 (local +0.655); reported hidden 0.155; covers any D/E/F panel from genome-wide sources | **highest near-term** |
| STATE | pretrained scores ≈ 0; trained ESM2-ST is a multi-week GPU project with uncertain gain over direct fusion when direct data exists for every target | low near-term |
| eventual hybrid | C0 backbone plus a validated improvement | **highest ceiling, conditional** on public-fold evidence |

## 12. What exact candidate should be built next?

One candidate: **C1 = the direct-evidence atlas fusion backbone, rebuilt as our own code,
with one predeclared, research-derived change.**

1. **Backbone = C0, unchanged.**
   * Per-source mean-centred effects from K562 GWPS, X-Atlas HCT116/HEK293T, H1 2025
     full and CD4.
   * Mass-based shrinkage (100k pseudocounts), two-space fusion (amplitude 0.6 / 0.3,
     clip 3), promoter prior, dual-moment integer generator.
   * Reimplemented in `src/virtual_cell/competition_v2/`, tested for **numerical
     equality** against the vendored C0 bundle before any change. That keeps the provenance
     boundary (MIT attribution) explicit, and lets the pipeline be re-prepared for the
     D/E/F panel on 2026-10-22 without editing vendored files.
2. **The single change: agreement-based per-target shrinkage.** Scale each target's
   fused effect by the raw cross-source agreement of its per-source effects. Our frozen
   Question-C result found raw source agreement to be the best transfer-quality predictor
   in all four held-out contexts. It was never usable on Arc because v1 targets had one
   or two sources; 281 targets now have three or more. This applies an existing finding
   and invents no new prior.
3. **Selection happens on public folds only, predeclared before running:**
   * leave-one-atlas-out: H1 2025 (all 300 targets, from train + validation + test),
     HCT116, HEK293T and K562 held out in turn;
   * local cell-eval2 with the corrected target indexing and local anchors;
   * C1 must beat C0 on the mean local Overall across folds, and on PDS in at least
     3 of 4 folds, or C1 is not submitted.
4. **Submission #2** is then C1 if it passes (3), otherwise C0 as a calibration of the
   public-to-hidden chain. Either way, compare the hidden result to its predeclared public
   expectation, as for v1.

Not in C1: a context main effect `m_hat`. It is a candidate *second* ablation, because
C0 omits the baseline's own mean response by design, but one change at a time. Also not
in C1: STATE, a learned model, or tuning on A/B/C.

## 13. Download and disk accounting (§10)

`df -h ~` before the phase: 1.1 TiB free (1,239,877,017,600 bytes), above the 400 GB rule.

| item | size |
|---|---|
| downloads (6 URL sources + GENCODE + 2 X-Atlas lines) | 251 GB planned; 233 GiB on disk |
| prepared statistics | 0.38 GB |
| C0 prediction / compact / `.vcc` | 4.8 / 0.9 / 3.4 GB |
| peak temporary (compact + pack, upstream estimate) | ≈ 39 + 29 GB |
| free after the phase | 921 GiB |

Nothing was duplicated. Our scPertEval `replogle22k562` is the essential screen, not
GWPS. One downloader defect was found and fixed: curl's internal `--retry` re-used a stale
`-C -` offset and truncated a partial file. Retries now happen in the outer loop, and the
final SHA-256 check guarantees integrity.

## 14. Licensing flags for the human

* **X-Atlas/Orion is CC BY-NC-SA 4.0**: non-commercial, ShareAlike.
* **VCC 2025 H1 and CD4 buckets carry no license file**; their terms are unverified.
* The K562 GWPS license is CC BY 4.0.
* The FAQ allows any data the entrant has rights to, and finalists must disclose data
  and models.

These need a human decision before a submission that uses them.

## 15. Relation to the research track (unchanged)

* **AtlasShift succeeds by adding direct measurements of the same knockdowns.** That is
  the conserved-effect (`beta`) transfer our research endorsed. It is not
  annotation-based unseen prediction, which the research found fails, and it makes no
  claim about `gamma`.
* **The promoter prior is a CRISPRi assay prior**, not a functional-annotation prior.
* **STATE would be pretrained perturbational transfer.** Its success or failure would not
  bear on the STRING/DepMap/pathway finding.

**None of the results in this phase contradict a frozen research conclusion.**

**DO NOT SUBMIT** was honoured. `c0_manifest.json` records `"submitted": false`.
