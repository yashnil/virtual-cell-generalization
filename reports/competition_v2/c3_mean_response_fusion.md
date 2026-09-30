# C3: mean response / source fusion

**Competition track, 2026-09-29. Nothing was submitted. No C3 bundle was built.**
Branch `competition/c1-license-clean`, HEAD `5a28314`, with C2/C3 work uncommitted.

## Result

**C3 FAILS its predeclared rule O. Recommendation: KEEP C1.**

The best candidate is **P**: per-target, source-only reliability weights.

* It raises public Overall by **+0.003** (0.0822 vs 0.0793), in both folds.
* It fails criterion 4: top-200 sign accuracy falls by 0.0001.

No weighting, scaling or consensus scheme over the three license-clean atlases improves
gene directions, and direction is the real bottleneck.

## Where things are

| artifact | path |
|---|---|
| C2 freeze (57 files) | `data/provenance/competition_v2/c2_freeze_sha256.txt` (auto-verified) |
| predeclaration | `reports/competition_v2/c3_predeclaration.md`, SHA-256 `2aa25c72…` (unchanged) |
| fusion module / tests | `src/virtual_cell/competition_v2/fusion_c3.py`; `tests/test_competition_v2_c3.py` (C1a is reproduced **bit for bit**) |
| scripts | `scripts/competition_v2/run_c3_fusion.py` (`--stage mean`, `--stage vcc`), `analyse_c3.py` |
| diagnostic object (§C) | `outputs/competition_v2/c3_fusion/diagnostic_{H1,K562,CD4}.npz`: per-source, fused and true responses (float16), masks, consensus classes, truth-noise fraction |
| tables | `…/c3_fusion/{source_transfer_matrix, magnitude_decomposition, sign_consensus, gene_calibration_*, per_target_source_cosine, source_reliability, candidates_mean_level, vcc_scaled}.csv`, `specs_*.json`, `triggers.json`, `perturbation_dependence.json`, `c3_decision.json` |
| figures | `reports/competition_v2/figures/c3_R{1..5}_*.{png,svg}` + `figures/sources/c3_R*.csv` |

## Setup

* **Folds:** leave-one-atlas-out, so every fold has exactly two predictors.
  * **H1** (predicted from K562 + CD4): VCC and mean level.
  * **K562** (predicted from H1 + CD4): VCC and mean level.
  * **CD4** (predicted from K562 + H1): mean level and effect PDS.
* **Scoring:** the C1 emitter (G0), promoter cap and frozen amplitude are unchanged,
  and the frozen local ruler is used. C1a reproduces C2's `G0_a1.00` exactly in both
  folds (`folds/*/reproduction.json`).
* **Implementation correction before any VCC scoring.**
  * Source-only reliability of the CD4 DE table is about 0: the median over 346 targets
    of `1 − SE²/var(lfc)` is 0.000.
  * So R1, R2 and P gave CD4 weight 0. The first mean-stage run then *removed* the
    CD4-only targets (K562-fold coverage 355 → 146), violating the predeclared property
    that "a target measured by one source keeps that source's response".
  * `combine` now falls back to the equal-weight value wherever the weighted mass is zero
    (test `test_zero_weight_keeps_single_source_responses`).
  * All reported numbers are from the corrected run.

## §D/E. Source quality: transfer matrix (log2fc, centred)

| donor → held out | cosine | top-200 sign acc. | norm ratio ‖p‖/‖t‖ | effect PDS | NMAE (top-200) |
|---|---|---|---|---|---|
| K562 → H1 | **0.085** | **0.571** | 2.50 | 0.887 | 1.15 |
| CD4 → H1 | 0.050 | 0.553 | 2.16 | 0.776 | 1.12 |
| H1 → K562 | **0.082** | **0.551** | 0.40 | 0.874 | 0.97 |
| CD4 → K562 | 0.029 | 0.521 | 0.62 | 0.687 | 0.99 |
| H1 → CD4 | **0.050** | **0.528** | 0.46 | 0.777 | 0.98 |
| K562 → CD4 | 0.029 | 0.518 | 1.62 | 0.687 | 1.03 |

Source-only reliability (median): **H1 0.787** (149 targets, ~1,000 cells each),
**K562 0.088** (384), **CD4 0.000** (346). Truth-noise energy fraction (median):
**H1 0.21**, **K562 0.90**, **CD4 ≥ 0.95**. Only H1 is a clean truth. On K562 and CD4 the
cosine ceiling is about 0.2–0.3.

## Answers

### 1. Why is C1's predicted mean "too small"?

**In effect space it is not.** The fused log2fc response is 2.24× (H1), 0.55× (K562) and
1.31× (CD4) the measured truth norm. Against the noise-corrected truth it is 2.7×, 1.7×
and 5.3×.

The 0.52× / 0.32× in C2 is measured on the emitted pseudobulk, and three things produce it:

1. **The frozen emission amplitude** scales the effect by 0.6 in log2fc and 0.3 in
   bulk_delta.
2. **Fusion averaging.** For two-source targets the fused norm is 0.91× / 0.77× / 0.56×
   the mean single-source norm, because the sources barely agree (median inter-source
   cosine 0.03–0.07).
3. **Different kinds of effect (an inference, not measured here).** C2's truth norm is
   the uncentred real profile, including its generic response and sampling noise. The
   prediction is a *centred*, target-specific effect.

What matters is that the predicted energy is **mostly orthogonal to the truth**. The
cosine is 0.03–0.085, and energy explained is negative (−2.0 H1, −0.28 K562, −0.95 CD4).
The MSE-optimal scale is ≈ cos · ‖t‖ / ‖p‖, which is well below 1. That is why
amplitude-up raised PDS but hurt MSE/NMAE in C2, and why norm-matching shrinks here (§5).

### 2. Which source atlases transfer best?

**The H1 ↔ K562 pair** (cosine 0.082–0.085, sign accuracy 0.55–0.57). **CD4 is the
weakest partner** in every column (0.029–0.050). Transfer is essentially symmetric
within a pair.

### 3. Is donor quality globally stable or context dependent?

**Globally consistent at the atlas level (A), but strongly perturbation-dependent (C).**

* **Atlas level:** the ranking never flips; CD4 is the worst donor wherever it appears.
* **Per target:** the better donor varies. K562 beats CD4 for 68 % of H1 targets, H1
  beats CD4 for 81 % of K562 targets, and H1 beats K562 for 63 % of CD4 targets.
* **Predictability:** the per-target winner is partly predicted by the per-target
  source-reliability difference (Spearman 0.36 on H1, 0.32 on K562, 0.10 on CD4).
* **Context-dependence (B)** could not be tested (§6).

### 4. Does reliability weighting help?

**Slightly, not enough.**

* **Global weights (R1 = R2):** they coincide because CD4's reliability is 0, so any
  exponent gives it zero weight where another source exists. Overall goes from 0.0793
  to 0.0808. PDS rises (H1 0.626 → 0.669) and JAC rises (0.024 → 0.036), but NMAE
  falls by 0.024 (material) and sign accuracy falls by 0.001. Fails criteria 4 and 8.
* **Per-target P:** Overall 0.0822, both folds +0.003, PDS 106 % / 101 %, cosine
  +0.002, no material member loss. Fails **only** criterion 4 (sign accuracy −0.0001).

### 5. Does source-specific scaling help?

**No, it hurts.**

* Nested inner-holdout scales are all below 1: `a_s` from 0.15 to 0.91. Norm-matching to
  a noise-corrected truth means shrinking.
* S1 and S2 collapse count-level PDS (H1 0.626 → 0.36 / 0.32), for Overall 0.038 / 0.036.
* This is the C1b lesson again: shrinking lets the shared pool-vs-reference offset
  dominate count-level PDS.

### 6. Does basal target-context similarity predict the best donor?

**Untestable with the license-clean sources.** C3b was declared infeasible in advance:
the CD4 DE source has no control profile, and CD4 is a predictor or the truth in every
fold. Only one context pair (H1, K562) has basal profiles.

### 7. Does source sign consensus predict correctness?

**Yes, reproducibly.** Top-200 sign accuracy by class:

| class | H1 | K562 | CD4 |
|---|---|---|---|
| agree | 0.670 | 0.587 | 0.570 |
| single source | 0.549 | 0.519 | 0.519 |
| conflict | 0.547 | 0.531 | 0.510 |

Agree beats conflict by 12.2 / 5.6 / 6.0 points, which triggered C3d. However, agreeing
cells are only 12–14 % (H1, K562) and 4 % (CD4) of the top truth genes. Most of the
signal rests on a single source.

### 8. Does consensus weighting improve gene directions?

**No, and by construction it cannot change a sign.** Down-weighting conflicting cells
(×0.75 / ×0.5) leaves sign accuracy identical. Cosine changes by +0.00001 / −0.0002.
Count-level PDS drops, giving Overall 0.0775 / 0.0758 against 0.0793. Consensus is
informative as a **trust signal**, not as an amplitude.

### 9. Is gene-wise bias reproducible?

**No, it is anti-reproducible.** The per-gene signed bias of C1a correlates *negatively*
across folds:

* H1~K562: ρ = −0.50;
* H1~CD4: ρ = −0.78;
* K562~CD4: ρ = −0.35.

An atlas's gene offsets flip sign depending on whether it is the truth or a predictor, so
this is a relative atlas offset, not a transferable gene bias. §K (gene calibration) was
**not built**, as predeclared.

### 10. Which C3 candidate performs best?

**P, per-target source-reliability weights** (primary ruler):

| arm | fold | PDS | MSE | NMAE | FID | REACH | JAC | Overall |
|---|---|---|---|---|---|---|---|---|
| C1a | H1 | 0.626 | −0.036 | −0.120 | −0.122 | −0.012 | 0.048 | 0.070 |
| **P** | H1 | 0.667 | −0.051 | −0.164 | −0.125 | −0.009 | 0.066 | **0.073** |
| C1a | K562 | 0.437 | −0.037 | 0.117 | −0.158 | 0.136 | 0.001 | 0.089 |
| **P** | K562 | 0.441 | −0.030 | 0.125 | −0.158 | 0.141 | 0.002 | **0.092** |
| C1a | mean | 0.531 | −0.037 | −0.002 | −0.140 | 0.062 | 0.025 | 0.0793 |
| **P** | mean | 0.554 | −0.041 | −0.019 | −0.142 | 0.066 | 0.034 | **0.0822** |

Every candidate, as fold means:

| arm | PDS | NMAE | REACH | JAC | Overall | Δ cosine | Δ sign acc. | rule O fails |
|---|---|---|---|---|---|---|---|---|
| R1 = R2 | 0.552 | −0.026 | 0.067 | 0.036 | 0.0808 | +0.0009 | −0.0010 | 4, 8 |
| Q | 0.531 | −0.002 | 0.060 | 0.024 | 0.0787 | −0.00004 | −0.00000 | 1, 3, 4, 5, 6, 8 |
| S1 | 0.342 | −0.016 | 0.056 | −0.023 | 0.0380 | −0.0031 | −0.00003 | 1–6, 8 |
| S2 | 0.329 | −0.008 | 0.054 | −0.028 | 0.0356 | −0.0043 | −0.0010 | 1–6, 8 |
| **P** | 0.554 | −0.019 | 0.066 | 0.034 | **0.0822** | **+0.0021** | −0.0001 | **4** |
| C3d-m | 0.526 | −0.004 | 0.061 | 0.022 | 0.0775 | +0.00001 | 0 | 1, 4, 5, 6, 8 |
| C3d-s | 0.520 | −0.004 | 0.064 | 0.018 | 0.0758 | −0.0002 | 0 | 1, 3, 4, 5, 6 |
| K | not built (§J did not trigger) | | | | | | | |

Q equals C1a to four decimals. The nested inner-holdout donor quality is symmetric
within each pair (K562↔CD4 0.029, H1↔CD4 0.050), so **global source weights are not
identifiable from three atlases**.

### 11. Does it pass the predeclared rule?

**No.** P meets criteria 1, 2, 3, 5, 6, 7 and 8. It fails criterion 4 (sign accuracy),
and the rule is not weakened. Its +0.003 gain is also small against the fold-to-fold
spread.

### 12. Is a submission #3 justified?

**No. KEEP C1.**

## Error sources

* **Magnitude: not the bottleneck.** Measured in the space where we fuse, predictions are
  as large as or larger than the truth. The energy that does not transfer is the problem,
  so scaling up or down trades PDS against MSE/NMAE.
* **Direction: the bottleneck.** Gene-level sign accuracy on the top truth genes is
  52–57 %, and response cosine is 0.03–0.085. Agreeing sources are right 57–67 % of the
  time but cover only 4–14 % of the top genes.
* **Donor selection: a small real lever.** Per-target reliability explains part of which
  donor is better (ρ ≈ 0.3), and exploiting it gives +0.003. Global weights cannot be
  learned from three atlases.
* **Measurement ceiling.** The K562 and CD4 truths are 90–95 % noise energy, so only the
  H1 fold measures direction cleanly. Local gains of this size cannot be resolved there.

## Validation

* `uv run pytest -q`: **602 passed, 1 skipped** (7 new C3 tests).
* `ruff check` and `ruff format --check`: clean.
* `uv build`: OK.
* `verify_c1_state.py --full` (`outputs/competition_v2/c3_fusion/state_full_end.json`):
  **0 failures**.
* C2 freeze manifest: 57 / 57 OK. C1 `.vcc` `fcc4e250…`, and the C1, C2 and C3
  predeclarations: unchanged.
* X-Atlas not used (permission still pending, none recorded); Kaden not used.
