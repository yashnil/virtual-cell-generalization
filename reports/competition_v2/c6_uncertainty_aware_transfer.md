# C6: uncertainty-aware conserved response estimation

**Competition track, 2026-10-02. Nothing was submitted. No bundle was built. C1 is
unchanged.**

## Result

**C6 FAILS. KEEP C1 UNTIL THE FINAL PANEL.**

* Gene-level uncertainty is real for the deep H1 source. It does **not** tell which
  source is right when sources disagree.
* No uncertainty-aware combination improves held-out direction over C1's equal fusion.

## Where things are

| artifact | path |
|---|---|
| predeclaration | `reports/competition_v2/c6_predeclaration.md`, SHA-256 `c0a0ba82…`; digest in `data/provenance/competition_v2/c6_predeclaration_digest.txt` (unchanged) |
| code | `src/virtual_cell/competition_v2/uncertainty.py`; `scripts/competition_v2/run_c6_uncertainty.py`; tests `tests/test_competition_v2_c6.py` (6) |
| outputs | `outputs/competition_v2/c6_uncertainty/`: `{moments_*.npz, sigma2.npz, calibration_deciles.csv, calibration_verdict.json, cross_source_diagnostic.csv, cross_source_diagnostic_raw_sigma.csv, sign_conflicts.csv, candidates_mean_level.csv, identity_protection.csv, branch_decision.json, spectrum.csv, spectrum_verdict.json, folds/*/scores_raw.csv, vcc_scaled.csv, c6_decision.json}` |

## C. Which sources support comparable uncertainty?

| source | estimator | basis | coverage |
|---|---|---|---|
| K562 GWPS | delta-method SE of C1's own log2fc effect, built from per-cell CPM moments | 82k cells (all retained targets), 8,000 controls | 90 % of cells finite |
| H1 2025 | the same estimator | all 490,973 cells of the three original files, each target against its own split's controls | 100 % |
| CD4 DE | `Σ lfcSE² / k²` from the supplied DESeq-style SEs | summary statistics; model-based, **not fabricated** | 79 % |

K562 and H1 use one estimator. CD4's is a different estimator of the same quantity.

## D. Does uncertainty predict reproducibility?

Split-half test: σ² is estimated on half A, and agreement with half B is measured per
σ² decile.

| source | sign agreement, decile 1 → 10 | Pearson, decile 2 → 10 | median \|A−B\|, decile 1 → 10 | verdict |
|---|---|---|---|---|
| **H1** | 0.884 → 0.615 (monotone) | 0.894 → 0.603 | 0.001 → 0.198 | **calibrated** (ρ = −1.0, gap 0.27) |
| **K562** | 0.569 → 0.508 (monotone) | 0.204 → 0.047 | 0.08 → 0.60 | **not calibrated** (gap 0.061 < 0.10) |
| **CD4** (proxy: condition pairs) | 0.558 → 0.571 (flat or reversed) | ~0.25 everywhere | 0.05 → 0.53 | **not calibrated** (ρ = +0.6) |

* **H1:** σ fully predicts reproducibility.
* **K562:** σ predicts the *size* of the error, but almost no direction. Every K562
  decile reproduces its own sign only 51–57 % of the time (split-half reliability
  ≈ 0.09).
* **CD4:** the proxy mixes biology with noise, so it is inconclusive by construction.
* **Consequence (predeclared):** K562 and CD4 carry no gene-level uncertainty into the
  candidates; they get the fold median. U3 shrinks only H1.

**Interpretation note, recorded after calibration and before scoring.** The predeclared
U3 rule ("requires the calibration of the sources it shrinks") was implemented as
"U3 shrinks only calibrated sources". The first code version also shrank uncalibrated
sources with their raw σ. That was corrected before any candidate was scored.

## E. When sources disagree, is the more precise source right? (the key test)

Population: conflicting-sign cells among each target's top-200 truth genes. 95 %
cluster-bootstrap CIs over targets.

| held out | sources | P(lower-σ source correct), predeclared σ | P(lower *raw* σ correct), secondary | conflict cells |
|---|---|---|---|---|
| H1 | K562 vs CD4 | 0.466 (0.440–0.491) | 0.489 (0.467–0.510) | 2,868 |
| K562 | H1 vs CD4 | 0.500 (0.486–0.513) | 0.507 (0.492–0.521) | 8,681 |
| CD4 | K562 vs H1 | 0.514 (0.490–0.539) | 0.515 (0.492–0.540) | 1,993 |

**No.** No fold beats 0.5, so the predeclared gate **stopped the weighting branch**. U1
and U2 are reported at mean level only.

*Exploratory, not a candidate:* picking the source with the larger `|e|/σ` is right
0.553 / 0.535 / 0.519 of the time. That is barely above what C1's equal average already
does implicitly in conflicts, since it follows the larger magnitude: **0.547 / 0.531 /
0.510** (§H). No unexploited signal is left there.

## H. Sign-conflict anatomy (top-200 truth cells)

| held out | agree / conflict / single | C1 | U1 | U2 | U3 | each source alone | median σ², correct vs wrong source |
|---|---|---|---|---|---|---|---|
| H1 | 3,531 / 2,868 / 21,801 | **0.547** | 0.526 | 0.545 | 0.547 | K562 0.534, CD4 0.466 | 0.035 vs 0.007 |
| K562 | 9,788 / 8,681 / 52,531 | **0.531** | 0.529 | 0.532 | 0.526 | H1 0.518, CD4 0.482 | 0.0039 vs 0.0039 |
| CD4 | 2,332 / 1,993 / 61,275 | 0.510 | **0.522** | 0.517 | 0.514 | K562 0.485, H1 0.515 | 0.029 vs 0.035 |

In conflicts, each source alone is close to a coin flip against the held-out truth.
Measurement uncertainty cannot resolve conflicts that neither source gets right.

## F/G/K. Candidates: mean level and identity protection

| held out | candidate | cosine | top-200 sign acc. | cosine to C1 | cells changed (\|Δ\| > 0.1) | CD4 effect PDS |
|---|---|---|---|---|---|---|
| H1 | **C1** | **0.0647** | **0.5640** | 1 | 0 | |
| H1 | U1 (IVW) | 0.0597 | 0.5618 | 0.963 | 11.0 % | |
| H1 | U2 (RE) | 0.0628 | 0.5638 | 0.991 | 2.6 % | |
| H1 | U3 (EB) | 0.0647 | 0.5640 | 1.000 | 0 % (no calibrated predictor) | |
| K562 | **C1** | **0.0451** | **0.5299** | 1 | 0 | |
| K562 | U1 | 0.0444 | 0.5297 | 0.958 | 2.0 % | |
| K562 | U2 | 0.0450 | 0.5299 | 0.998 | 0.1 % | |
| K562 | U3 | 0.0444 | 0.5293 | 0.994 | 0.9 % | |
| CD4 | C1 | 0.0309 | 0.5202 | 1 | 0 | 0.6989 |
| CD4 | U1 | 0.0320 | 0.5206 | 0.909 | 9.0 % | 0.6955 |
| CD4 | U2 | 0.0319 | 0.5204 | 0.981 | 2.9 % | 0.6990 |
| CD4 | U3 | 0.0304 | 0.5203 | 0.994 | 1.3 % | 0.6976 |

No candidate collapses toward control: the norm ratio to C1 is 0.98–1.00.

**VCC, U3 only** (the only candidate past the gates), frozen C1 generator, primary ruler.
C1 reproduces C2's `G0_a1.00` to ≤ 1.1e-16 in both folds.

| fold | arm | PDS | MSE | NMAE | FID | REACH | JAC | Overall |
|---|---|---|---|---|---|---|---|---|
| H1 | C1 | 0.6256 | −0.0356 | −0.1203 | −0.1224 | −0.0120 | 0.0477 | 0.0698 |
| H1 | U3 | 0.6256 | −0.0356 | −0.1203 | −0.1223 | −0.0120 | 0.0477 | 0.0698 |
| K562 | C1 | 0.4366 | −0.0374 | 0.1166 | −0.1578 | 0.1362 | 0.0013 | 0.0888 |
| K562 | U3 | 0.4302 | −0.0524 | 0.1007 | −0.1656 | 0.1357 | 0.0009 | **0.0836** |

In the H1 fold U3 equals C1 by construction. The 1e-4 FID difference comes from a
different floating-point summation order.

**Rule L for U3:** fails criteria 1 (cosine), 2 (sign accuracy), 3 (≥ 2 folds),
4 (Overall ≥ +0.005; actual −0.0026) and 6. **C6 FAILS.**

## I. Is the response matrix low-rank?

| source | components above 2× the noise ceiling | signal energy in them | verdict |
|---|---|---|---|
| H1 (299 targets × 10,780 genes) | **20** | **67 %** | clearly low-rank |
| K562 (384 × 8,246) | 1 | 23 % | noise-dominated beyond one axis |

The predeclared requirement was strong separation in *both* cell sources, so **L1 was
not run.**

## M. Scientific answers

1. **Does estimated measurement uncertainty predict reproducibility?**
   * **H1:** yes, strongly.
   * **K562:** for error size, yes; for direction, barely.
   * **CD4:** untestable within source; the cross-condition proxy shows nothing.
2. **When sources disagree, is the more precise source more often right?**
   * **No.** It is right 0.47–0.51 of the time (CIs include 0.5 or lie below it), with
     raw or gated σ, in every fold.
3. **Is source disagreement mostly noise or biological context?**
   * **H1 vs CD4** (the only pair without K562): mostly **context**. The pooled
     between-source variance is τ² ≈ 0.006, about 4× H1's median σ² of 0.0014, and in
     conflicts the more precise source is not more often right.
   * **Pairs with K562** (K562 vs CD4: τ² = 0; K562 vs H1: τ² ≈ 0.0006): the two cannot
     be separated. K562's measurement noise (σ² ≈ 0.035) swamps any heterogeneity.
4. **Does inverse-variance meta-analysis beat equal averaging?** **No.** U1 lowers
   cosine in H1 and K562 and raises it only in CD4. The mean cosine and sign accuracy
   both fall.
5. **Does random-effects modelling help separate conserved and context-specific
   responses?** **No measurable gain.**
   * τ² absorbs the precision differences, so U2 reverts almost to equal weights: 0.1–3 %
     of cells changed, cosine within ±0.002 of C1.
   * It measures heterogeneity; it does not resolve it.
6. **Is the perturbation-response matrix measurably low-rank?**
   * **In H1, yes** (20 components, 67 % of signal energy).
   * **In K562, no** (1 component).
7. **Can cross-perturbation denoising improve gene direction?** **Not tested.** The
   predeclared spectral gate failed on K562. H1's structure makes it a candidate idea
   for a deep source only, not for the current source mix.

## Leakage audit (§J)

* **Uncertainty inputs:** σ², τ² and the EB τ² use predictor-source cells and summary
  statistics only. H1's σ² is computed from H1's own files, and it is used only in the
  K562 and CD4 folds, where H1 is a predictor.
* **Truth:** it enters only the scoring, the §E gate and the diagnostics. The only
  decision it touched was the predeclared §E gate, which disabled candidates and
  selected none.
* **Calibration:** it uses within-source split halves only.

## Validation

See the research log entry of 2026-10-02.
