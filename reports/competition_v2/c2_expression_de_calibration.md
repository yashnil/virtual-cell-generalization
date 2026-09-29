# C2: expression / DE realisation calibration

**Competition track, 2026-09-28/29. Nothing was submitted. No C2 bundle was built.**
Branch `competition/c1-license-clean`, HEAD `5a28314`, with C2 work uncommitted.

## Result

**C2 FAILS its predeclared selection rule. Recommendation: KEEP C1.**

No change to the cell generator or the effect amplitude improves the public folds in a
way that survives rule J:

* The only arms that raise public Overall are **G0 (the C1 emitter) at a = 1.25 or 1.5**.
  They fail criterion 4, because the null test disqualifies the C1 generator. Even
  setting that aside, their gain is confined to the K562 fold, whose DE truth is nearly
  empty; H1 is flat.
* Realistic generators (G1*) fix the DE over-calling and raise PDS. They then lose more
  on FID, JAC and MSE than they gain, because those members, as scored, *reward* the C1
  emitter's over-calling and its low-noise pseudobulk.

## Where things are

| artifact | path |
|---|---|
| C1 official result (frozen) | `reports/competition_v2/c1_official_result.md`, `outputs/…/c1_license_clean/official_score_parsed.json` |
| predeclaration + 2 amendments | `reports/competition_v2/c2_predeclaration{,_amendment_1,_amendment_2}.md`; SHA-256 in `data/provenance/competition_v2/c2_predeclaration_digest.txt` (all unchanged at the end) |
| generators | `src/virtual_cell/competition_v2/realisation.py` (G1/G2/G3, ideal DE, cell diagnostics); tests `tests/test_competition_v2_c2.py` |
| fold runner / decisions | `scripts/competition_v2/run_c2_calibration.py`, `analyse_c2.py` |
| candidate builder (not run: rule failed) | `scripts/competition_v2/build_c2_candidate.py` |
| outputs | `outputs/competition_v2/c2_calibration/{folds/<fold>/phase*_*.csv, generator_decision.json, amplitude_decision.json, c2_decision.json, scores_scaled_all.csv, scores_scaled_fold_mean.csv, null_qualification.csv}` |
| figures | `reports/competition_v2/figures/c2_{A_official_scorecard, B_mean_vs_generator, C_amplitude, D_null_de}.{png,svg}`, sources in `figures/sources/` |

## Method, briefly

* **Folds:** the C1 public leave-one-atlas-out folds, with the same cells, splits, seeds
  and gene axes. H1 has 150 targets (predicted from K562 + CD4). K562 has 390 targets
  (predicted from H1 + CD4). CD4 is effect-level only, and a per-row cosine is invariant
  to amplitude, so it cannot inform this phase.
* **Primary ruler:** the frozen C1 local anchors, recomputed and matching to 1e-16. The
  hidden scale is affine in each raw member, so any raw improvement is a hidden
  improvement. The ruler only weights members inside Overall.
* **Sensitivity ruler:** the oracle mean response emitted through G1c.
* **Evidence the primary ruler is right for FID:** back-solving V1's official FID
  (raw 0.407 → −0.349) gives an official FID baseline anchor of about 0.54. That sits in
  the over-calling regime of our G0 anchor (0.556–0.572), not the realistic-cell regime
  (0.20–0.29).
* **Reproduction gate:** the anchors and the null reproduce to 1e-16. C1a reproduces to
  5.7e-7, because the C1 code was edited after its fold scores were written (amendment 2).
  The C1 baseline below is this run's `G0_a1.00`.
* **Process notes:** two runs were restarted. One stopped at the reproduction gate; the
  other was swapping (13 GB of swap at 9 workers on H1) and was rerun with 5 workers.
  Their logs are kept as `phase1_attempt*.log`. G0 at a ≥ 1.75 cannot be emitted: the
  frozen dual-moment fit misses its 1e-3 tolerance. Those arms are recorded as
  infeasible in `phase2_infeasible.json`, not tuned around.

C1 baseline on the primary ruler (fold mean):

| PDS | MSE | NMAE | FID | REACH | JAC | Overall |
|---|---|---|---|---|---|---|
| 0.531 | −0.037 | −0.002 | −0.140 | +0.062 | +0.025 | **+0.079** |

## C. Public anatomy: is the hidden pattern reproduced?

Hidden C1: PDS 0.602, with every other member ≤ 0.12.

* **H1 reproduces the hidden shape.** It has PDS 0.626; MSE < 0 (clamped); FID −0.12;
  REACH −0.01; JAC +0.05; NMAE −0.12.
* **K562** has lower PDS (0.437) and higher REACH and NMAE (+0.136, +0.117). Its truth
  calls a median of **3** DE genes per target, so its FID, JAC and REACH mostly measure
  calls against an almost empty reference.

Raw ingredients, C1 at a = 1 (H1 / K562):

| ingredient | H1 | K562 |
|---|---|---|
| predicted DE genes / perturbation (median) | 5,712 | 1,378 |
| true DE genes / perturbation | 1,296 | 3 |
| directional precision (pooled) | 0.545 | 0.547 |
| DE overlap with truth (median genes) | 751 | 1 |
| predicted / true pseudobulk effect norm | **0.52** | **0.32** |
| library-size bias (median) | +2.5 % | +1.0 % |
| genes-detected bias | **+10.7 %** | **+19.9 %** |
| gene-wise variance, generated / real | **0.50** | **0.64** |
| log-CPM variance ratio | 0.37 | 0.79 |
| cell–cell cosine distance, generated vs real | 0.036 vs 0.093 | 0.231 vs 0.338 |
| pseudobulk (jackknife) dispersion ratio | 1.01 | 0.28 |

## D. Mean response vs generator

"IDEAL" scores the predicted mean directly, with no cells.

* Profiles are `log1p(5e4 · p_bulk)` with zero sampling noise.
* DE is an analytic Welch test at n = 400 under the real reference's per-gene noise.

Primary ruler:

| arm | fold | PDS | MSE | NMAE | FID | REACH | JAC | Overall | pred. DE |
|---|---|---|---|---|---|---|---|---|---|
| IDEAL (mean ceiling) | H1 | 0.625 | −0.083 | −0.121 | −0.588 | +0.018 | −0.181 | −0.041 | 1,399 |
| G0 = C1 cells | H1 | 0.626 | −0.036 | −0.120 | −0.122 | −0.012 | +0.048 | +0.070 | 5,712 |
| IDEAL | K562 | 0.435 | −0.217 | +0.112 | −0.874 | +0.275 | +0.125 | +0.012 | 17 |
| G0 = C1 cells | K562 | 0.437 | −0.037 | +0.117 | −0.158 | +0.136 | +0.001 | +0.089 | 1,378 |

### 1. Which C1 weaknesses come from the mean prediction?

* **PDS is entirely mean-driven.** It is identical with and without cells
  (0.625 / 0.626; 0.435 / 0.437).
* **Magnitude is under-predicted.** The predicted pseudobulk effect is only 52 % (H1)
  and 32 % (K562) of the true effect norm.
* **MSE and NMAE are already weak at pseudobulk level**: IDEAL MSE −0.08 / −0.22, and
  H1 NMAE −0.12, identical with cells. Expression accuracy and log-FC accuracy are
  **magnitude/direction errors of the mean**, not generator artifacts.
* **Directional precision is ≈ 0.55–0.57 even for the exact mean.** The direction of
  the transferred signature is only modestly right gene by gene. This caps FID and JAC
  whatever generator is used.

### 2. Which come from the generator?

* **The DE call volume.** The mean alone implies about 1,400 calls on H1 (real 1,296)
  and 17 on K562 (real 3). The G0 cells produce 5,712 and 1,378.
* **The heterogeneity deficit**: half the gene-wise variance, +11–20 % genes detected,
  and cells about 2.6× too similar to one another.
* **But the generator artifacts are score-positive on the members that matter most:**
  * **FID** = correct calls / max(n_pred, n_real). At ≈ 54 % precision, calling more
    genes maximises it. G0 lifts FID from −0.59 to −0.12 (H1) and from −0.87 to −0.16
    (K562) over the realistic-count ceiling.
  * **MSE** credits `min(C_pred, C_real)` of sampling noise. G0's cells *look* noisy per
    cell, but its pseudobulk is a 1,600-donor average, so it collects credit for noise
    it barely carries. It beats the noiseless ideal (−0.036 vs −0.083).
  * **JAC** on H1 also favours G0 (+0.05 vs −0.18), since more calls mean more overlap
    with a 1,300-gene truth.
  * **Only K562 JAC and REACH penalise the over-calling.** That is the almost-empty-truth
    regime.

## E. Global amplitude (grid 0.50–2.00, before generation)

Fold mean, primary ruler (`figures/c2_C_amplitude.png`):

| arm | PDS | MSE | NMAE | FID | REACH | JAC | Overall | rule J fails |
|---|---|---|---|---|---|---|---|---|
| G0 a 0.50 | 0.408 | −0.062 | −0.011 | −0.128 | 0.049 | −0.002 | 0.052 | 1, 2, 3, 4, 5 |
| G0 a 0.75 | 0.478 | −0.043 | +0.001 | −0.138 | 0.046 | 0.013 | 0.067 | 1, 3, 4, 5 |
| **G0 a 1.00 (C1)** | 0.531 | −0.036 | −0.002 | −0.140 | 0.062 | 0.025 | **0.079** | — |
| G0 a 1.25 | 0.568 | −0.036 | −0.017 | −0.144 | 0.081 | 0.033 | 0.087 | **4 only** |
| G0 a 1.50 | 0.596 | −0.041 | −0.045 | −0.150 | 0.086 | 0.040 | **0.090** | **4 only** |
| G0 a 1.75 / 2.00 | infeasible: the emitter's moment fit fails on H1 | | | | | | | |
| G1c a 0.50 | 0.498 | −0.768 | −0.034 | −1.013 | 0.078 | −0.143 | −0.102 | 1, 2, 3, 5 |
| G1c a 1.00 | 0.604 | −0.841 | −0.026 | −0.425 | 0.107 | −0.056 | 0.034 | 1, 2, 5 |
| G1c a 1.50 | 0.644 | −1.142 | −0.071 | −0.286 | 0.109 | −0.012 | 0.064 | 1, 2, 5 |
| G1c a 1.75 (a\*) | 0.656 | −1.597 | −0.107 | −0.243 | 0.115 | +0.002 | 0.070 | 1, 2, 5 |
| G1c a 2.00 | 0.661 | −2.128 | −0.157 | −0.227 | 0.113 | +0.011 | 0.067 | 1, 2, 5 |

### 3. Does a global amplitude adjustment help?

**Not in a way the rule accepts, and not robustly.**

* Under G0, a = 1.5 raises the fold-mean Overall to 0.090 (+0.011). PDS goes to 0.596,
  REACH to 0.086 and JAC to 0.040, while MSE, NMAE and FID get slightly worse.
* **Per fold the gain is K562-only.** H1 Overall is 0.070 → 0.071 → 0.068 at
  a = 1 / 1.25 / 1.5. K562 goes 0.089 → 0.103 → 0.112, driven by REACH on the
  3-DE-gene truth.
* On the realistic-count ceiling (IDEAL), amplitude helps everywhere: at a = 2, FID goes
  from −0.59 to −0.28 on H1. This confirms that the mean is under-sized (§1). Under G0,
  though, FID is saturated by the over-calling, so amplitude cannot buy FID.
* Nested G1c selection (pick on one fold, score on the other) picks a = 1.75 both ways:
  * K562 0.121 vs C1 0.089;
  * **H1 0.020 vs C1 0.070**.

### 4. Does per-target amplitude improve further?

**No.** It ran because G1c's a\* = 1.75 beats a = 1 in both folds (the predeclared
trigger). With `a_t = 1.75 · (k_t / k̄)^β`, where `k_t` is the number of GREEN predictor
sources:

* β = −0.5 gives Overall 0.061;
* β = +0.5 gives 0.067;
* the global a\* gives 0.070.

Both are worse.

## I. Null DE qualification (zero effect vs the real reference)

| generator | H1 median FP DE (of 10,776) | K562 (of 7,681) | fraction up | qualifies (≤ 1 %) |
|---|---|---|---|---|
| **G0 (C1)** | **3,837 (35.6 %)** | **1,053 (13.7 %)** | 0.75 / 0.87 | **no** |
| G1c / G1b (shared real donors) | 0 | 1 | — | yes |
| G1ci (per-target real donors) | 0 | 0 | — | yes |
| G2 (multinomial) | 200 | 59 | 0.73 / 0.82 | no (H1) |
| G3 (NB) | not run: G2 null variance ratio 0.87 / 0.93 ≥ 0.8 gate | | | |

### 5. Does the current generator overproduce DE under the null?

**Yes, massively.** From zero effect it manufactures a median of 3,837 DE genes (H1) and
1,053 (K562), mostly "up". The cause is its 4-cell-averaged template: fewer zeros, so
higher Wilcoxon ranks. The null Jaccard is 0 in both folds. Real-donor generators call
0–1.

## G/H. Generator realism

### 6. Which generator best matches real cell heterogeneity?

**G1ci.** It uses real control donors, drawn per target, with binomial thinning and
depth-proportional Poisson additions.

| ratio / value vs real held-out cells | H1 | K562 |
|---|---|---|
| gene-wise variance ratio | 0.997 | 1.010 |
| Fano ratio | 0.993 | 1.011 |
| zero fraction, generated vs real | 0.516 vs 0.517 | 0.622 vs 0.622 |
| cell–cell cosine distance | 0.090 vs 0.093 | 0.336 vs 0.338 |
| genes-detected bias | +0.3 % | +0.1 % |

G1c is as realistic per cell. Its shared donors, however, make the MSE noise credit
collapse (the ρ mechanism, amendment 1). G2 is under-dispersed (variance 0.87) and fails
the null. G0 is the least realistic generator. **Realism does not pay under this
scorer.**

### 7. Which candidate improves hidden-relevant metrics publicly?

Only G0 at a = 1.25 / 1.5, and only through PDS, REACH and JAC, and only in K562. Every
realistic-generator arm lowers public Overall. The best, G1c at a = 1.75, reaches 0.070
against 0.079.

### 8. How much C1 PDS is retained?

Every realistic arm raises PDS:

| arm | H1 | K562 |
|---|---|---|
| G1ci a 1 | 115 % | 115 % |
| G1c a 1.75 | 119 % | 130 % |
| G0 a 1.5 | 112 % | 113 % |

PDS was never the limiting criterion. The only arms that fail it are G1b, G2 and a = 0.5.

## J. Rule J verdict

| arm (best of family) | 1 Overall ↑ | 2 ≥ 2 members ↑ | 3 PDS ≥ 95 % | 4 null-qualified | 5 sensitivity ↑ | pass |
|---|---|---|---|---|---|---|
| G0 a 1.50 | ✓ (0.090) | ✓ (REACH, JAC) | ✓ | **✗** | ✓ | **no** |
| G1c a 1.75 | ✗ (0.070) | ✗ | ✓ | ✓ | ✗ | no |
| G1ci a 1.00 | ✗ (0.021) | ✗ | ✓ | ✓ | ✗ | no |
| G1c a 1.75, β ±0.5 | ✗ | ✗ | ✓ | ✓ | ✗ | no |

### 9. Is C2 justified?

**No.** No arm passes, and the rule is not weakened. The one near-miss (G0 a = 1.5) fails
qualification by construction, and its public gain rests on the fold whose DE truth is
uninformative.

### 10. Exact candidate for submission #3

**None from C2. KEEP C1** (`c1_license_clean_val.vcc`, `fcc4e250…`, already scored). No
C2 bundle was generated or packaged.

## Interpretation (binding for the next phase)

* The remaining hidden deficit is **not a cell-generation problem that a better generator
  fixes**. Under this scorer:
  * FID rewards over-calling at ≈ 55 % precision;
  * MSE credits apparent per-cell noise.

  So the C1 emitter's artifacts are, on balance, locally optimal.
* **The real bottleneck is in the mean.** The transferred signature is gene-by-gene only
  ≈ 55 % directionally right, and about 2–3× too small. Scaling it up amplifies wrong
  components as much as right ones, so MSE and NMAE worsen while PDS rises.
* The **hidden-vs-public gap on FID** (hidden −0.02 vs local −0.12 to −0.16) is
  consistent with a hidden FID anchor near 0.54 and a hidden truth with more DE genes.
  Local FID is therefore conservative for G0.
* A human-level decision, **not taken here**, is whether to predeclare a separate
  "G0 amplitude only" candidate without the null-qualification criterion. The public
  evidence for it is +0.011 mean Overall, entirely from K562. I do not recommend spending
  a submission on it.

## N. Validation

* `uv run pytest -q`: **595 passed, 1 skipped** (8 new C2 tests).
* `ruff check`: clean. `ruff format --check`: clean.
* `uv build`: sdist + wheel OK.
* `verify_c1_state.py --full` (`outputs/competition_v2/c2_calibration/state_full_end.json`):
  **0 failures**. This covers 16 freeze manifests / 287 digests, 7 raw-checksum files /
  20 digests, the C0 and V1 bundles, and the C1 `.vcc` SHA unchanged (`fcc4e250…`).
* C1 predeclaration (`fd4efa6f…`), C2 predeclaration (`dcb0df82…`) and amendments 1
  (`02c360b9…`) and 2 (`2c32a248…`): unchanged.
* X-Atlas: not used (still pending permission). Kaden: not used.
