# C2 predeclaration: expression / DE realisation calibration

Written 2026-09-28, **before any C2 arm was scored**. Its SHA-256 is stored in
`data/provenance/competition_v2/c2_predeclaration_digest.txt`, and
`scripts/competition_v2/run_c2_calibration.py` refuses to run if the file has changed.

## Fixed

* The C1 atlas mean predictor: equal-weight centred fusion of K562 GWPS + VCC 2025 H1 +
  CD4 DE, frozen amplitude log2fc 0.6 / bulk_delta 0.3, clip 3, promoter cap.
  License-clean sources only. **No X-Atlas, no Kaden.**
* Public folds: H1 (150 targets, truth = H1 2025 cells, predictors K562 + CD4) and K562
  (390 targets, truth = K562 GWPS cells, predictors H1 + CD4). Same cells, reference/pool
  split, seeds and gene axes as the C1 folds. CD4 is effect-level only. Its PDS is a
  per-row cosine, which is invariant to any per-target amplitude, so CD4 cannot inform
  E or F and is not rerun.
* Reproduction gate: G0 at amplitude 1.00 must reproduce the frozen C1a fold raw scores.
  PDS must match to ≤ 1e-9, and the other members to ≤ 1e-9 (same seeds).

## Rulers

* **Primary ruler:** the frozen C1 fold anchors (`ANCHOR_mean_response` emitted through
  G0, and `ANCHOR_split_half`), recomputed here and required to equal the C1 values.
  The ruler is the same for every candidate.
* **Sensitivity ruler:** the same oracle mean response emitted through G1c instead of G0.
  It is also fixed across candidates. The hidden scale is affine in each raw member with
  fixed organiser anchors, so the ruler only changes how members are weighted inside
  Overall. A raw-metric improvement is an improvement under any ruler.
* Overall uses MSE clamped at 0, as in C1, because the official V1 MSE of exactly 0
  indicates a clamp. The member-level MSE comparison uses the **unclamped** scaled value,
  which is equivalent to raw MSE.

## Generators (the mean predictor is identical for all of them)

* **G0**: the C1 dual-moment emitter on the 4-cell-averaged template (unchanged).
* **G1c**: 400 real single control cells from the template pool, drawn once per fold
  without replacement and shared by every target (like G0's shared template). Per gene,
  `r_g = p_g / c_g`, where `p` is the target's mean-CPM composition and `c` is the pool's
  mean-CPM composition.
  * Down (`r < 1`): binomial thinning, `Bin(x, r)`.
  * Up (`r > 1`): `x + Poisson((r - 1) · d_i · c_g)`, so a gene can switch on in a cell
    that had zero counts.
  * At `r = 1` the donor cell is returned unchanged. The library depth is preserved in
    expectation.
* **G1b**: G1c with `p` and `c` taken from the pseudobulk (depth-weighted) composition.
* **G2**: `Multinomial(d_i, p_bulk)`, with `d_i` the depths of the same 400 donors.
* **G3** (negative binomial): `Poisson(d_i · p_g · Γ(1/φ_g, φ_g))`, with per-gene `φ_g`
  from the pool controls by the method of moments (`var = μ + φ μ²`, `μ_ig = d_i c_g`,
  floored at 0). **Run only if justified:** G2's null arm must have a median gene-wise
  raw-count variance ratio (generated / real reference, over tested genes) below 0.8 in
  both folds. The number is recorded either way.

No neural generators.

## Null qualification (I)

For each generator, emit a zero effect: `p` equals the pool composition, with no promoter
cap. Score it against the reference with the exact local DE pipeline. A generator
**qualifies** if its median null false-positive DE count per perturbation is ≤ 1 % of
the tested genes in **both** folds (H1: ≤ 107; K562: ≤ 76). G0 is expected to fail
(its C1 null called a median of 3,837 on H1).

## Idealised pseudobulk (D)

No cells are generated.

* Profiles are `log1p(5e4 · p_bulk)` with zero sampling dispersion.
* The DE table is an analytic Welch test at n = 400 against the real reference: the mean
  CPM is `1e6 · p_cpm`, and the variance is the reference's per-gene CPM variance scaled
  by `(m_pred / m_ctrl)²`, followed by BH.
* The log fold change is taken from the means.

This is a mean-response ceiling under real-control noise. It is scored on the primary
ruler.

## Global amplitude (E)

Grid `a ∈ {0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00}`. The multiplier is applied to the
fused effect in both spaces *before* the frozen amplitude and clip (the existing
`shrink` argument).

* Run under G0 and under the selected generator (below).
* Among grid values that pass rule J against C1 (G0, a = 1), choose the one with the
  highest mean primary Overall. If none passes, report the full tradeoff and select
  nothing.
* A nested estimate (select on one fold, score on the other) is reported alongside.

## Generator selection

Among **qualified** generators at a = 1, choose the one with the highest mean primary
Overall that passes rule J. If none passes rule J at a = 1, choose the qualified
generator with the highest mean primary Overall, for use in the amplitude study only.

## Per-target amplitude (F): conditional

Run only if the global study is **consistent**: the selected `a*` ≠ 1.00 improves mean
primary Overall over a = 1 under the selected generator **in both folds**.

If it runs, it is one family with one parameter:
`a_t = a* · (k_t / k̄)^β`, `β ∈ {−0.5, +0.5}`, where

* `k_t` is the number of usable GREEN predictor sources for target `t` in the fold
  (source-side, inference-available);
* `k̄` is its mean over covered targets;
* targets with `k_t = 0` keep `a*`.

It is compared against the global `a*`. No target outcome is used and no regressor is
fitted.

## Rule J (C2 selection, not to be weakened)

A candidate passes only if, as means over the H1 and K562 folds on the primary ruler:

1. Overall is greater than C1's; **and**
2. at least **two** of MSE (unclamped), NMAE, FID, REACH, JAC are greater than C1's;
   **and**
3. PDS is ≥ 95 % of C1's scaled PDS **in each fold**; **and**
4. the generator qualifies under the null test; **and**
5. robustness: Overall is also greater than C1's on the sensitivity ruler.

A candidate that improves DE members but loses more than 5 % of PDS is rejected.

## Build

Only if a candidate passes: emit 3 contexts × 300 targets × 400 cells × 18,533 genes on
the Arc controls, with the same sources and seeds convention. Run all C1 checks and
`vcc prep --dry-run`. Package. **Do not submit.**
