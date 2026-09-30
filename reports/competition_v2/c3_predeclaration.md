# C3 predeclaration: mean response / source fusion

Written 2026-09-29, **before any C3 diagnostic or candidate was computed**. Its SHA-256 is
appended to `data/provenance/competition_v2/c3_predeclaration_digest.txt`, and
`scripts/competition_v2/run_c3_fusion.py` refuses to run if the file has changed.

## Fixed (hard rules)

* **Unchanged:** the C1 generator and count emission (G0 dual-moment emitter, frozen
  amplitude 0.6 / 0.3, clip 3); the license-clean sources (K562 GWPS, VCC 2025 H1,
  CD4 DE); the promoter cap; the panel; the output format.
* **Not used:** X-Atlas, Kaden, STRING/DepMap, gamma/pathway, STATE, neural models, or
  any hidden or leaderboard outcome.
* The only object changed is the fused mean effect `beta_hat[p, g]`, in both effect
  spaces.

## Folds (leave-one-atlas-out)

With three GREEN atlases, every fold has exactly two predictor sources:

| fold | truth | predictors | scored with |
|---|---|---|---|
| H1 | H1 2025 | K562 + CD4 | mean-level + full VCC (G0 cells) |
| K562 | K562 GWPS | H1 + CD4 | mean-level + full VCC |
| CD4 | CD4 DE | K562 + H1 | mean-level + effect PDS |

Cells, splits, seeds and anchors are the C1/C2 ones. The C1 baseline is this run's C1a,
which must reproduce the C2 run's `G0_a1.00` raw members to ≤ 1e-12. With only two
predictors per fold, "majority" sign consensus cannot occur. Consensus classes are
**agree** (both sources, same sign), **conflict** and **single-source**.

## Mean-level definitions (log2fc, centred, as C1 fuses)

* **Truth:** the held-out atlas's own centred log2fc effect through the same code: H1
  and K562 from their statistics, CD4 via `cd4_effect`.
* **Genes:** the fold genes minus the panel target genes, per target restricted to genes
  where both truth and prediction are defined.
* **Cosine and Pearson:** per target, averaged over targets with a prediction.
* **Norm ratio:** median of `‖p‖/‖t‖`. It is also reported against a noise-corrected
  truth norm:
  * for H1 and K562, a split-half of the truth cells (`‖A − B‖²/4` is the noise
    energy);
  * for CD4, `Σ lfcSE² / k`.
* **Sign accuracy:** per target, over the 200 genes with the largest `|t|` where
  `p ≠ 0`, the fraction with `sign(p) = sign(t)`, averaged over targets.
* **Energy explained:** `1 − Σ‖t − p‖² / Σ‖t‖²`, pooled.

## Candidate family (capacity ladder)

All candidates redistribute weight among the available direct sources. A target
measured by one source keeps that source's response. Nothing shrinks the fused response
toward zero except the conditional consensus factor (C3d).

| rung | id | definition |
|---|---|---|
| 1 | **C1a** | equal weights, no scaling (baseline) |
| 2 | **R1** | `w_s ∝ rel_s`, where `rel_s` is source-only reliability (below) |
| 2 | **R2** | `w_s ∝ rel_s^0.5` |
| 2 | **Q (C3a)** | `w_s ∝ max(q_s, 0)`, where `q_s` is the mean per-target transfer cosine of source `s` to the **inner** held-out atlas (the third atlas, never the outer one). Equal weights if every `q ≤ 0` |
| 3 | **S1** | per-source scale `a_s` with equal weights. `a_s` = median over targets of (noise-corrected inner-truth norm / `‖e_s‖`) on common genes, from the inner held-out atlas only. The same `a_s` is applied in both spaces |
| 3 | **S2** | `a_s` scales with R1 weights |
| 4 | **C3b** | **declared infeasible.** Basal-control similarity needs a control profile for every source and target context, and the CD4 DE source has none. CD4 is a predictor or the truth in every fold |
| 5 | **P (C3c)** | per-target weights `w_{s,t} ∝ rel_{s,t}`, per-target source-only reliability. Falls back to `rel_s` where a target's reliability is not computable. Linear, no fitted parameters |
| 5 | **C3d-m / C3d-s** | **conditional on §L.** The fused response times `f`, where `f = 1` for agree and single-source genes and `f = 0.75` (moderate) or `0.5` (strong) for conflicting `(target, gene)` cells. The sign comes from the log2fc components. Applied in both spaces |
| 6 | **K** | **conditional on §J.** `beta_cal = beta + λ_g · b_g`, where `b_g` is the per-gene signed bias from the **training** atlases only, and `λ_g = n_g / (n_g + 100)` with `n_g` the number of training observations of gene `g` |

### Reliability (source data only)

* **K562 and H1** (cells): per target, the Spearman–Brown split-half Pearson,
  `2r/(1 + r)`. The halves are fixed-seed random splits. The measure uses centred
  log2fc on genes with control mean CPM ≥ 5. A target needs ≥ 40 cells. `rel_s` is
  the median over targets.
* **CD4** (DE table): per usable (target, condition),
  `1 − mean(lfcSE²) / var(lfc)`, clipped to [0, 1]. The per-target value is the mean
  over its usable conditions. `rel_s` is the median.

### Conditional triggers

* **§L → C3d:** the top-200 sign accuracy of *agree* genes must exceed that of
  *conflict* genes by ≥ 5 percentage points in **every** fold (all three folds have
  both classes).
* **§J → K:** the per-gene signed bias of C1a, `mean_t (t − p)`, must have Spearman
  ρ ≥ 0.3 across **every** pair of folds on common genes (the genes measured in both
  folds).

### Ladder stop rule

Mean-level metrics are computed for every candidate. VCC scoring proceeds rung by rung,
2 → 3 → 5 → 6. The first rung containing a candidate that passes rule O is the selected
rung. Within it, the passing candidate with the highest mean Overall is chosen. Higher
rungs are **not** VCC-scored once a rung passes.

## Rule O (C3 may replace C1 only if all hold; not to be weakened)

1. The mean public Overall (H1, K562, primary ruler, MSE clamped) is greater than C1's.
2. PDS is ≥ 95 % of C1's in H1 and K562 (VCC) **and** CD4 (effect PDS).
3. The mean-response cosine, averaged over the three fold means, is greater than C1's.
4. Sign accuracy, averaged over the three fold means, is greater than C1's.
5. Overall improves in **both** H1 and K562.
6. Neither fold accounts for more than 75 % of the summed Overall gain.
7. No BLOCKED or UNKNOWN source is used.
8. At least one of MSE (unclamped), NMAE, FID or REACH improves (fold mean), and none of
   the four drops by more than 0.02 scaled.

## Build

Only if a candidate passes: emit the Arc bundle with the unchanged C1 emitter, run all
C1 checks and `vcc prep --dry-run`, and package. **Do not submit.**
