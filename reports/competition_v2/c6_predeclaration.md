# C6 predeclaration: uncertainty-aware conserved response estimation

Written 2026-10-02, **before any uncertainty, calibration, diagnostic or candidate number
was computed**. Its SHA-256 is stored in
`data/provenance/competition_v2/c6_predeclaration_digest.txt`, and
`scripts/competition_v2/run_c6_uncertainty.py` refuses to run if the file has changed.

## Fixed

* **Sources:** the GREEN C1 sources only: K562 GWPS, VCC 2025 H1, CD4 DE.
* **Not used:** X-Atlas, KOLF, Jurkat, VIPerturb, STRING, DepMap, pathways, STATE,
  neural networks, gene embeddings, or any hidden or leaderboard outcome.
* **Unchanged:** the C1 source statistics, effects and centering, the promoter cap, the
  G0 generator, the panel and the output format.
* **The only object changed** is the fused conserved mean `beta_hat[p, g]`, in both
  effect spaces.
* **Folds:** the same as C1–C4, with two predictors in each:
  * H1 (K562 + CD4): VCC and mean level;
  * K562 (H1 + CD4): VCC and mean level;
  * CD4 (K562 + H1): mean level and effect PDS.
* **Baseline:** C1 is the frozen `fusion.fused_effects` (equal weights). Its VCC raw
  members must reproduce C2's `G0_a1.00` to ≤ 1e-12.

## C. Measurement uncertainty `sigma2[s, p, g]` (log2fc space, the space C1 fuses)

**Cell sources (K562, H1)** use a delta-method SE of C1's own log2fc effect, built from
per-cell CPM moments. With the target-cell mean `m_t`, per-cell variance `v_t` and
`n_t` cells, the control `c`, `v_c`, `n_c`, the C1 shrink fraction `f = N / (N + 1e5)`
(where `N` is the target's total counts) and `m' = f·m_t + (1 − f)·c`:

```
sigma2 = [f / ((m'+1) ln2)]² · v_t/n_t  +  [((1−f)/(m'+1) − 1/(c+1)) / ln2]² · v_c/n_c
```

* The evaluation points are the frozen C1 statistics; the variances come from cells.
* **K562 cells:** `fold_cells_K562.npz`, which holds every retained target's cells
  (capped at 1,000) and 8,000 controls.
* **H1 cells:** all cells of the three original files (`adata_Training`,
  `adata_Validation`, `adata_Test`), streamed once, with controls = their
  non-targeting cells.
* The centering-vector variance is neglected (it is a mean over hundreds of targets).

**CD4 (summary statistics only)** uses the supplied DESeq-style `lfcSE`. C1's CD4 effect
is the mean of the centred log2FC over the usable conditions, so
`sigma2 = Σ_c lfcSE²_c / k²`. This is a model-based SE of the same quantity, *not*
fabricated. It is reported as a different estimator from the cell sources.

**Missing SE:** where a source effect exists but its SE is not finite, that
`(s, p, g)` keeps C1's equal treatment. Its weight is the median weight of that source
in that fold.

## D. Calibration test (per source; an estimator that fails is rejected)

* **Cell sources:**
  * Random split halves of each target's cells, with seed `SEED + 61`.
  * Half-level effects use C1's centering vector; `sigma2_A` uses half A's cells only.
  * `(p, g)` cells are binned into deciles of `sigma2_A`.
  * Per decile: sign agreement between halves (the direction-agreement measure),
    Pearson of `e_A` vs `e_B`, and the median `|e_A − e_B|`.
* **CD4 (proxy; no within-condition replicates):** for targets with ≥ 2 usable
  conditions, take each condition's centred lfc as the halves and that condition's SE
  as the uncertainty. This mixes biology with noise and is reported as a proxy.
* **A source's sigma2 is CALIBRATED** if Spearman(decile, sign agreement) ≤ −0.8 **and**
  sign agreement in the lowest-sigma decile minus the highest is ≥ 0.10.
* **An uncalibrated source** gets its fold-median `sigma2` for every gene: no gene-level
  information.

## E. Cross-source diagnostic (the key test)

* **Population:** the cells where both predictors measured the gene and disagree in
  sign, restricted to each target's top-200 `|truth|` genes (the sign-accuracy
  population).
* **Statistic:** P(the source with the lower `sigma2` has the held-out sign), with a
  95 % CI from 2,000 cluster-bootstrap resamples over targets.
* **Also reported:** each source's own conflict accuracy, and the result restricted to
  the H1 fold, the only clean truth.
* **The weighting branch (U1, U2) proceeds only if** the CI lower bound exceeds 0.5 in
  ≥ 2 of the 3 folds. Otherwise U1 and U2 are **stopped**: they are reported at mean
  level only and not VCC-scored.

## F. Candidates (no tuned functions)

* **U1, fixed-effect inverse variance:** `w[s,p,g] = 1 / sigma2[s,p,g]`. It is a
  per-cell weighted mean over the measuring sources, with CD4 folded in the same way.
* **U2, random effects:** `w = 1 / (sigma2 + tau2_g)`. `tau2_g` is the DerSimonian–Laird
  moment estimate **pooled per gene over all perturbations** measured by both
  predictors of the fold:
  ```
  tau2_g = max(0, (Σ_p Q_pg − Σ_p (k_pg − 1)) / Σ_p (S1_pg − S2_pg / S1_pg))
  ```
  with `w = 1/sigma2`, `S1 = Σ w` and `S2 = Σ w²`. This uses source-side data only.
* **U3, empirical Bayes:**
  * Each source's effect is shrunk toward zero: `e* = e · τ²_{s,g} / (τ²_{s,g} + sigma2)`.
  * `τ²_{s,g} = max(0, mean_p e² − mean_p sigma2)` over the source's usable retained
    targets.
  * The same multiplier is applied to the source's bulk_delta effect; CD4 is shrunk in
    log2fc before its bulk conversion.
  * Then C1 equal fusion is applied.
* The weights for both spaces come from the log2fc uncertainty.
* **U3 requires** the calibration (§D) of the sources it shrinks. It does not depend
  on §E.

## G. Identity protection (reported for every candidate)

* the response norm ratio;
* PDS (VCC folds and CD4 effect PDS);
* the per-target cosine to C1;
* the fraction of measured `(p, g)` cells changed materially from C1, i.e.
  `|Δ| > 0.1` log2.

## H. Sign-conflict analysis

Every top-200 truth cell is classified as agree, conflict or single-source. Within the
conflict cells, report the accuracy of C1, U1, U2, U3 and each source, the median
`sigma2` of correct vs wrong source calls, and the counts.

## I. Optional low-rank diagnostic (only after U1–U3)

* **Spectrum:** for each cell source, compute the spectrum of the centred effect matrix
  (usable targets × genes with control CPM ≥ 5). Compare it with the noise spectrum of
  the split-half matrix `(A − B)/2`.
* **Strong separation** means `r* = #{i : s_i² > 2 · max noise s²}` is ≥ 3 in both cell
  sources, **and** those components hold ≥ 20 % of the signal energy above noise.
* **Only then is L1 tested:**
  * a basis is learned from the predictor sources' targets **outside** the fold panel;
  * each source response is projected onto it, followed by the best of U1–U3;
  * the rank `r ∈ {5, 10, 20, 50}` is chosen by inner cross-source reconstruction
    cosine on those training targets.
* No autoencoder, no neural model.

## J. Leakage

* `sigma2`, `tau2` and `τ²` use predictor-source data only.
* A held-out atlas never contributes to its own fold's predictors, uncertainties or
  hyperparameters.
* The truth enters only the scoring and the §E / §H diagnostics, which select nothing
  except the predeclared §E gate.

## L. Pass rule (not to be weakened)

C6 succeeds only if the best candidate, chosen as the one with the highest mean public
Overall among those VCC-scored, satisfies **all** of the following:

1. The mean-response cosine (mean of the 3 fold means) improves over C1.
2. Top-200 sign accuracy (mean of the 3 fold means) improves over C1.
3. Both cosine and sign accuracy improve in ≥ 2 of the 3 folds.
4. The mean public Overall (H1, K562) improves by **≥ +0.005**. A +0.001-scale gain
   does not justify a new model.
5. PDS is ≥ 98 % of C1 in H1, K562 and CD4 (effect).
6. No fold contributes more than 75 % of the summed Overall gain.
7. No member (MSE, NMAE, FID, REACH, JAC; fold mean, primary ruler) drops by more than
   0.02.

**If C6 passes:** build and package the candidate, report the SHA, then stop for a human
decision. **If it fails:** keep C1; no C7.
