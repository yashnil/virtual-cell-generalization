# Reliability-aware transferability / confidence model — v1

**Result: the simple statistic wins.** Raw source agreement predicts transfer
quality well in all four held-out contexts, with monotone risk-coverage and
monotone calibration everywhere including HepG2. Neither monotone calibration
(M1) nor a regularised multivariable model (M2) beat it. **Per the predeclared
decision rule, the final confidence estimator is raw source agreement.**

A second, non-obvious finding: **trustworthiness and error magnitude are
different objectives.** Targeting the least-confident perturbations captures
*less* reproducible error than random selection, because low confidence tracks
low reproducible signal. Ranking by *expected error magnitude* instead captures
roughly twice random.

**The point predictor was not altered.** No gamma model v3, no neural network,
no foundation embeddings, no Arc predictions, no generative model.

Date: 2026-09-20 · runtime 1.0 min
Reproduce: `uv run python scripts/run_transferability_confidence.py`
then `uv run python scripts/plot_transferability_confidence.py`

All ten freeze records verified before and after, including the permanent
pathway-termination freeze. Confidence code is a new module
(`virtual_cell.modelling.transferability`); nothing frozen was modified.

---

## 1. What exact reliability-aware target was used?

With independent target halves `h1 = L + e1`, `h2 = L + e2` and a source-only
prediction `B` (independent of target noise by construction):

```
signal_energy   = <h1, h2>        unbiased for ||L||^2
residual_energy = <h1-B, h2-B>    unbiased for ||L - B||^2
D               = residual_energy / signal_energy
```

Every noise cross-term has zero expectation. Validated on synthetic data with a
known latent signal: `signal_energy` recovers `||L||²` to 5 % at three noise
levels while the naive `||h||²` is inflated by >20 %; `residual_energy` recovers
the true squared error to 5 %; and `D` recovers planted values of 0.25, 0.0625
and 1.0 for predictors at 0.5×, 0.75× and 0×. **`D` is never clipped** — the
`D > 1` regime (worse than predicting zero) is real and was observed.

* **Primary target: reproducible residual energy** — no denominator, no
  instability.
* **Secondary target: `D`**, only where the denominator is estimable.
* **Retained: reliability-normalised directional similarity** (`r/sqrt(rho)`).

## 2. How was near-zero signal handled?

By a stability rule **derived from synthetic simulation before any outer fold
was evaluated**, then frozen (`outputs/transferability_v1/stability_rule.json`,
hashed in `data/provenance/`).

The simulation used the empirical gene-level half-noise sd measured from source
halves (**0.0600**, 6,640 genes) and asked: at what estimated signal energy does
`D` land within 0.25 of its true value in ≥90 % of replicates? The transition is
sharp — within-tolerance rises 0.14 → 0.77 → 1.00 as signal energy goes
0.16 → 1.13 → 3.48 — giving

**min signal energy = 1.9674.**

My first attempt at this returned the grid's own lower bound because the grid
never reached the unstable regime; it was re-derived over a wider range so the
threshold is genuinely estimated rather than an artefact.

In practice the rule excludes very little: the stable fraction is **0.999 / 0.995
/ 0.999 / 0.996** (K562 / RPE1 / HepG2 / Jurkat). Real perturbations sit far
above the threshold, so `D` is usable almost everywhere — but the rule was fixed
in advance regardless of that outcome.

## 3. How strongly does raw source agreement predict transfer failure, by fold?

Spearman between the score and `−D` (higher = better ranking):

| held out | **M0 source agreement** | M1 | learned (M2) | source reliability | source magnitude | random |
|---|---:|---:|---:|---:|---:|---:|
| K562 | **0.554** | 0.531 | 0.531 | 0.504 | −0.492 | 0.021 |
| RPE1 | **0.786** | 0.786 | 0.786 | 0.793 | −0.760 | 0.036 |
| HepG2 | **0.656** | 0.654 | 0.654 | 0.713 | −0.685 | 0.002 |
| Jurkat | **0.556** | 0.538 | 0.544 | 0.521 | −0.504 | −0.007 |

Strong and consistent in every context. Source reliability is comparable (and
slightly better in HepG2 and RPE1), which is expected — the two are related.

## 4. Does simple calibration improve it?

**No.** M1 (isotonic calibration of source agreement, fitted by inner LOCO) is
identical to M0 on ranking in RPE1/HepG2 and slightly *worse* in K562 (0.531 vs
0.554) and Jurkat (0.538 vs 0.556). A monotone transform cannot change rank
order; it only affects the calibrated value, and the inner-fitted transform did
not transfer better than the raw score.

## 5. Does the multivariable model improve it?

**No.** Inner LOCO selected M1 in three folds and
`M2:agreement_plus_magnitude:a=0.1` in Jurkat. On the outer targets:

* M0 ≥ learned on ranking in **all four** contexts (deltas +0.0005 to +0.0229).
* M0 ≥ learned on risk-coverage similarity in **16 of 20** (context, coverage)
  cells; the largest single difference is K562 at 10 % coverage, where M0 is
  **+0.057 better**.

## 6. Which features add information beyond source agreement?

Ridge coefficients (target `D`; negative raises confidence), by mean |coef|:

| feature | mean abs | sign consistent? |
|---|---:|---|
| source_magnitude | 0.215 | **no** |
| basal_similarity_mean | 0.166 | **no** |
| source_reliable_energy | 0.163 | **no** |
| source_reliability | 0.143 | yes |
| source_cells_min | 0.102 | yes |
| source_sd_ratio | 0.092 | no |
| source_sign_consistency | 0.061 | no |
| source_agreement | 0.043 | no |
| target_gene_basal | 0.021 | yes |

Only three of twelve features have a sign-consistent coefficient, and the two
largest are sign-inconsistent — the model is not finding stable structure.
`source_agreement` itself has a small coefficient because it is collinear with
`source_min_agreement` (their fitted coefficients are identical to four decimal
places, a clear collinearity signature) and with `source_reliability`.

**Answer: none of them, reliably.** The ablations confirm it: no feature group
produced a model that beat M0 on the outer targets.

## 7. Are risk-coverage curves monotonic in every context?

**Yes, for M0, in all four contexts and on both metrics.**

Median reliability-normalised similarity, by coverage:

| held out | 100 % | 75 % | 50 % | 25 % | 10 % |
|---|---:|---:|---:|---:|---:|
| K562 | 0.484 | 0.504 | 0.542 | 0.600 | **0.646** |
| RPE1 | 0.611 | 0.634 | 0.663 | 0.699 | **0.718** |
| HepG2 | 0.677 | 0.680 | 0.688 | 0.694 | **0.709** |
| Jurkat | 0.505 | 0.527 | 0.586 | 0.660 | **0.733** |

Random ranking is flat or non-monotone, as it must be.

## 8. How much does reliable risk decrease at 50 %, 25 %, 10 % coverage?

Median `D` (lower is better):

| held out | 100 % | 50 % | 25 % | 10 % |
|---|---:|---:|---:|---:|
| K562 | 1.231 | 0.950 (**−22.8 %**) | 0.783 (**−36.4 %**) | 0.728 (**−40.9 %**) |
| RPE1 | 0.805 | 0.732 (−9.1 %) | 0.694 (−13.8 %) | 0.664 (−17.5 %) |
| HepG2 | 0.895 | 0.688 (**−23.1 %**) | 0.650 (−27.4 %) | 0.612 (**−31.6 %**) |
| Jurkat | 1.011 | 0.827 (−18.3 %) | 0.723 (−28.5 %) | 0.675 (**−33.3 %**) |

Note K562 and Jurkat sit **above D = 1 at full coverage** — on average the
transfer is worse than predicting zero there — and selection brings them below
1. That is the operational value: confidence identifies the subset where the
predictor is actually worth using.

## 9. Experiment prioritisation — a negative result worth stating plainly

Fraction of total reproducible error captured by the **least-confident** budget:

| held out | 5 % | 10 % | 20 % | 30 % | 50 % |
|---|---:|---:|---:|---:|---:|
| K562 (M0) | 0.027 | 0.052 | 0.108 | 0.170 | 0.321 |
| RPE1 (M0) | 0.011 | 0.025 | 0.055 | 0.109 | 0.267 |
| HepG2 (M0) | 0.029 | 0.059 | 0.120 | 0.185 | 0.332 |
| Jurkat (M0) | 0.026 | 0.047 | 0.101 | 0.159 | 0.308 |
| **random** | ~0.05 | ~0.10 | ~0.20 | ~0.30 | ~0.50 |

**Targeting the least-confident perturbations captures roughly half of what
random selection would.** The reason is structural: low source agreement tracks
low *reproducible signal*, and a perturbation with no reproducible response
contributes little absolute error no matter how badly it is predicted.

Ranking instead by **expected error magnitude** (source response magnitude)
captures far more than random:

| held out | 5 % | 10 % | 20 % | 30 % | 50 % |
|---|---:|---:|---:|---:|---:|
| K562 | **0.110** | **0.214** | **0.375** | **0.506** | **0.697** |
| RPE1 | 0.108 | 0.218 | 0.407 | 0.564 | 0.767 |
| HepG2 | 0.121 | 0.240 | 0.414 | 0.523 | 0.693 |
| Jurkat | 0.101 | 0.191 | 0.376 | 0.510 | 0.717 |

So: *"experimenting on the 20 % highest-expected-error perturbations captures
about 38–41 % of the reproducible prediction error, roughly twice random."*

But that same magnitude ranking is a **terrible** trustworthiness score — its
risk-coverage similarity at 10 % coverage is 0.23–0.41, far *below* the
all-perturbation baseline. **The two objectives are close to opposites and must
not be conflated.**

## 10. Does the method work in HepG2?

**Yes — and this is the first phase in which HepG2 has behaved normally.**
Spearman 0.656 (second highest), monotone risk-coverage, monotone calibration,
and the second-largest D reduction at 50 % coverage (−23.1 %). HepG2 was the
known-negative for pathway gamma; confidence estimation is a different question
and it works there.

## 11. Is confidence meaningfully calibrated across contexts?

**Yes, monotonically, in all four.** Mean reliability-normalised quality by
equal-count source-agreement quintile:

| held out | Q1 | Q2 | Q3 | Q4 | Q5 |
|---|---:|---:|---:|---:|---:|
| K562 | 0.294 | 0.402 | 0.452 | 0.540 | 0.606 |
| RPE1 | 0.343 | 0.489 | 0.566 | 0.633 | 0.702 |
| HepG2 | 0.493 | 0.590 | 0.607 | 0.662 | 0.693 |
| Jurkat | 0.303 | 0.410 | 0.472 | 0.541 | 0.665 |

Strictly increasing in every context. **This is not a probabilistic
calibration** and is not claimed as one — the target is a correlation-like
quality, not an event probability. It is an interpretable monotone mapping from
source agreement to expected prediction quality.

## 12. Is a learned confidence model justified?

**No. Source agreement remains the estimator.**

The predeclared rule was: *a learned model is justified only if it improves
meaningfully over raw source agreement on outer-context selective-prediction or
experiment-prioritisation metrics; if M1/M2 do not beat source agreement
consistently, STOP.* M0 is ≥ the learned model on ranking in 4/4 contexts and on
risk-coverage in 16/20 cells. The rule fires: **stop, and adopt the simple
statistic.**

This is a scientifically valid outcome, not a failure. The confidence estimator
is one number computable from the source contexts alone, with no fitting, no
hyperparameters and nothing to leak.

## 13. What does this imply for the project's scientific framing?

The arc of this project has been a sequence of honest negatives converging on a
positive:

* The interaction γ **exists**, is **reproducible**, and is genuinely
  **recoverable** at pathway resolution from biology rather than aggregation.
* But **predicting γ does not improve response prediction** — established twice,
  with the contaminating mechanism identified and removed.
* What *does* work is **knowing when the simple conserved predictor can be
  trusted**. That is reliable, monotone, calibrated, works in every context
  including the one that failed everywhere else, and needs no learned model.

The defensible framing is therefore **not** "we predict context-specific
responses" but **"we predict conserved responses, and we can say in advance,
per perturbation and per unseen context, how much to trust each prediction."**
For the Arc setting that is directly actionable: a confidence score can gate or
weight submitted perturbations without any additional measurement.

The complementary result — that capturing absolute error requires ranking by
expected magnitude, not by confidence — means a practical system should carry
**two** scores for two different questions, and should never use one for the
other's job.

## Limitations

* Four contexts; every per-context number is n = 1.
* `source_agreement` and `source_reliability` perform comparably and are
  correlated; this study does not establish which is primary, only that neither
  needs a learned wrapper.
* The learned model's ridge shows clear collinearity (identical coefficients for
  `source_agreement` and `source_min_agreement`), so its coefficients should not
  be read as feature importances.
* `D > 1` at full coverage in K562 and Jurkat means the point predictor is, on
  average, worse than predicting zero in those contexts — a property of the
  frozen predictor, not of the confidence estimator.
* Calibration is monotone but not probabilistic; no probabilistic claim is made.

## Standing statement

Independent four-context study on public scPertEval data. **Not** a Molina &
Zhang reproduction. Arc contexts A/B/C were not used anywhere. The point
predictor was not altered. No gamma model v3, no neural network, no
foundation-model embeddings, no Arc submission predictions, no single-cell
generative model.
