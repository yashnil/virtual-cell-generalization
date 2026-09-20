# Transferability foundations — v1

**Diagnostic only. No transferability model, context encoder, GNN, foundation
model, generative model or Arc submission model was built.** Independent
four-context study on standardized public scPertEval data — **not** a Molina &
Zhang reproduction.

Date: 2026-09-20 · runtime 0.8 min · peak RSS 7.7 GB
Reproduce: `uv run python scripts/run_transferability_foundations.py`
then `uv run python scripts/plot_transferability_foundations.py`

All prior freeze records maintained and re-verified: protocol 4/4, canonical v1
12/12, decomposition phase 21/21, zero-shot v1 14/14. Nothing earlier was
modified.

Every new estimator carries a **leakage test** that replaces the held-out
context's entire response matrix with noise and requires the estimate to be
bit-identical (`tests/test_foundations.py`).

---

## A. Basal template / scale recoverability

### A.1 Descriptive — basal expression does *not* encode the response template

`alpha_c` is the context-level response template (evaluation-only, from the
four-context decomposition); the basal deviation is
`control_mean[c] - mean_{all} control_mean`, computable at inference.

| held out | r(alpha, basal dev) | cosine | Spearman(\|alpha\|,\|dev\|) | fraction of ‖alpha‖² inside the source basal span | fitted global scalar k |
|---|---:|---:|---:|---:|---:|
| K562 | −0.056 | −0.036 | +0.290 | 0.067 | −0.0071 |
| RPE1 | −0.046 | −0.020 | +0.265 | 0.082 | −0.0078 |
| HepG2 | −0.263 | −0.272 | +0.336 | **0.0015** | **+0.0025** |
| Jurkat | +0.065 | +0.068 | +0.333 | 0.011 | −0.0134 |

Three independent signals all say the same thing:

* **Direction is absent.** Gene-wise correlation is ≈0 and flips sign across
  folds; one fold is meaningfully *negative*.
* **The subspace is wrong.** Only **0.15 %–8.2 %** of `alpha`'s energy lies in
  the span of the source basal deviations. There is essentially nothing for a
  linear map from basal space to reach.
* **The scalar is not stable.** The fitted global coefficient changes sign
  across folds, so even a one-parameter mapping does not generalise.

The only positive signal is a weak magnitude association (Spearman +0.27 to
+0.34 between `|alpha|` and `|basal deviation|`): genes that differ a lot
basally also tend to have larger template components, but **without a
consistent direction that is not usable for prediction.**

### A.2 LOCO performance — response energy explained

The key metric. `zero` is the uncorrected source mean; `scale_only` shrinks the
perturbation-specific component by a scalar fitted by **inner leave-one-source-out**
on source contexts only; `oracle_*` rows use the true `alpha` and are ceilings,
not methods.

| estimator | K562 | RPE1 | HepG2 | Jurkat | available? |
|---|---:|---:|---:|---:|:--:|
| zero (source mean) | −0.467 | +0.098 | −0.065 | −0.184 | yes |
| direct_basal | −17.19 | −6.03 | −10.51 | −10.62 | yes |
| global_scalar | −0.465 | +0.098 | −0.066 | −0.187 | yes |
| ridge_subspace | −0.442 | +0.133 | −0.098 | −0.158 | yes |
| **scale_only** | **−0.211** | **+0.189** | **+0.094** | **−0.012** | **yes** |
| oracle_alpha | −0.424 | +0.378 | −0.058 | −0.133 | no (ceiling) |
| **scale + oracle_alpha** | **+0.046** | **+0.425** | **+0.119** | **+0.117** | no (ceiling) |

Fitted shrinkage: **0.455, 0.440, 0.428, 0.469** — remarkably consistent across
folds, and all well below 1.

**This overturns the hypothesis carried in from the previous phase.** I expected
template offset to be the main lever. It is not:

* Every basal-derived estimator of `alpha` fails. `direct_basal` is
  catastrophic (energy −6 to −17); `global_scalar` changes nothing; the
  rank-≤2 `ridge_subspace` map helps marginally in three folds and **hurts** in
  HepG2 — i.e. it is not generalising, it is fitting three points.
* **Even the oracle template leaves energy explained negative in three of four
  folds.** With the true `alpha` the residual is exactly `(4/3) gamma` (proved
  as a test), so the template was never the dominant error term.
* **Scale calibration alone — fully inference-available — does most of the
  work**, lifting every fold (K562 −0.467→−0.211, HepG2 −0.065→**+0.094**,
  Jurkat −0.184→**−0.012**, RPE1 +0.098→+0.189).
* Scale *and* template together make energy explained **positive in all four
  folds**, but the template half of that is not currently obtainable.

The reason the shrinkage is ≈0.44 rather than ≈1 is structural: the source mean
carries `beta_p − gamma[c*,p]/3` and cannot predict `gamma` at all, so the
variance-optimal point prediction is heavily shrunk. This is a property of the
problem, not a modelling choice.

---

## B. Source agreement validation

`source_agreement[p]` = mean pairwise Pearson among the three source responses
for perturbation `p`. Source contexts only.

Per held-out context, against reliability-normalised transfer success
(95 % bootstrap CI over perturbations):

| held out | vs raw success (Spearman) | vs normalised success (Spearman) | Pearson (normalised) | n |
|---|---|---|---|---:|
| K562 | +0.857 [0.841, 0.872] | **+0.729** [0.689, 0.763] | +0.715 | 924 |
| RPE1 | +0.884 [0.871, 0.895] | **+0.838** [0.816, 0.855] | +0.806 | 1,101 |
| HepG2 | +0.828 [0.810, 0.846] | **+0.622** [0.582, 0.663] | +0.594 | 842 |
| Jurkat | +0.801 [0.777, 0.823] | **+0.704** [0.661, 0.740] | +0.687 | 972 |

**The relationship holds in every fold independently**, with CIs far from zero.
It is not a pooling artefact.

### Confound control

| controlling for | K562 | RPE1 | HepG2 | Jurkat |
|---|---:|---:|---:|---:|
| source magnitude | +0.594 | +0.595 | +0.362 | +0.571 |
| source reliability | +0.530 | +0.293 | +0.150 | +0.540 |
| source cells/pert | +0.728 | +0.838 | +0.622 | +0.704 |
| target-gene basal expression | +0.716 | +0.817 | +0.602 | +0.679 |
| **all four jointly** | **+0.525** | **+0.301** | **+0.142** | **+0.547** |

**Source agreement adds information beyond simply detecting weak or noisy
perturbations**, but the margin is real and modest, and it is weakest exactly
where the target is noisiest (HepG2, +0.142). Source reliability and source
magnitude carry much of the shared signal; cell count and target-gene basal
expression carry none.

> **An over-adjustment I corrected.** The first confound set included
> `source_min_pair_agreement`, which drove the joint partial down to +0.12–0.19.
> That variable is the *minimum* of the very pairwise correlations whose *mean*
> is the exposure — adjusting for it removes the signal itself, not a confound.
> It is reported separately (+0.04 to +0.48) and excluded from the joint set.

---

## C. Pathway-level gamma recoverability

Gene sets: **MSigDB 2024.1.Hs**, Hallmark (50 sets) and C2:CP:REACTOME (1,736),
downloaded 2026-09-20, SHA-256 recorded in `data/provenance/msigdb/`. Used
exactly as released. Predeclared aggregation: **unweighted mean of the response
over each set's genes present in the frozen 6,640-gene space**, sets with <10
genes present dropped. 45 Hallmark and 874 Reactome sets qualify. No pathway
definition or aggregation choice was tuned.

Aggregation is linear, so decomposing then scoring equals scoring then
decomposing — asserted as a test, which keeps the comparison to gene level
exact rather than approximate.

Median `r(gamma_true, gamma_hat)` for the `basal_affine` baseline:

| held out | gene level | Hallmark | Reactome | Hallmark ceiling | Hallmark normalised |
|---|---:|---:|---:|---:|---:|
| K562 | +0.185 | **+0.570** | +0.437 | 0.874 | **0.652** |
| RPE1 | +0.002 | **+0.182** | +0.175 | 0.937 | 0.194 |
| HepG2 | +0.033 | +0.026 | +0.066 | 0.801 | 0.032 |
| Jurkat | +0.217 | **+0.510** | +0.369 | 0.852 | **0.599** |

**Yes — coarser biological resolution reveals transferable interaction structure
that is largely invisible gene by gene.** Hallmark gamma recovery is roughly
**3× the gene-level value** for K562 and Jurkat, and it lifts RPE1 from
essentially nothing (0.002) to a real if modest +0.182. Normalised against the
pathway measurement ceiling, K562 and Jurkat reach **0.60–0.65 of what is
attainable**.

Two honest qualifications:

* Part of the gain is reliability: pathway gamma reliability (`rho_full`
  0.64–0.88) is far above gene-level (0.28–0.55), because averaging ~100 genes
  suppresses noise. But the **ceiling-normalised** values also rise sharply
  (K562 0.32 → 0.65), so this is not purely a noise effect.
* **HepG2 remains at zero at every resolution.** Pathway aggregation does not
  manufacture transferable structure where none exists; it reveals it where it
  does. That is reassuring about the result rather than worrying.

Hallmark (45 coarse sets) outperforms Reactome (874 finer sets) consistently,
which is consistent with coarser aggregation buying more noise suppression.

---

## D. Candidate D definitions — characterised, not selected

Let `h1, h2` be independent target half responses and `A` a source-only
prediction. Under `h_i = L + e_i` with `e_i` independent of `L`, of each other,
and of `A`:

* `E<h1, h2> = ||L||²` — **reliable energy**, unbiased and *not* inflated by
  measurement noise.
* `E<h1 − A, h2 − A> = ||L − A||²` — **reliable residual energy**.

| candidate | formula | meaning of low / high | verdict |
|---|---|---|---|
| **D_unexplained_fraction** | `<h1−A, h2−A> / <h1, h2>` | 0 = transfer explains all reproducible signal; 1 = explains none; **>1 = worse than predicting zero** | **best behaved** |
| D_explained_fraction | `1 −` the above | success orientation | equivalent |
| D_corr_gap | `1 − r(A, h)/sqrt(rho_half)` | correlation shortfall | usable, but ignores scale errors entirely |
| D_raw_residual_fraction | same ratio on raw noisy quantities | — | **biased upward by noise** |
| D_raw_gamma_norm | `‖pooled − A‖` | — | **rejected: depth-biased** |

Synthetic validation (all in `tests/test_foundations.py`):

* Reliable energy recovers the planted `‖L‖²` across noise levels while the
  naive `‖h‖²` is inflated by >10 %.
* `D_unexplained_fraction` recovers planted values of `(1−share)²` for
  share = 0.25/0.5/0.9 to ±0.05.
* It correctly **exceeds 1** for an anti-correlated, over-scaled prediction —
  the regime we actually observe, so the target must not be clipped to [0,1].
* `D_raw_gamma_norm` **grows monotonically with noise** at fixed latent signal
  (0.3 → 1.0 → 2.5 noise), while the reliable fraction stays flat. This is the
  decisive reason not to use a raw residual norm as D.

**Edge case found and handled.** The per-pair ratio divides by a *noisy*
estimate of `‖L_p‖²`. For unreliable perturbations that denominator approaches
zero and the ratio destabilises — its mean even goes negative at high noise.
`pooled_unexplained_fraction` (sum of numerators over sum of denominators) is
stable and is the correct headline statistic; the per-pair value needs a
reliability filter before it is used to rank perturbations.

Observed values for source-mean transfer:

| held out | pooled unexplained | per-pair median | fraction >1 (worse than zero) |
|---|---:|---:|---:|
| K562 | 1.296 | 1.503 | 83.8 % |
| RPE1 | 0.752 | 0.900 | 39.8 % |
| HepG2 | 0.798 | 1.070 | 56.2 % |
| Jurkat | 1.007 | 1.196 | 66.5 % |

**No target is selected.**

---

## E. Decision report

### 1. Can basal controls recover alpha/template information?

**No, not usefully.** Direction is absent (r ≈ 0, sign-inconsistent), only
0.15–8.2 % of `alpha` lies in the source basal span, and the fitted global
scalar flips sign across folds. The only signal is a weak *magnitude*
association (Spearman +0.27–0.34). Every basal-derived estimator either fails
outright or fails to generalise.

### 2. Does template calibration fix negative energy explained?

**No — but scale calibration does, and it is available.** Even the oracle
template leaves energy negative in three of four folds. `scale_only`, fitted on
source contexts alone, lifts every fold and turns HepG2 positive and Jurkat to
≈0. Scale plus oracle template is positive in all four, but that template is not
obtainable today. **Scale is the lever; template offset was a red herring.**

### 3. Does source agreement predict transferability within each held-out context?

**Yes, in all four independently**: Spearman +0.622 to +0.838 against
reliability-normalised success, with bootstrap CIs far from zero.

### 4. Does it survive source magnitude/reliability confounding?

**Yes, attenuated and variable**: joint partial +0.525 (K562), +0.301 (RPE1),
+0.142 (HepG2), +0.547 (Jurkat). It adds information beyond detecting weak or
noisy perturbations, but the margin is modest and weakest in the noisiest
context.

### 5. Is pathway-level gamma materially more recoverable than gene-level?

**Yes, substantially — for three of four contexts.** Hallmark recovery is ~3×
gene-level for K562 (0.185 → 0.570) and Jurkat (0.217 → 0.510), and rescues RPE1
from nothing (0.002 → 0.182). Ceiling-normalised recovery reaches 0.60–0.65.
HepG2 stays at zero at every resolution.

### 6. Which candidate D is scientifically best behaved?

**`D_unexplained_fraction` = `<h1−A, h2−A> / <h1, h2>`**, reported pooled for
headline numbers and per-pair (with a reliability filter) for ranking. It is
unbiased under a stated noise model, validated against planted values,
represents the important >1 regime rather than clipping it, and is immune to the
depth bias that disqualifies a raw residual norm. Not selected as final.

### 7. What is actually available at Arc inference time?

Per the audited Arc bundle: **control cells for each unseen context**
(18,400 per context), **perturbation identity** (300 target genes), and the
18,533-gene panel. Nothing else about the target. So of everything studied here,
the usable inputs are: source-context responses (from public data), basal
control expression of the target, perturbation identity, and external gene
priors. **Source agreement and the fitted scale are both computable; the target
template is not.**

### 8. What is the lowest-capacity model justified by the evidence?

A **scale-calibrated conserved-transfer predictor plus a transferability
estimate**: source-mean response, shrunk by one scalar fitted by inner
leave-one-source-out, with per-perturbation confidence from source agreement.
That is two scalars and one descriptive feature — no learned context function.
Four contexts do not justify more, and the `ridge_subspace` result (helps in
three folds, hurts in the fourth) is a concrete demonstration of what
over-capacity looks like here.

### Recommendation

**B — pathway-level gamma model, combined with the scale calibration from A and
the transferability signal from B.**

The evidence, not preference, points here:

* Pathway gamma is the only place where a **large** amount of previously
  invisible transferable structure appeared (3× gene-level, 0.60–0.65 of
  ceiling). Nothing else in this study moved a number that far.
* Scale calibration is a cheap, fully available fix for the largest
  point-prediction failure mode, and should be folded in regardless of what else
  is built.
* Transferability (source agreement) is real and confound-robust but modest
  after adjustment; it is better as a **confidence output alongside** a
  predictor than as the primary research target.
* Template recovery from basal expression is **ruled out** by this study and
  should not be pursued further with four contexts.
* Exact gene-level gamma remains the weakest option and is now clearly
  dominated by the pathway-level route.

Concretely: predict gamma at Hallmark resolution, apply scale-calibrated
conserved transfer at gene level, and carry source agreement as a per-
perturbation confidence. Evaluate against pathway-level reliability ceilings.

## Caveats

* Four contexts, three sources per fold. Every context-level claim rests on
  n = 4. The `ridge_subspace` inconsistency is a live warning about capacity.
* HepG2 is unrecoverable at every resolution tested and drags all pooled
  statistics; it is also the shallowest context.
* The oracle-template rows are ceilings computed from held-out responses and are
  not achievable methods.
* Pathway gains are partly a reliability effect; the ceiling-normalised numbers
  are the honest comparison and they still rise sharply.
* The fitted shrinkage is slightly conservative by construction (inner folds
  have two sources, not three), documented in `fit_response_shrinkage`.

## Standing statement

Independent four-context study on public scPertEval data. **Not** a Molina &
Zhang reproduction. Arc contexts A/B/C were not used anywhere: not trained on,
not intersected with the response space, not used to choose preprocessing, and
no identity inference was attempted. **No D predictor, neural context encoder,
GNN, foundation-model predictor, generative model or Arc submission model was
trained.**
