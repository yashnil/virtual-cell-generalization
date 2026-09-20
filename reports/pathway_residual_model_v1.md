# Pathway residual model — v1

**Result: NEGATIVE.** The learned pathway correction does **not** improve
zero-shot pathway-response prediction beyond scale-calibrated conserved
transfer, in any of the four held-out contexts. The model's own inner selection
declined to correct in two of four folds, and an oracle sweep shows the maximum
attainable gain was **+0.003** Pearson. The result is reported in full.

**No deep model was built.** No Arc predictions, no gene-level lift-back, no D
predictor, no neural network, no foundation-model embeddings.

Date: 2026-09-20 · runtime 2.4 min + diagnostics · peak RSS 11.5 GB
Reproduce: `uv run python scripts/run_pathway_residual_model.py`
then `uv run python scripts/plot_pathway_residual_model.py`

All eight freeze records verified **before and after** modelling: protocol 4/4,
canonical 12/12, decomposition 21/21, zero-shot 14/14, foundations 12/12,
discovery 26/26, MSigDB 2/2, scPertEval data 4/4. Modelling code lives in a new
namespace (`virtual_cell.modelling`); no frozen analysis module was touched.

---

## Design

* **Outer folds**: the four immutable LOCO splits. The outer target's responses
  were used exactly once per fold, for evaluation.
* **Inner pseudo-LOCO**: each of the three source contexts held out in turn,
  trained on the other two. Model family, ridge penalty and shrinkage λ were all
  selected here. Predeclared inner criterion: median per-perturbation Pearson.
* **Representation**: frozen Hallmark 2024.1.Hs, 45 sets, the same predeclared
  unweighted-mean aggregation as the falsification study.
* **Baseline `B`**: scale-calibrated conserved transfer, scale fitted by inner
  leave-one-source-out on source contexts only (fitted values 0.52–0.59).
* **Grids**: λ ∈ {0, 0.25, 0.5, 0.75, 1.0}; ridge α ∈ {1, 10, 10², 10³, 10⁴}.
* **17 inference-available features**, standardised on training rows only.

### Leakage

Enforced by tests, not assertion: `test_features_never_read_the_target_responses`,
`test_inner_selection_never_touches_the_outer_target` and
`test_outer_prediction_uses_the_target_only_through_basal` each overwrite the
outer target's entire response matrix with noise and require **bit-identical**
features, selection and predictions. A companion test confirms the features *do*
move when a source changes, so the leakage tests cannot pass trivially.

### The residual is not gamma — and this matters

With scale `s`, source set of size 3:

```
R = Y - B = (4/3)·alpha_c + (1-s)·beta_p + (1 + s/3)·gamma[c,p]
```

asserted exactly in `test_residual_matches_the_documented_algebra`. Centring over
perturbations removes the `alpha` term exactly but **leaves `(1-s)·beta`**, which
at `s ≈ 0.55` is nearly half of the conserved effect. So a model can appear to
"work" by undoing the shrinkage rather than learning any interaction. This is why
the correlation of `R_hat` with `gamma` and with `beta` is reported separately
below — and it turned out to be the decisive diagnostic.

---

## 1. Which model family was selected in each outer fold?

| held out | family | α | λ | fitted scale |
|---|---|---:|---:|---:|
| K562 | **M0** (no correction) | — | **0.00** | 0.591 |
| RPE1 | M1 (deterministic) | — | 0.75 | 0.570 |
| HepG2 | M2 (pooled ridge) | 10⁴ | 0.75 | 0.524 |
| Jurkat | **M0** (no correction) | — | **0.00** | 0.589 |

**M3 (pathway-specific ridge) was not evaluated**, because the precondition —
"M2 shows genuine outer-context evidence" — was not met. M2 was selected in one
fold and made it slightly worse.

## 2. Which shrinkage was selected?

λ = 0 in two folds (K562, Jurkat) and λ = 0.75 in two (RPE1, HepG2). **In half
the folds the model used its honest option to decline to correct.**

## 3–6. Did the correction improve each context?

| held out | baseline r | corrected r | Δr | baseline energy | corrected energy | Δenergy | frac improved |
|---|---:|---:|---:|---:|---:|---:|---:|
| **K562** | 0.5004 | 0.5004 | **0.0000** | −1.448 | −1.448 | 0.000 | 0.00 |
| **RPE1** | 0.7058 | 0.6905 | **−0.0153** | +0.348 | +0.334 | −0.014 | 0.40 |
| **HepG2** | 0.6504 | 0.6471 | **−0.0034** | +0.352 | +0.334 | −0.018 | 0.51 |
| **Jurkat** | 0.5066 | 0.5066 | **0.0000** | −0.081 | −0.081 | 0.000 | 0.00 |

**3. K562 — no.** Declined to correct (λ=0).
**4. RPE1 — no, it hurt** (−0.0153 Pearson, −0.014 energy); only 40 % of
perturbations improved.
**5. HepG2 — no, it hurt slightly** (−0.0034, −0.018 energy).
**6. Jurkat — no.** Declined to correct.

### Capacity ladder (each family at its own inner-best λ)

| held out | M0 r | M1 r | M2 r | M0 energy | M1 energy | M2 energy |
|---|---:|---:|---:|---:|---:|---:|
| K562 | 0.5004 | 0.5004 | 0.5004 | −1.448 | −1.448 | −1.448 |
| RPE1 | **0.7058** | 0.6905 | 0.6904 | **+0.348** | +0.334 | +0.333 |
| HepG2 | **0.6504** | 0.6453 | 0.6471 | 0.352 | **+0.362** | +0.334 |
| Jurkat | 0.5066 | 0.5066 | 0.5066 | −0.081 | −0.081 | −0.081 |

**M0 is at least as good as M1 and M2 on Pearson in every fold.** More capacity
never helped.

### Is this a selection failure or a real absence of signal?

An **oracle λ sweep** (evaluation-only, uses the outer target, not a usable
method) settles it:

| held out | λ = 0 | best oracle | best λ | **max attainable gain** |
|---|---:|---:|---:|---:|
| K562 | 0.5004 | 0.5037 | 1.00 | **+0.0033** |
| RPE1 | 0.7058 | 0.7058 | 0.00 | **0.0000** |
| HepG2 | 0.6504 | 0.6529 | 0.50 | **+0.0024** |
| Jurkat | 0.5066 | 0.5066 | 0.00 | **0.0000** |

**λ = 0 is essentially optimal everywhere.** Nested LOCO did not fail to find a
gain; there was no gain to find. This is the single most important control in
the study, because "the inner folds cannot see the K562–Jurkat pairing" was the
obvious alternative explanation and it is now excluded.

## 7. Did the model shrink appropriately in the known-negative setting?

**Partially — and this is a genuine criticism.** HepG2 is the known-negative
(falsification showed Hallmark *below* its own null there). The model selected
λ = 0.75 for HepG2, i.e. it *did* apply a large correction where no recoverable
interaction exists, and performance degraded (−0.018 energy). It shrank correctly
in K562 and Jurkat but not in HepG2.

The forced-M2 diagnostic shows what it latched onto in HepG2:
`r(R_hat, gamma) = −0.054` (nothing) but `r(R_hat, beta-like) = +0.767`. **In the
known-negative context the correction is essentially pure beta re-prediction.**

## 8. Did Hallmark beat matched-random modelling?

**No.** 20 replicates, full nested selection re-run on each (so the null is not
handicapped):

| held out | Hallmark Δr | matched-random Δr | p | z |
|---|---:|---|---:|---:|
| K562 | +0.0000 | +0.0048 ± 0.0046 | 0.905 | −1.02 |
| RPE1 | −0.0153 | −0.0189 ± 0.0172 | 0.476 | +0.21 |
| HepG2 | −0.0034 | −0.0008 ± 0.0046 | 0.667 | −0.55 |
| Jurkat | +0.0000 | +0.0011 ± 0.0025 | 1.000 | −0.43 |

Hallmark is indistinguishable from random aggregation — **because neither
produces a gain**. This does *not* contradict the falsification battery, which
compared γ *recovery*, not response *prediction*. See §12.

## 9. Did Reactome reproduce the pattern?

Yes, qualitatively:

| held out | family | λ | Δr |
|---|---|---:|---:|
| K562 | M2 (α=1) | 0.25 | **+0.0120** |
| RPE1 | M2 (α=10⁴) | 0.50 | −0.0529 |
| HepG2 | M2 (α=10⁴) | 0.75 | −0.0164 |
| Jurkat | M0 | 0.00 | 0.0000 |

One small positive (K562 +0.012), two negatives, one abstention. No systematic
improvement, and the one positive is well inside the Hallmark null spread.

## 10. Which features mattered?

Standardised M2 coefficients, ranked by mean |coefficient|:

| feature | mean abs | sign consistent across folds? |
|---|---:|---|
| **source_mean** | 0.0815 | yes (all negative) |
| **weighted_source** | 0.0723 | no |
| source_mean_centred | 0.0070 | yes |
| weighted_source_centred | 0.0066 | yes |
| weighted_minus_mean | 0.0058 | no |
| pathway_size | 0.0030 | yes |
| source_reliability | 0.0027 | yes |
| target_basal_pathway | 0.0022 | yes |
| source_agreement | 0.0021 | **no** |
| basal_pathway_deviation | 0.0019 | **no** |
| source_magnitude | 0.0017 | yes |
| gene_in_pathway | 0.0008 | yes |
| source_cells | 0.0008 | yes |
| target_gene_basal | 0.0007 | yes |
| source_sd / range / sign_agreement | ≤0.0006 | no |

Two features dominate by an order of magnitude — `source_mean` and
`weighted_source`, i.e. **the raw conserved response itself**. Every biological
and target-context feature is an order of magnitude smaller.

* **Does source agreement add beyond source magnitude?** As a *model feature*,
  no — coefficient 0.0021, sign-inconsistent. (As a *confidence* signal it
  remains strong; see §11. These are different questions.)
* **Does target-gene basal expression contribute?** No — 0.0007, the smallest
  non-trivial coefficient.
* **Does basal-context similarity contribute?** No — `basal_pathway_deviation`
  is 0.0019 and sign-inconsistent.

No causal biological interpretation is drawn from any coefficient.

## 11. Did source-agreement confidence stay calibrated?

**Yes, cleanly, in all four folds.** Spearman(source agreement, corrected r):
K562 **+0.602**, RPE1 **+0.534**, HepG2 **+0.529**, Jurkat **+0.572**.

By quartile the ordering is strictly monotone in every context — e.g. HepG2
0.076 → 0.592 → 0.748 → 0.789, Jurkat 0.065 → 0.462 → 0.604 → 0.727 — and the
correction leaves it essentially unchanged. **The confidence signal is the one
component of this system that works as intended.**

## 12. Is pathway residual modelling justified beyond a conserved predictor?

**On this evidence, no — and the reason is specific and informative.**

The forced-M2 diagnostic (M2 fitted in every fold so `R_hat` is non-degenerate)
shows the model *does* find the interaction where it exists:

| held out | r(R_hat, gamma) | r(R_hat, beta-like) | deterministic r(gamma-hat, gamma) | ‖R_hat‖/‖R_true‖ |
|---|---:|---:|---:|---:|
| K562 | **+0.469** | −0.424 | +0.570 | 0.563 |
| RPE1 | +0.132 | **+0.838** | +0.182 | 0.201 |
| HepG2 | −0.054 | **+0.767** | +0.026 | 0.228 |
| Jurkat | **+0.428** | −0.491 | +0.510 | 0.465 |

So in K562 and Jurkat the learned correction recovers γ at r ≈ 0.43–0.47,
approaching the deterministic 0.51–0.57 from the falsification battery. **The γ
signal is real and the model finds it.** But:

1. **γ is a small share of the pathway response.** Improving a component worth
   ~21 % of energy, recovered at r ≈ 0.45, moves the total prediction by
   almost nothing.
2. **The residual target is contaminated by β.** `R` contains `(1-s)·beta`, and
   the model predicts that β component with the **wrong sign** in exactly the two
   contexts where it finds γ (−0.424, −0.491). The harmful β term cancels the
   helpful γ term, netting ≈ 0.
3. In RPE1 and HepG2 the correction is **almost entirely β** (+0.84, +0.77) with
   no γ, so it is pure noise-plus-harm.

This is the third distinct statement in a chain that must not be collapsed:
*γ exists and is reproducible* ≠ *γ is recoverable zero-shot* ≠ **predicting γ
improves the response prediction**. The first two are established; the third is
now falsified for this construction.

## 13. Next justified step

**Not** gene-level lift-back (nothing to lift back), **not** richer perturbation
priors (the biological features contributed ~0), and **not** a deeper model
(capacity never helped).

Two options are supported by the diagnostics:

**(a) Fix the residual target, then retest — cheap and well-motivated.** Define
the learning target against the *unshrunk* transfer (`s = 1`), where the algebra
gives `centred R = (4/3)·gamma` exactly and the β contamination vanishes by
construction. Then add the γ correction to the *scale-calibrated* baseline,
keeping the two calibrations decoupled. This directly removes the mechanism
identified in §12, and it is a one-line change to the target definition. It is
the only cheap way to find out whether the β contamination was the whole
problem.

**(b) Otherwise, stop pathway modelling and pivot to transferability / D.**
Source-agreement confidence is the only component here that demonstrably works
(§11, Spearman +0.53 to +0.60 across all four contexts, monotone by quartile).
A calibrated "will conserved transfer work for this perturbation in this
context" predictor is well-posed and immediately useful for the Arc setting,
where knowing which predictions to trust has direct value.

**Recommendation: run (a) once as a narrow, pre-specified retest. If it does not
produce an outer-fold gain in at least K562 or Jurkat, stop pathway modelling
and move to (b).** Do not iterate on (a) beyond a single attempt — the oracle
sweep caps any achievable gain from this family of corrections at ~+0.003 unless
the target definition genuinely changes what is being learned.

## Limitations

* **Training/prediction mismatch inherent to n = 4.** Inner pseudo-targets are
  built from **two** source contexts while the outer prediction uses **three**,
  so features like `source_mean` have different noise and attenuation at fit and
  at application time. Unavoidable with four contexts; partially absorbed by
  per-fold standardisation; a real limitation.
* The fitted scale for a two-source inner fold rests on single-source inner-inner
  folds and is correspondingly crude.
* Baseline energy explained is strongly negative for K562 (−1.45) even at
  pathway level, so that fold's energy metric is dominated by a scale/template
  mismatch the correction was never going to fix.
* 20 null replicates (not 100) because each requires a full nested-LOCO refit.
* Four contexts. Every per-context number is n = 1.

## Standing statement

Independent four-context study on public scPertEval data. **Not** a Molina &
Zhang reproduction. Arc contexts A/B/C were not used anywhere. No Arc
predictions were generated, no single-cell generative model, no neural network,
no foundation-model embeddings, no D predictor, and no gene-level lift-back.
