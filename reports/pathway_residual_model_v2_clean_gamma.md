# Pathway residual model v2 — clean (beta-free) gamma target

**Result: NEGATIVE. The predeclared stopping rule has fired. Pathway modelling
is TERMINATED. There will be no v3.**

Removing the beta contamination worked exactly as the algebra predicted — gamma
prediction improved substantially (K562 r = 0.47 → **0.78**, Jurkat 0.43 →
**0.62**). It still did not improve the actual zero-shot response prediction in
any context. At the theory-motivated coefficient it made every context
materially worse.

Date: 2026-09-20 · runtime ~3 min · Hallmark 45 sets
Reproduce: `uv run python scripts/run_pathway_residual_model_v2.py`

v1 is permanently frozen (`pathway_residual_v1_freeze.txt`, 15 digests) and its
conclusion is retained verbatim: *"Pathway residual model v1 did not improve
outer zero-shot response prediction."* All nine freeze records verify before and
after v2.

**Exactly one thing changed from v1: the training target.** Splits, features,
baseline, scale fitting, leakage contract, families (M0/M1/M2 only), grids and
metrics were reused verbatim — asserted by
`test_features_are_identical_to_v1`, which requires the v2 design matrix to be
**byte-identical** to v1's, paired with `test_v2_target_differs_from_v1_target`
so the claim cannot be vacuous.

---

## The algebra, verified numerically

With unshrunk transfer `A[p] = mean_{s∈S} Y[s,p]`, `beta_p` appears with
coefficient 1 in both `Y[c]` and `A` and **cancels exactly, before any
centring, for any source-set size**:

```
R1 = Y[c] - A = (alpha_c - mean_S alpha) + (gamma[c,p] - mean_S gamma[·,p])
centre_p(R1)  = gamma[c,p] - mean_S gamma[·,p]
```

* **Outer case** (3 sources): `= (4/3)·gamma[c,p]` — verified to 1e-10 on
  planted balanced data (`test_unshrunk_centred_residual_is_exactly_four_thirds_gamma`).
* **Inner pseudo-target case** (2 sources): `= gamma[c] − (gamma[a]+gamma[b])/2`
  — beta-free, but a *different* linear combination and a different scale. This
  mismatch is inherent to four contexts and is documented, not hidden
  (`test_two_source_clean_target_is_beta_free_but_not_a_multiple_of_gamma`).

Multiplying beta by 50 leaves the clean target unchanged to 1e-9
(`test_clean_target_removes_beta_exactly`), while the v1 target shifts by
exactly `(1−s)·Δbeta` (`test_scaled_residual_still_contains_beta_when_s_is_not_one`).

**Theory coefficient.** `B = template + s·centre(A)` leaves a gamma deficit of
`(1 + s/3)·gamma`, so a correction predicting `(4/3)·gamma` should enter at
`lambda_theory = (3 + s)/4`. Verified exactly in
`test_lambda_theory_exactly_corrects_the_gamma_deficit`. With s ≈ 0.52–0.59 this
is **0.88–0.90**.

---

## 1. Did removing beta contamination improve gamma prediction itself?

**Yes, substantially — the fix worked.**

| held out | r(R̂, γ) v1 | r(R̂, γ) **v2** | r(R̂, β-like) v1 | r(R̂, β-like) v2 | ‖R̂‖/‖R_true‖ v2 |
|---|---:|---:|---:|---:|---:|
| K562 | +0.469 | **+0.778** | −0.424 | −0.927 | 0.98 |
| Jurkat | +0.428 | **+0.623** | −0.491 | −0.917 | 0.77 |
| RPE1 | +0.132 | +0.182 | +0.838 | +0.537 | 0.22 |
| HepG2 | −0.054 | +0.125 | +0.767 | −0.920 | 0.37 |

The correction now also has close to the right **magnitude** for K562
(‖R̂‖/‖R_true‖ = 0.98) rather than v1's 0.56. The strongly negative β-like
correlations are expected, not a defect: `beta-like = beta − gamma/3`, so a
predictor that is nearly pure gamma anti-correlates with it.

## 2–5. Did it improve actual response prediction?

**No, in all four contexts.** Hallmark, per outer context, selected λ:

| held out | sel. family | sel. λ | λ_theory | baseline r | corrected r | **Δr** | bootstrap 95 % CI | Δenergy | frac improved |
|---|---|---:|---:|---:|---:|---:|---|---:|---:|
| **K562** | M2 (α=10) | 0.25 | 0.90 | 0.5004 | 0.4968 | **−0.0036** | [−0.0134, +0.0069] | −0.319 | 0.50 |
| **RPE1** | M1 | 0.25 | 0.89 | 0.7058 | 0.6885 | **−0.0173** | **[−0.0243, −0.0107]** | −0.016 | 0.40 |
| **HepG2** | M2 (α=10³) | 0.50 | 0.88 | 0.6504 | 0.6462 | **−0.0042** | [−0.0195, +0.0039] | −0.026 | 0.42 |
| **Jurkat** | M2 (α=10⁴) | 0.25 | 0.90 | 0.5066 | 0.5006 | **−0.0060** | [−0.0133, +0.0051] | −0.006 | 0.52 |

**2. K562 — no** (−0.0036, CI includes zero).
**3. RPE1 — no, significantly harmful** (−0.0173, CI **excludes** zero).
**4. HepG2 — no** (−0.0042, CI includes zero).
**5. Jurkat — no** (−0.0060, CI includes zero).

Three folds are statistically indistinguishable from no change; one is
significantly worse. **Not a single fold improved.**

## 6. Selected λ vs λ_theory — the decisive diagnostic

| held out | selected λ | λ_theory | r at selected λ | r at λ_theory | **λ_theory vs baseline** |
|---|---:|---:|---:|---:|---:|
| K562 | 0.25 | 0.90 | 0.4968 | 0.4829 | **−0.0175** |
| RPE1 | 0.25 | 0.89 | 0.6885 | 0.6055 | **−0.1003** |
| HepG2 | 0.50 | 0.88 | 0.6462 | 0.6244 | **−0.0260** |
| Jurkat | 0.25 | 0.90 | 0.5006 | 0.4599 | **−0.0467** |

Inner selection consistently chose λ ≈ 0.25–0.5, **far below** the
mathematically correct 0.88–0.90 — and applying the theoretically correct amount
makes every context substantially worse, RPE1 catastrophically so.

This is the key finding. If `R̂` were an accurate estimate of `(4/3)γ`, λ_theory
would be optimal by construction. That it is strongly harmful means the
correction is **directionally right but not accurate enough per element**: at
r(R̂, γ) = 0.78, roughly 40 % of the correction's variance is error, and that
error is injected into a response whose γ component is small. The inner folds
correctly detected this and shrank λ — they simply could not shrink it to a
value that helped, because no positive λ helps.

## 7. Did Hallmark beat matched-random modelling?

**No.** 20 replicates, full nested selection re-run on each:

| held out | Hallmark Δr | matched-random Δr | p | z |
|---|---:|---|---:|---:|
| K562 | −0.0036 | −0.0011 ± 0.0034 | 0.952 | −0.76 |
| RPE1 | −0.0173 | −0.0139 ± 0.0126 | 0.714 | −0.27 |
| HepG2 | −0.0042 | −0.0034 ± 0.0090 | 0.667 | −0.09 |
| Jurkat | −0.0060 | −0.0009 ± 0.0030 | 0.905 | −1.74 |

Hallmark is indistinguishable from random aggregation, and in every fold the
point estimate is *worse* than the null mean.

## 8. Did Reactome support the same conclusion?

**Yes, fold for fold.**

| held out | sel. family | sel. λ | Δr | bootstrap CI |
|---|---|---:|---:|---|
| K562 | M1 | 0.25 | −0.0004 | [−0.0016, +0.0006] |
| RPE1 | M1 | 0.25 | −0.0169 | **[−0.0251, −0.0092]** |
| HepG2 | M1 | 0.25 | −0.0050 | [−0.0144, +0.0030] |
| Jurkat | **M0** | **0.00** | 0.0000 | [0, 0] |

No improvement anywhere; RPE1 again significantly harmful; Jurkat abstained
entirely. An independent ontology reaches the same verdict.

## 9. Did source-agreement confidence remain useful?

**Yes — unchanged and still the one thing that works.** Spearman(source
agreement, final prediction quality): K562 **+0.605**, RPE1 **+0.538**, HepG2
**+0.498**, Jurkat **+0.580**, monotone across quartiles in every context (e.g.
Jurkat 0.061 → 0.452 → 0.597 → 0.727). The correction leaves it essentially
untouched.

## 10. Is pathway modelling justified, or terminated?

**Terminated.** The predeclared stopping rule in §G fires on all three of its
disjunctive conditions:

| predeclared trigger | fired? |
|---|---|
| K562/Jurkat show essentially no usable response gain | **yes** (−0.0036, −0.0060; both CIs include zero) |
| gains are tiny / within bootstrap uncertainty | **yes** (three of four CIs include zero; the fourth is negative) |
| matched-random pathways perform equivalently | **yes** (p = 0.67–0.95, all z ≤ 0) |

Any one would have been sufficient. All three fired.

The scientific reading is precise and worth stating carefully, because it is
**not** "pathway biology is irrelevant". The chain established across phases is:

1. γ exists and is reproducible — **true** (decomposition + sensitivity).
2. γ is recoverable zero-shot at pathway resolution, genuinely from biology
   rather than aggregation — **true** (falsification battery, p = 0.010 against
   100 geometry-preserving nulls).
3. A learned γ correction improves the actual response prediction — **false**,
   now tested twice, with the contamination mechanism identified and removed.

Step 3 fails for a quantitative reason, not a conceptual one: γ is ~21 % of
response energy and is recovered at r ≈ 0.6–0.8, so the error injected exceeds
the signal added. **A correction must be far more accurate than "substantially
correlated" before it pays for itself.**

## 11. Next primary direction

**Transferability / confidence, not another gamma model.**

Stated explicitly as required: *the next primary direction is a
transferability/confidence model built on scale-calibrated conserved transfer
plus source agreement — not another gamma model.*

The justification is empirical and consistent across every phase:

* **Scale-calibrated conserved transfer** is the best available point predictor
  and is fully inference-available. It was never beaten by any correction.
* **Source agreement** predicts transfer success at Spearman +0.50 to +0.61
  within every outer context, survives joint confound control (+0.14 to +0.55),
  is monotone by quartile, and is computable at inference from source contexts
  alone. It is the only component of this system that has worked in every
  experiment.
* Knowing *which* predictions to trust has direct value in the Arc setting,
  where a per-perturbation confidence can gate or weight submissions.

**No pathway model v3.** Revisiting pathway-level γ would require new external
evidence — more contexts, or a materially more accurate γ estimator — not
another iteration on these four cell lines.

## Limitations

* Four contexts; every per-context number is n = 1.
* Inner pseudo-targets use two sources while the outer prediction uses three, so
  the training target is `gamma[c] − (gamma[a]+gamma[b])/2` while the outer
  target is `(4/3)·gamma[c]` — different linear combinations and scales. λ
  absorbs some of this, imperfectly. Unavoidable at n = 4 and a genuine
  limitation on what any nested-LOCO model here can learn.
* One inner-selection miss worth recording: for RPE1 the capacity ladder shows
  M2 at λ=0.25 would have scored 0.7206 versus the 0.7058 baseline (+0.0148),
  but inner selection chose M1 (0.6885). So in one of four folds a better option
  existed and was not selected. This does not change the verdict — a single
  +0.015 in one fold, unreplicated and not selectable, is exactly the kind of
  result the predeclared rule exists to discount — but it is the strongest
  counter-evidence available and is recorded rather than omitted.
* 20 null replicates (not 100), because each requires a full nested-LOCO refit.

## Standing statement

Independent four-context study on public scPertEval data. **Not** a Molina &
Zhang reproduction. Arc contexts A/B/C were not used anywhere. No Arc
predictions, no gene-level lift-back, no D predictor, no neural network, no
foundation-model embeddings, no generative model. M3 was not evaluated.
