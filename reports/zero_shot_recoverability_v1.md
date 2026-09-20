# Zero-shot recoverability — diagnostic v1

**Diagnostic only. No model was trained.** Every baseline is a closed-form
function of source-context responses and basal (control) expression. This is an
independent four-context study, **not** a Molina & Zhang reproduction.

Date: 2026-09-19 · runtime 7.4 min · peak RSS 18.3 GB
Reproduce: `uv run python scripts/run_zero_shot_recoverability.py`
then `uv run python scripts/plot_zero_shot_recoverability.py`

All canonical and robustness artifacts were frozen first
(`data/provenance/scperteval/decomposition_phase_freeze.txt`, **21/21 verify**)
and none were modified.

---

## 1. Frozen LOCO splits

`data/splits/loco_v1/` — four folds, each with a SHA-256, plus
`loco_manifest.yaml` recording the leakage contract verbatim.

| fold | sources | target |
|---|---|---|
| 1 | RPE1 + HepG2 + Jurkat | K562 |
| 2 | K562 + HepG2 + Jurkat | RPE1 |
| 3 | K562 + RPE1 + Jurkat | HepG2 |
| 4 | K562 + RPE1 + HepG2 | Jurkat |

**Allowed as prediction input:** basal/control cells of the held-out context,
perturbation identity, source-context perturbation responses, external priors
free of target outcomes.
**Forbidden:** any target perturbation-response measurement; feature selection,
normalisation or hyperparameters fitted to target responses; **and the
four-context `beta` and `gamma`**, because both are computed from all four
contexts and therefore contain the held-out response.

Enforcement is tested, not just asserted:
`test_source_mean_never_reads_the_target_context` replaces the entire target row
with noise and requires every baseline to be bit-identical.

## 2/3. Algebra of source-only transfer, and the three quantities

With `delta[c,p] = mu + alpha_c + beta_p + gamma[c,p]` and the balanced-design
side conditions `sum_c alpha_c = 0`, `sum_c gamma[c,p] = 0`, averaging over the
three sources gives `mean_S alpha = -alpha_{c*}/3` and
`mean_S gamma = -gamma[c*,p]/3`, so

```
A[p] = mean_{c in S} delta[c,p] = mu - alpha_{c*}/3 + beta_p - gamma[c*,p]/3
```

Two consequences, both proved as tests:

* **Identity 1** — `delta[c*,p] - A[p] = (4/3)(alpha_{c*} + gamma[c*,p])`, and
  centring the residual over perturbations leaves **exactly `(4/3) gamma[c*,p]`**.
  So "predict the centred transfer residual" and "predict gamma" are the *same
  problem* up to a constant. This is why gamma may never be a predictor input.
* **Identity 2** — the source mean carries `-gamma[c*,p]/3`. Conserved transfer
  is **mildly anti-correlated with the held-out interaction by construction**;
  it is not gamma-neutral.

The three quantities the brief asks to separate:

| | definition | uses target response? | role |
|---|---|---|---|
| **A** | four-context `gamma[c*,p]` | **yes** | evaluation / descriptive **only** |
| **B** | `A[p] = mean_S delta[c,p]` and the other baselines | no | available at inference |
| **C** | `delta[c*,p] - A[p]`, centred = `(4/3) gamma[c*,p]` | yes | evaluation target |

## 4. Reliability-aware evaluation

Additive independent-noise model `h = L + e`, with `e` independent of `L` and
between halves:

* `corr(h1, h2) = Var(L)/(Var(L)+Var(e)) = rho_half` — the split-half
  correlation **is** the reliability of one half.
* A *perfect* latent predictor achieves `corr(L, h1) = sqrt(rho_half)`.
  **The ceiling is the square root of reliability, not reliability.**
* Therefore `rho_latent = corr(P, h) / sqrt(rho_half)`.
* Pooling both halves halves the noise variance: `rho_full = 2 rho_half/(1+rho_half)`
  (Spearman-Brown).

Each step is validated on synthetic data with a known latent signal
(`tests/test_loco.py`): reliability recovers the planted noise ratio to ±0.02 at
three noise levels; the achieved ceiling matches `sqrt(rho)` and is shown to
differ from `rho` by >0.05; and disattenuation recovers planted latent
correlations of 0.3/0.6/0.9 to ±0.05. Disattenuation returns NaN below
`rho_half = 0.05`, where the denominator is both tiny and badly estimated.

Measured target reliabilities (10 independent split-half repeats):

| held-out | rho_half | rho_full | ceiling = sqrt(rho_full) | usable fraction |
|---|---:|---:|---:|---:|
| K562 | 0.241 | 0.388 | 0.623 | 73.1 % |
| RPE1 | **0.571** | **0.727** | **0.852** | 87.1 % |
| HepG2 | 0.226 | 0.368 | 0.606 | 66.6 % |
| Jurkat | 0.276 | 0.432 | 0.657 | 76.9 % |

## 5. Zero-shot baselines

Basal similarity (control profiles, Pearson) and the derived weights — **no
target response anywhere**:

| held-out | basal similarity to sources | nearest | basal-affine weights |
|---|---|---|---|
| K562 | RPE1 0.892, HepG2 0.908, **Jurkat 0.934** | Jurkat | RPE1 0.085, HepG2 0.328, **Jurkat 0.588** |
| RPE1 | K562 0.892, **HepG2 0.915**, Jurkat 0.898 | HepG2 | K562 0.110, **HepG2 0.554**, Jurkat 0.337 |
| HepG2 | K562 0.908, **RPE1 0.915**, Jurkat 0.891 | RPE1 | K562 0.387, **RPE1 0.506**, Jurkat 0.108 |
| Jurkat | **K562 0.934**, RPE1 0.898, HepG2 0.891 | K562 | **K562 0.626**, RPE1 0.277, HepG2 0.097 |

Median per-perturbation Pearson against the held-out response (95 % bootstrap CI
over perturbations):

| baseline | K562 | RPE1 | HepG2 | Jurkat | pooled |
|---|---:|---:|---:|---:|---:|
| **source mean** | 0.277 [0.254,0.298] | 0.357 [0.332,0.376] | 0.326 [0.300,0.354] | 0.283 [0.260,0.299] | **0.306** |
| basal affine | 0.277 | 0.355 | 0.316 | **0.292** | 0.310 |
| basal simplex | 0.278 | 0.355 | 0.316 | 0.292 | 0.310 |
| nearest basal | 0.212 | 0.269 | 0.269 | 0.212 | 0.239 |
| best single source | 0.212 | 0.269 | 0.269 | 0.212 | — |

**Averaging beats picking.** The source mean and the basal-weighted combinations
are indistinguishable from each other and clearly better than any single source,
including the basal-nearest one. With only three source contexts there is no
evidence that a learned function of context would add anything on raw response
prediction.

### The most important negative result: energy explained is negative

| baseline | K562 | RPE1 | HepG2 | Jurkat | pooled |
|---|---:|---:|---:|---:|---:|
| source mean | **−0.467** | +0.098 | −0.065 | −0.184 | **−0.138** |
| basal affine | −0.438 | +0.075 | −0.119 | −0.113 | −0.146 |
| nearest basal | −0.979 | −0.158 | −0.539 | −0.236 | −0.472 |

Conserved transfer has a solidly positive **direction** (r ≈ 0.31, cosine ≈ 0.31)
but as a *point prediction of the raw response it is worse than predicting
zero* in three of four folds. The held-out template `alpha_{c*}` is not
recoverable from source responses, and the source mean substitutes
`-alpha_{c*}/3` for it, so the prediction sits in the wrong place with roughly
the wrong scale. Any future model that is scored on raw magnitude — Arc's
metrics included — must address the template and the scale, not just direction.

## 6. Gamma recoverability — the central question

Gamma reliability is much lower than the response as a whole, as expected:

| held-out | gamma rho_half | gamma rho_full | ceiling |
|---|---:|---:|---:|
| K562 | 0.225 | 0.368 | 0.607 |
| RPE1 | 0.381 | 0.552 | 0.743 |
| HepG2 | 0.166 | 0.285 | 0.534 |
| Jurkat | 0.206 | 0.342 | 0.585 |

> **A bug found and fixed during this run.** The first implementation replaced
> only the *target* row with a half while leaving the three source contexts at
> full data. Because `gamma[c,p]` is a function of all four contexts, the source
> contribution was then identical in both "halves" and the correlation was
> inflated to 0.90–0.94 — irreconcilable with the canonical 49.5 % gamma
> reproducibility, which is what exposed it. It also averaged gamma across
> repeats before correlating, averaging the noise away a second time. The fixed
> version splits **every** context and averages per-repeat correlations. The
> numbers above are from the fixed version.

Recovery of gamma by the zero-shot baselines, `r(gamma_true, gamma_hat)`:

| held-out | basal affine | nearest basal | ceiling | disattenuated (affine) | frac > 0 |
|---|---:|---:|---:|---:|---:|
| K562 | **+0.185** [0.174,0.196] | +0.175 | 0.607 | 0.323 | 85.8 % |
| RPE1 | +0.002 [−0.009,0.016] | +0.023 | 0.743 | 0.004 | 50.8 % |
| HepG2 | +0.033 [0.026,0.044] | **−0.022** | 0.534 | 0.068 | 59.9 % |
| Jurkat | **+0.217** [0.210,0.229] | **+0.241** | 0.585 | 0.395 | 93.3 % |

**Gamma is partially recoverable zero-shot for K562 and Jurkat, and not at all
for RPE1 and HepG2.** The K562/Jurkat effect is not marginal — the bootstrap CI
excludes zero by a wide margin and 86–93 % of individual perturbations have
positive recovery. For RPE1 the median is 0.002 and the CI straddles zero.

### Why: gamma is shared only between genuinely similar contexts

Cross-context gamma correlation. **The null is not zero:** `sum_c gamma[c,p] = 0`
forces an expected pairwise correlation of `-1/(C-1) = -0.333`.

| pair | median r(gamma) | excess over null | basal similarity |
|---|---:|---:|---:|
| **K562–Jurkat** | **−0.137** | **+0.196** | **0.934** |
| K562–HepG2 | −0.300 | +0.033 | 0.908 |
| RPE1–HepG2 | −0.345 | −0.011 | 0.915 |
| HepG2–Jurkat | −0.408 | −0.074 | 0.891 |
| K562–RPE1 | −0.411 | −0.077 | 0.892 |
| RPE1–Jurkat | −0.414 | −0.081 | 0.898 |

Only **one** pair exceeds the null appreciably: K562–Jurkat, which is also the
most basally similar pair. Spearman(basal similarity, excess) = **+0.714**
across the six pairs — suggestive, but **n = 6 and it rests almost entirely on
one point**, so it is a hypothesis, not a finding.

This explains the fold pattern exactly: K562 and Jurkat each have a genuinely
gamma-correlated partner in their source set, and recover gamma. RPE1 and HepG2
do not, and recover nothing.

## 7. Transferability — descriptive only, no target defined

Spearman against conserved-transfer success, pooled over folds:

| candidate feature | vs raw r | vs reliability-normalised r | available at inference? |
|---|---:|---:|:--:|
| **source agreement** (mean pairwise r among the 3 source responses) | **+0.828** | **+0.726** | **yes** |
| target reliability rho_half | +0.864 | +0.472 | no |
| response magnitude ‖delta‖ | +0.730 | +0.512 | no |
| ‖gamma‖ | +0.477 | +0.303 | no |
| target-gene basal expression | +0.173 | +0.094 | yes |
| basal context distance | +0.094 | −0.036 | yes |
| target cells per perturbation | −0.023 | −0.174 | no |

**Source agreement is the standout.** It is the strongest single predictor
(+0.828), it survives reliability normalisation nearly intact (+0.726 — so it is
not merely tracking measurement noise), and unlike reliability or ‖delta‖ **it is
computable at inference time from source contexts alone**. Target reliability
drops from +0.864 to +0.472 under normalisation, confirming that most of its raw
association was the measurement-noise artefact the correction is designed to
remove.

Basal context distance is useless for *per-perturbation* transferability
(+0.094 → −0.036), even though basal similarity matters at the *context* level
(§6). Cells per perturbation is flat — above the inherited 30-cell floor, depth
does not predict transfer success.

By reliability stratum (source mean):

| rho_half | n | raw r | normalised r | source agreement | ‖gamma‖ |
|---|---:|---:|---:|---:|---:|
| < 0.1 | 1,341 | 0.041 | 0.288 | 0.036 | 3.33 |
| 0.1–0.3 | 1,120 | 0.234 | 0.423 | 0.168 | 3.50 |
| 0.3–0.6 | 1,213 | 0.407 | 0.523 | 0.262 | 3.88 |
| > 0.6 | 1,382 | 0.592 | 0.637 | 0.407 | 4.90 |

### Three distinct classes

**Transferable** (reliable and well predicted) — POLRMT, SMG5, TFAM: normalised
r 0.86–0.92 with source agreement 0.56–0.85. Core mitochondrial-transcription and
mRNA-surveillance machinery, behaving the same way everywhere.

**Reliably context-specific** (highly reliable, near-zero or negative transfer) —
BCR in K562 (rho 0.917, r −0.074), NSMCE2 and EXOSC1 and SHQ1 and GRWD1 in RPE1,
GAB2 in K562, IPO7 in Jurkat. These are the scientifically interesting cases:
the response is genuinely reproducible and genuinely not shared. Source agreement
is near zero for all of them, so they are *flaggable in advance*.

**Large but unreliable** — PIAS4, PSMC1, ANAPC1, TPRKB, CHEK1, CDC16, TTK, SCFD1,
all in Jurkat with 31–36 cells and ‖delta‖ ≈ 5.8–6.8. Apparent effect near the
top of the distribution, reliability ≈ 0.00–0.05. These are the depth-biased
‖delta‖ artefact from the v1/sensitivity phases showing up concretely. Without
the reliability machinery they would be mistaken for strong context-specific
biology.

Counts: 1,768 of 5,056 pairs reliable (rho > 0.5); 1,341 unreliable (rho < 0.1).

## 8. Figures

`outputs/zero_shot_v1/` (git-ignored, regenerable):

| figure | content |
|---|---|
| `fig1_transfer_by_context.png` | transfer performance by held-out context and baseline, three metrics |
| `fig2_observed_vs_normalised.png` | observed vs reliability-normalised, and the `sqrt(rho)` ceiling |
| `fig3_error_vs_reliability.png` | transfer success by reliability stratum; magnitude vs reliability |
| `fig4_gamma_recoverability.png` | gamma recovery by context against its ceiling, and the distribution |
| `fig5_source_agreement.png` | source agreement vs transfer success, before and after normalisation |
| `fig6_basal_vs_gamma_similarity.png` | basal similarity vs cross-context gamma sharing, against the forced null |
| `fig7_example_classes.png` | the three perturbation classes with named examples |

---

## 9. Decision point

### A. How much held-out response does conserved transfer capture?

Median per-perturbation Pearson **0.306** pooled (0.277–0.357 by fold); after
reliability normalisation, **0.49–0.64** of the attainable latent correlation.
So conserved transfer captures roughly **half of what is measurable in
direction** — but **explains negative response energy** (−0.138 pooled) because
the held-out template and scale are wrong. Direction yes; calibrated magnitude no.

### B. How much reproducible signal remains unexplained?

Most of it. Even at the >0.6 reliability stratum, normalised r is 0.637, leaving
~60 % of reliable variance unexplained. Gamma holds 21 % of total response energy
(canonical v1) of which ~50 % is reproducible, and conserved transfer recovers
essentially none of it by construction (Identity 2 makes it slightly
anti-correlated).

### C. Is exact gene-level gamma measurably predictable zero-shot?

**Partially, and only when a genuinely similar context is available.** K562
r = 0.185 and Jurkat r = 0.217 (disattenuated 0.32 and 0.40) against ceilings of
0.61 and 0.59 — real, well above zero, but recovering only about a third to a
half of what the ceiling allows. RPE1 (0.002) and HepG2 (0.033) show nothing.
**This is not "gamma is unpredictable", and it is not "gamma is predictable" —
it is context-dependent**, and with four contexts we cannot say more.

### D. Does context similarity help?

**At the context level, yes; at the perturbation level, no.**
Basal similarity correctly identifies K562↔Jurkat as the one gamma-sharing pair
(Spearman +0.714 over 6 pairs, driven by one point — weak evidence). But
basal-weighted baselines do **not** beat the plain source mean on raw response
(0.310 vs 0.306), and basal context distance is useless for per-perturbation
transferability (−0.036 normalised). Similarity tells you *whether* gamma will
transfer, not *which* perturbations.

### E. Are there observable features that identify where beta transfer succeeds?

**Yes — source agreement.** Spearman +0.828 raw and **+0.726 after reliability
normalisation**, and it is computable at inference from source contexts alone.
Target-gene basal expression (+0.094) and basal context distance (−0.036) are
not useful. This is the single most actionable result in the study.

### F. What is the strongest next research target?

Ranked by what the data support, not by novelty:

1. **Transferability / D prediction — strongest.** Source agreement already
   delivers +0.726 normalised association with transfer success using only
   inference-time information, with no model at all. A calibrated predictor of
   "will conserved transfer work for this perturbation in this context" is
   well-posed, evaluable, and immediately useful. The reliability machinery built
   here is exactly what is needed to define its target without the noise artefact
   that would otherwise dominate it.
2. **Template and scale calibration — underrated, and arguably first.** Negative
   energy explained is the largest single failure mode found, it is
   context-level rather than perturbation-level (so only four effective data
   points, but the target's own control cells are available at inference), and it
   blocks any metric that scores magnitude. Cheap to attempt, high leverage.
3. **Pathway-level gamma — plausible.** Gene-level gamma is near its noise floor
   for most pairs; aggregating to pathways should raise reliability and may
   expose structure that gene-level correlation cannot resolve. Untested here.
4. **Exact gene-level gamma prediction — weakest.** It works only where a similar
   partner context exists (2 of 4 folds), recovers a third to a half of a ceiling
   that is itself only 0.53–0.74, and we have four contexts. High risk of fitting
   context idiosyncrasy.

**Recommended: (1) with (2), and (3) as a cheap diagnostic before committing to
(4).**

## Caveats

* **Four contexts, three sources per fold.** Every context-level statement rests
  on n = 4 (or n = 6 pairs). The basal-similarity/gamma-sharing link leans on a
  single pair. Do not fit anything of appreciable capacity to context.
* K562–Jurkat are both suspension leukaemia lines and RPE1/HepG2 are adherent
  epithelial/hepatocyte; the "similar partner" effect is biologically plausible
  but confounded with dataset of origin (Replogle vs Nadig) and cannot be
  separated with four contexts.
* Disattenuation is undefined at low reliability; 13–33 % of pairs per fold are
  excluded from normalised numbers. Medians are over the usable subset.
* Reliability was measured with 10 split-half repeats, fewer than the canonical
  50; the resulting medians are stable but individual pair estimates are noisy.

## Standing statement

Independent four-context diagnostic on standardized public scPertEval data.
Not a Molina & Zhang reproduction. Arc contexts A/B/C were not used anywhere.
**No model was built** — no Ridge, MLP, gamma predictor, D predictor, context
encoder, GNN, foundation-model embedding, conditional generative model, or Arc
submission model.
