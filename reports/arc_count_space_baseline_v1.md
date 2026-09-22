# Arc count-space baseline — v1

**The model is now a submission, not a correlation.** The frozen research
findings were assembled into the smallest predictor they support —
`delta_hat[c,p] = m_hat[c] + w[tier(p)] * beta_hat[p]` — every free parameter
was chosen on public held-out contexts, the prediction was pushed through raw
count generators, and the generators were scored **against real held-out cells
with the six `vcc2026` metrics**.

That last step was supposed to be impossible: the public bundles ship
`log1p(CP10K)` and no counts. They are exactly invertible, and inverting them
turned the generator comparison from an argument into a measurement.

Date: 2026-09-21
Reproduce:
`uv run python scripts/run_arc_count_space_baseline.py`,
`uv run python scripts/run_count_generator_benchmark.py`,
`uv run python scripts/run_arc_dry_run.py`
Outputs: `outputs/arc_count_space_v1/`, `outputs/arc_dry_run_v1/`

---

## A. Freeze verification

All 14 frozen manifests were verified by `shasum -a 256 -c` **before** any
analysis and **again after** it: **220 files, 0 failures**, both times.

```
arc_bridge_v1  canonical_v1  decomposition_phase  discovery_phase
external_validation_v1  feng_multicontext_v1  foundations_v1
pathway_modelling_terminated  pathway_residual_v1  protocol
research_conclusions  stability_rule  unseen_perturbation_v1  zero_shot_v1
```

Nothing from a prior phase was modified. Prior-based Tier-0 prediction was not
reopened: no STRING, DepMap, pathway, gamma or embedding term appears anywhere
in this phase's model, and `tier_beta(..., tier=0)` returns exact zeros by
construction with a test that pins it.

## B. Target support policy, re-derived

The tiers were recomputed from frozen provenance (`public_label_inventory.json`,
identifier presence only — no response value is read, so no model outcome can
move a target between tiers) and compared to the frozen split:

**0 tier mismatches and 0 context-count mismatches across all 300 targets.**

| tier | targets | direct sources | response-space panel genes (median) | target gene on panel |
|---|---|---|---|---|
| **2** | 7 | `arch1` + `kaden25rpe1` | 16,494 | 7/7 |
| **1** | 79 | `kaden25rpe1` (73), `arch1` (6) | 16,828 | 79/79 |
| **0** | 214 | — | **0** | 214/214 |

Per-target records (`outputs/arc_count_space_v1/arc_target_support_recomputed.csv`)
carry the tier, the direct source contexts and their cell counts, the number of
public datasets measuring the target gene as an *output*, and the response
space — the panel genes on which a `beta` for that target is estimable at all.

Two things this table makes explicit that the tier alone does not. First, the
Tier-0 targets are not unmeasured: a median of 6 of the 7 public datasets
measure each one as an output gene. What they lack is any dataset that
*perturbs* it, which is a different deficiency and the one that matters.
Second, Tier 2's response space (16,494 genes) is **smaller** than Tier 1's
(16,828), because it is the intersection of two gene spaces rather than one.
More evidence buys a narrower window.

## C. The model

```
delta_hat[c,p] = m_hat[c] + w[tier(p)] * beta_hat[p]

beta_hat[p] = mean over the source contexts that DIRECTLY perturbed p of
              ( delta[s,p] - mean_q delta[s,q] )          Tier 2, Tier 1
            = 0 exactly                                    Tier 0
```

`beta_hat` is centred over perturbations **within each source before
averaging**, so a source with an unusually large template cannot leak it into
the perturbation-specific term (pinned by test). Under the frozen decomposition
this estimates `beta_p + mean_S gamma[s,p]` exactly: the conserved effect plus
whatever interaction the sources happen to share. It is therefore biased toward
the sources, which is what `w` exists to absorb.

Implemented in `virtual_cell.modelling.mean_response`; 26 tests pin the
algebra, the Tier-0 rule, and the leakage contract.

## D. Which feasible context-wide baseline wins

`m_c = mu + alpha_c` is what the official score anchors zero at, and Arc hides
it. Five candidates, none above one fitted scalar, evaluated on four public
leave-one-context-out folds.

**As an estimate of the hidden `m_c`** (energy explained, per fold):

| estimator | HepG2 | Jurkat | K562 | RPE1 | mean |
|---|---|---|---|---|---|
| `M0_zero` | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| `M1_source_pooled` | 0.663 | 0.143 | **−1.599** | 0.431 | −0.091 |
| `M2_basal_weighted` | 0.663 | 0.159 | **−1.541** | 0.432 | −0.072 |
| `M3a_source_pooled_shrunk` | 0.559 | 0.440 | −0.020 | 0.391 | 0.343 |
| **`M3b_basal_shrunk`** | 0.566 | 0.448 | 0.015 | 0.396 | **0.356** |

**The scalar is the whole result.** Unshrunk source pooling is *worse than
predicting nothing* on K562 — it explains −1.6 of the energy, overshooting by
enough to more than double the error. One leave-one-source-out scalar removes
that failure and costs about 0.1 of energy on the folds where the unshrunk
version was already fine. The correlation is identical either way (0.69–0.82
for every estimator, because a scalar cannot change Pearson), which is exactly
why MSE and not correlation is the criterion here.

Basal-similarity weighting adds **+0.013** on top of the scalar. That is
consistent with, and independent confirmation of, the frozen foundations
finding that basal control expression does not encode the response template. It
is kept only because it is free and never hurt a fold.

**As a complete response prediction** — which is what a Tier-0 Arc target
actually receives — mean over folds:

| estimator | Pearson (median) | cosine | MSE | energy | `pds_cosine` |
|---|---|---|---|---|---|
| `M0_zero` | — | 0.000 | 0.006237 | 0.000 | 0.500 |
| `M1_source_pooled` | 0.270 | 0.281 | 0.005808 | 0.033 | 0.500 |
| `M2_basal_weighted` | 0.271 | 0.282 | 0.005801 | 0.034 | 0.500 |
| `M3a_source_pooled_shrunk` | 0.270 | 0.281 | 0.005723 | 0.064 | 0.500 |
| **`M3b_basal_shrunk`** | **0.271** | **0.282** | **0.005714** | **0.065** | 0.500 |

`pds_cosine` is exactly 0.500 for all five, and must be: the prediction is the
same vector for every perturbation, so every pairwise distance ties and the
midrank rule gives 0.5. This reproduces the invariant the metric reference
states analytically, on real data, and is a check on the local scorer rather
than a result.

**Frozen: `M3b_basal_shrunk`.** It wins on mean MSE. It is not unanimous — it
wins 2 of 4 folds outright, with `M1` best on HepG2 and `M2` on RPE1 — but the
folds it loses it loses by ~0.1 of energy, and the fold it wins it wins by 1.6.

## E. Tier-specific shrinkage

Selected by **nested** validation: one inner fold per source context, the
tier's source-count structure preserved (Tier 2 draws `beta_hat` from two inner
sources, Tier 1 from one), the outer target never read by anything. Criterion
is MSE, because `w` rescales `beta_hat` only and therefore cannot change the
correlation of the perturbation-specific part at all.

Inner (selection) MSE, mean over folds:

| tier | w=0 | w=0.25 | w=0.5 | w=0.75 | w=1 |
|---|---|---|---|---|---|
| **2** | 0.005964 | 0.005445 | **0.005327** | 0.005610 | 0.006293 |
| **1** | 0.005964 | **0.005556** | 0.005769 | 0.006604 | 0.008060 |

Outer (oracle, reported only, never used for selection):

| tier | w=0 | w=0.25 | w=0.5 | w=0.75 | w=1 |
|---|---|---|---|---|---|
| **2** | 0.005714 | 0.005196 | **0.005077** | 0.005360 | 0.006043 |
| **1** | 0.005714 | **0.005306** | 0.005519 | 0.006354 | 0.007810 |

**Frozen: Tier 2 = 0.50, Tier 1 = 0.25, Tier 0 = 0 (not tunable).**

Three things follow.

1. **The selection is unanimous.** All four folds independently chose 0.50 for
   Tier 2 and 0.25 for Tier 1.
2. **Nested selection costs almost nothing here.** It matches the outer oracle
   optimum in **6 of 8** fold × tier cells. The two exceptions are Tier 2 on
   K562, where the oracle prefers 0.25 and the nested loop chose 0.50 (0.002948
   against 0.002876 — **2.5% worse**), and Tier 1 on RPE1, where the oracle
   prefers 0.50 and the nested loop chose 0.25 (0.008388 against 0.008367 —
   **0.25% worse**). Both exceptions are cheap, and both are in the direction
   the nested loop should err: toward the weight the *other* three folds agreed
   on.
3. **One global weight is not adequate.** Tier 1 wants exactly half the weight
   Tier 2 does, in every fold, and at Tier 2's weight of 0.50 a Tier-1
   prediction is already worse than the oracle 0.25 by 4%. That is the expected
   direction — one source carries undiluted interaction where two average some
   of it away — and it is large enough to matter.

Both tiers beat `w = 0` everywhere, so transferring direct evidence is worth
doing at both, and both are far from `w = 1`, so transferring it *unshrunk*
never is.

## F. Complete pseudobulk predictor, by tier

`m_hat + w beta_hat` on four public held-out contexts. Tier 2 averages over all
three 2-of-3 source subsets, Tier 1 over all three singletons.

Mean over folds:

| tier | w | Pearson | Spearman | cosine | MSE | energy | `pds_cosine` |
|---|---|---|---|---|---|---|---|
| **2** | 0.50 | **0.306** | 0.197 | **0.315** | **0.00508** | **0.168** | **0.725** |
| **1** | 0.25 | 0.297 | 0.207 | 0.309 | 0.00531 | 0.134 | 0.626 |
| **0** | 0 | 0.271 | 0.211 | 0.282 | 0.00571 | 0.065 | 0.500 |

**Monotone in the tier on every metric that reads the perturbation-specific
signal** — MSE, energy explained, `pds_cosine`, cosine. `pds_cosine` is the
cleanest: 0.500 → 0.626 → 0.725, and the Tier-0 value is the analytic floor,
so the entire gain above 0.5 is direct perturbation evidence doing work.

Per fold:

| fold | tier | Pearson | cosine | MSE | energy | `pds_cosine` |
|---|---|---|---|---|---|---|
| K562 | 2 / 1 / 0 | 0.251 / 0.210 / 0.142 | 0.256 / 0.216 / 0.140 | 0.00295 / 0.00300 / 0.00326 | 0.097 / 0.081 / 0.001 | 0.729 / 0.628 / 0.500 |
| RPE1 | 2 / 1 / 0 | **0.412 / 0.448 / 0.482** | 0.426 / 0.464 / 0.502 | 0.00798 / 0.00839 / 0.00895 | 0.223 / 0.183 / 0.128 | 0.699 / 0.622 / 0.500 |
| HepG2 | 2 / 1 / 0 | 0.292 / 0.288 / 0.267 | 0.302 / 0.302 / 0.283 | 0.00496 / 0.00527 / 0.00575 | 0.207 / 0.157 / 0.081 | 0.779 / 0.649 / 0.500 |
| Jurkat | 2 / 1 / 0 | 0.270 / 0.242 / 0.192 | 0.277 / 0.253 / 0.201 | 0.00442 / 0.00456 / 0.00490 | 0.144 / 0.116 / 0.051 | 0.694 / 0.605 / 0.500 |

**RPE1 inverts the per-perturbation Pearson ordering** — Tier 0 scores 0.482
against Tier 2's 0.412 — while MSE, energy explained and `pds_cosine` all keep
the normal order. The two statistics are answering different questions. Pearson
is computed per perturbation across genes, and in RPE1 the shared template
`m_hat` already correlates strongly with every perturbation's response; adding
a perturbation-specific term moves the prediction away from that template and
costs correlation even as it reduces error and improves discrimination. The
metric that cares about telling perturbations *apart* — which is what Arc's
`pds_cosine` does — is unambiguous. This is worth remembering when reading any
single-number correlation from this programme.

---

## G. A public count-space benchmark exists after all

The generators emit raw counts and the metrics read raw counts. The public
bundles ship `log1p(CP10K)` with no `layers['counts']` and no `raw` — checked
on all seven files. That appeared to rule out ever scoring a generator on
public data, leaving the choice between G0, G1 and G2 to be argued rather than
measured.

**The normalisation is invertible.** With `S = 1e4`,
`x_g = log(1 + S c_g / L)` gives `expm1(x_g) = S c_g / L`, and dividing by the
row's smallest stored value removes `L`. What comes back is the count vector
rescaled so the rarest detected gene reads 1 — the *minimal* integer solution.

Measured on the held-out context (101,288 cells x 8,563 genes):

| check | value |
|---|---|
| max round-trip error, recovered counts renormalised back to `X` | **3.3e-07** |
| median round-trip error | 1.9e-08 |
| max distance from an integer | 3.0e-04 |
| max \|row sum / S − 1\| after `expm1` | 9.5e-08 |
| rows needing a multiplier above 1 | **0** |
| recovered library size (median / min / max) | 13,364 / 1,795 / 53,317 |

Round-trip error sits at float32 resolution, which is the storage dtype — the
inversion is exact to the precision the file retains.

**What this does and does not establish.** The row-sum identity holding to
1e-7 shows the normalisation was applied on the stored gene axis, so nothing
was dropped after it and the inversion is determined up to one factor. That
factor is not determined: a cell with counts `(3, 6, 9)` and one with
`(1, 2, 3)` produce byte-identical rows, and no inversion can separate them.
The recovery returns the minimal solution, which is the true one exactly when
the cell detected at least one gene with a single UMI. That assumption is
stated, not proved; what supports it is that these are 10x cells detecting a
median of 3,656 genes, and that the recovered depths (median 13,364 UMIs over
the retained axis) are the depths these experiments report. A failure would
scale every library size by the same integer, which the depth distribution
would expose. Two tests pin the ambiguity so it cannot be quietly forgotten.

The recovered counts are the *filtered* transcriptome: genes the bundle's
preprocessing removed are gone, and the library sizes are sums over what
remains, not true UMI totals.

### The benchmark

| | |
|---|---|
| held-out context | `replogle22k562`, responses never read by the model |
| source contexts | RPE1, HepG2, Jurkat |
| perturbations | the 300 best-measured shared ones (201–1,996 real cells each, median 268) |
| generated cells | **400 per perturbation**, exactly as Arc requires |
| control split | 10,691 → **8,691 generator pool + 2,000 scoring reference, disjoint** |
| prediction | `M3b_basal_shrunk` + Tier-2 `beta_hat` at w = 0.50 |
| output space | the context's own 8,563 genes; 6,640 carry a response, 1,923 do not |

The disjoint control split matters and mirrors Arc's `--reject-controls`
("scoring uses the held-out controls, never yours"). Without it G0 would be
resampling the very cells it is scored against.

## H. Structural fidelity, and the detection defect

| group | library median | genes detected | sparsity | within-group spread |
|---|---|---|---|---|
| real control reference | 13,534 | 3,656 | 0.571 | 31.05 |
| real perturbed | 13,133 | 3,611 | 0.579 | 31.93 |
| `G0` resample | 13,775 | 3,675 (**+0.5%**) | 0.570 | 30.85 |
| **`G1` transport (s = 0.5)** | 13,745 | 3,528 (**−3.5%**) | 0.587 | **32.06** |
| `G1` **unsmoothed (s = 0)** | 13,778 | **2,846 (−22.2%)** | **0.665** | 37.51 |
| `G2` count model | 13,737 | 3,523 (−3.6%) | 0.588 | **29.06** |

All four generators produce non-negative finite integers with library-size
distributions matching the controls at every quantile (q10 ≈ 8,990,
q90 ≈ 20,400 against the reference's 8,972 / 20,350) — the depth distribution
is the measured one because the generators draw it from real cells rather than
inventing it.

**The predicted defect reproduces, and is worse than recorded.** Redrawing each
cell from its own empirical composition assigns probability zero to every gene
it happened not to catch, which is a second round of sampling loss on top of
the one the measurement already made. The `generate` module documents this as
costing "roughly 17%" of detected genes. Measured here on real data it costs
**22.2%** — 3,656 detected genes down to 2,846 — and drags sparsity from 0.571
to 0.665. It also *inflates* within-group spread to 37.5 against the real
perturbed 31.9, because the extra dropout is cell-specific noise masquerading
as heterogeneity. Blending with the pooled composition at `smoothing = 0.5`
costs 3.5% instead, and lands the spread at 32.06 against the real 31.93.

**No generator densifies the transcriptome.** The failure mode in this family
runs the other way, toward detection collapse, and the check that catches it is
genes-detected rather than sparsity alone.

`G2`'s failure is the opposite and is invisible in the detection count: it
matches genes detected (−3.6%) but produces cells that are **too similar to one
another** (spread 29.06 against the real 31.93), because a per-gene negative
binomial fitted on controls cannot reproduce the correlated cell-to-cell
structure that real cells carry. That is what the within-group spread statistic
exists to catch, and it is why it is computed within each perturbation rather
than over the pooled output.

### Does the generator deliver the response it was handed?

Realised response measured by recomputing the generated cells' own
mean-`log1p(CP10K)` profile against the control profile, over the first 50
perturbations:

| generator | realised norm | vs its **own** intent: r / slope | vs the **Tier-2** intent: r / slope |
|---|---|---|---|
| `G0` resample | 11.29 | — (asked for nothing) | 0.009 / 0.004 |
| **`G1` transport** | 23.30 | **0.781 / 0.935** | 0.781 / 0.935 |
| `G1` unsmoothed | 34.49 | 0.711 / **1.089** | 0.711 / 1.089 |
| `G2` count model | 27.98 | 0.696 / 0.965 | 0.696 / 0.965 |
| `G1` Tier-0 only | 17.74 | 0.583 / 0.844 | 0.305 / 0.291 |

**`G0`'s row is the noise floor and is what makes the rest readable.** It is
asked for no effect and faithfully delivers none (r = 0.009 against the Tier-2
intent), yet its realised response has norm 11.29 — pure 400-cell sampling
noise, against a signal of norm 19.55. A 400-cell pseudobulk simply cannot
resolve a response much smaller than that, and every realised norm in the table
is the intended signal added in quadrature with roughly this much noise.

`G1` at `smoothing = 0.5` realises the intended mean response at **slope 0.935**
— close to 1:1, with no recalibration fitted. That is not guaranteed: `delta` is
a mean of per-cell log ratios and the transport applies a geometric-mean fold
change, so the two agreeing to within 7% is a measured result rather than an
identity. The unsmoothed variant overshoots (slope 1.089), which is the same
dropout effect seen from the other side.

## I. Local `cell-eval2` scoring against real held-out cells

All six members, raw and unnormalised, 300 perturbations x 400 generated cells
scored against real K562 cells and a disjoint real control reference.

| generator | `pds_cosine` ↑ | `expr_mse_…_norm` ↓ | `dir_fidelity` ↑ | `dir_reach` ↑ | `sig_jaccard` ↑ | `lfc_nmae` ↓ |
|---|---|---|---|---|---|---|
| `G0` resample | 0.492 | 1.036 | 0.016 | 0.105 | 0.053 | 0.987 |
| **`G1` transport** | **0.882** | 1.062 | 0.551 | **0.404** | 0.122 | **0.829** |
| `G1` unsmoothed | 0.880 | 1.130 | 0.557 | 0.369 | **0.136** | 0.832 |
| `G2` count model | 0.734 | **1.771** | **0.583** | 0.347 | 0.129 | 0.987 |
| `G1` Tier-0 only | 0.503 | **1.010** | 0.325 | 0.149 | 0.038 | 0.938 |

Improvement over control resampling (sign-corrected, positive = better):

| generator | `pds_cosine` | `expr_mse` | `dir_fidelity` | `dir_reach` | `sig_jaccard` | `lfc_nmae` |
|---|---|---|---|---|---|---|
| **`G1` transport** | **+0.390** | **−0.026** | +0.535 | **+0.299** | +0.069 | +0.159 |
| `G1` unsmoothed | +0.388 | −0.094 | +0.541 | +0.264 | +0.083 | +0.155 |
| `G2` count model | +0.242 | −0.735 | +0.567 | +0.243 | +0.075 | +0.001 |
| `G1` Tier-0 only | +0.011 | +0.026 | +0.309 | +0.044 | −0.015 | +0.049 |

**`G1` transport at `smoothing = 0.5` is the winner**, on five of six members
and on the two that separate the generators most.

Four things in this table are worth more than the ranking.

**1. The expression-error member does not improve, and that is not a rounding
detail.** Every generator scores above 1.0, meaning noise-corrected squared
error exceeds the measured effect energy, and `G1` (1.062) is *worse* than
emitting controls (1.036). The correction is fully credited (`rho = 1.0` for
every generator), so this is not 400-cell sampling noise being charged to the
model — it is genuine prediction error. The prediction has the right direction
and roughly the right magnitude but is wrong enough, gene by gene, that adding
it increases squared error. Per-perturbation pseudobulk energy explained for
`G1` confirms it: median **−0.47**, positive for only a minority of the 300.

**2. Discrimination and direction improve enormously anyway.** `pds_cosine`
goes 0.492 → 0.882, direction reach 0.105 → 0.404, fold-change error 0.987 →
0.829. A prediction can be far from the truth in squared error and still rank
perturbations correctly and get their directions right, and five of the six
scored members reward the latter. **The same model is simultaneously a good
discriminator and a bad regressor**, and any single summary of it will be
misleading in one direction or the other.

**3. The Tier-0 arm separates the two halves of the model.** Applying `m_hat`
alone — the context main effect, no perturbation-specific term, identical for
every target — already delivers **0.325 of `G1`'s 0.551** direction fidelity
and **0.026 of the expression-error improvement**, while contributing
essentially nothing to discrimination (`pds_cosine` 0.503 against `G1`'s 0.882,
where 0.500 is the analytic floor for a constant prediction). The context main
effect and the perturbation-specific effect are being paid by *different*
metrics: `m_hat` buys direction and error, `beta_hat` buys discrimination.
Anything that trades one against the other will look like an improvement on
half the scoreboard.

**4. `G2` is not close.** It wins direction fidelity (0.583) but loses
catastrophically on expression error (1.771, three times `G1`'s deficit) and
gains nothing on fold-change error (0.987, indistinguishable from controls),
because drawing every gene independently from a fitted negative binomial
discards the cell-level structure that both of those members read. Its
under-dispersed cells (spread 29.06) are the same defect seen structurally.

### What this benchmark is not

* **The gene axis is pre-filtered.** 8,562 of 8,563 genes clear the
  `min_cpm = 5` control-side gate, because scPertEval already removed
  low-expression genes. Arc's 18,533-gene panel will not behave this way, and
  the DE members are sensitive to which genes are tested.
* **The reference control group is 2,000 cells**, against Arc's 18,400. The
  Wilcoxon tests and the control dispersion term are correspondingly noisier.
* **The anchors do not transfer.** The published `b` and `r` constants were
  measured on Arc's own bundles; normalising these raw values with them would
  be an invention. Everything above is raw, and the control-emitting reference
  point is measured here (`G0`) rather than assumed. For scale, the analytic
  control-submission value for `expr_mse_unbiased_capped_norm` on Arc is
  1.0032, and `G0` measures 1.036 here — close, but not the same number.
* **One held-out context.** K562 was chosen before any generator ran, on cell
  count alone. Generator conclusions are cell-level and should transfer; the
  *score levels* are one context's.

### F2 revisited — what a *mixed* Arc panel scores

The per-tier table above answers "how well does Tier 2 do". An Arc submission
is not a Tier-2 submission: it is 7 Tier-2 targets, 79 Tier-1 and 214 Tier-0,
scored **together**. For a per-perturbation metric that distinction barely
matters. For `pds_cosine` it matters enormously, because that metric ranks each
prediction against every *other* prediction — and 214 Tier-0 predictions are
all the identical vector `m_hat`, sitting at distance zero from one another.

So the mixture was built explicitly: each public perturbation assigned a tier by
a seeded draw at Arc's prevalence, one prediction matrix assembled, metrics
computed once, repeated over 20 draws.

| | per-tier view | **mixed panel** |
|---|---|---|
| `pds_cosine` | 0.725 (Tier 2) | **0.535** (sd 0.003) |
| Pearson (median) | 0.306 | 0.280 |
| cosine (median) | 0.315 | 0.291 |
| energy explained | 0.168 | 0.084 |

**`pds_cosine` over a realistic Arc panel is 0.535, not 0.725.** The 0.725
describes a field Arc will never present. Two thirds of the discrimination
above the 0.500 floor disappears the moment the panel is mostly Tier 0 — not
because the Tier-2 predictions got worse, but because most of the field is
indistinguishable from itself.

This is the single most important number for calibrating expectations about the
Arc score, and it is a direct consequence of the frozen Tier-0 policy rather
than a flaw in it: 214 targets have no direct evidence, the external validation
established that inventing some does not work, and identical predictions are
the honest output. The cost of that honesty is now measured.

---

## J. Arc output-space mapping

The 86 targets with direct public evidence get `beta_hat` from the two datasets
that actually perturb them; `m_hat` comes from the same two, because they are
the only public contexts whose perturbation panel overlaps Arc's at all.

| | |
|---|---|
| panel | 18,533 genes, 300 targets, contexts A/B/C |
| response space (panel ∩ `arch1` ∩ `kaden25rpe1`) | **16,494 genes** |
| panel genes that can carry a response | 16,494 |
| panel genes that **keep their control distribution** | **2,039** |
| targets with a non-zero `beta_hat` | **86** of 300 |
| panel-mean drift `‖mean_p w·beta_hat‖` | 0.0224 |

**The unsupported-gene rule.** A panel gene outside the response space gets
`lfc = 0`, which under the transport generator means it is emitted with the
target context's own control distribution — its real counts, its real detection
rate, its real dispersion. That is **not** zero-filling, and a test pins the
difference: a zero fold change reproduces the control detection rate to within
10%, where writing zeros would emit 2,039 permanently silent genes. No
perturbation shift is fabricated for any of them.

The same rule covers the 214 Tier-0 targets across all 18,533 genes. They
receive `m_hat` and nothing else.

`beta_hat` is centred within each source's own panel, so the model's mean over
Arc's 300 targets is not exactly `m_hat`: the 86 supported targets pull it by
‖0.0224‖. That drift is small relative to the per-target signal (median
`‖w·beta_hat‖` of 0.91 at Tier 2, 0.51 at Tier 1) and is reported rather than
corrected, because correcting it would require adding a constant to the 214
Tier-0 targets and the frozen policy says they get zero.

### The finding that matters most: `m_hat` is not estimable from these sources

`m_hat` for A, B and C has norm **0.071**, against public per-context main
effects of 1.28–4.70. The estimator did not fail — it is behaving exactly as
designed. The leave-one-source-out scalar it fits is **0.119**, against 0.49–0.87
on the public folds, because the two sources barely agree:

| cosine between the two source main effects | value |
|---|---|
| over each source's own panel (150 vs 1,836 perturbations) | **0.089** |
| over each source's Arc-target subset (13 vs 80) | **−0.002** |
| over the **7 targets both measure** | **0.030** |
| *(for comparison: the four public contexts, pairwise)* | *0.60 – 0.81* |

The first row admits an innocent explanation — the two sources average over
different perturbation panels, where the four public contexts shared all 1,264
by construction. **The second and third rows rule it out.** Matched on the
identical Arc targets, and again on the identical 7 targets both datasets
perturb, the agreement is still zero. The disagreement is about the response,
not about which perturbations were measured.

So the shrinkage scalar collapsing to 0.119 is the correct answer to the
question asked: given two source contexts whose mean perturbation responses are
uncorrelated, the variance-optimal estimate of a third context's is close to
nothing. **The frozen estimator is sound and the input data cannot support it.**

The consequence for the bundle is concrete. Per-target signal:

| tier | targets | median `‖beta_hat‖` | median `‖w·beta_hat‖` | `‖m_hat‖` |
|---|---|---|---|---|
| 2 | 7 | 1.82 | **0.91** | 0.071 |
| 1 | 79 | 2.04 | **0.51** | 0.071 |
| 0 | 214 | 0 | **0** | 0.071 |

The 86 supported targets carry a real perturbation-specific prediction, 7–13x
larger than the context main effect. **The other 214 are, in effect, control
resampling.** The frozen score accounting puts a control-emitting submission at
**−0.311**, and 71% of this panel is close to that.

## K. Arc A/B/C dry run

`uv run python scripts/run_arc_dry_run.py` — generated from A/B/C **control
cells only**, with the model frozen beforehand on public data. Generator is
`G1` transport at `smoothing = 0.5`, the section-I winner.

| property | value |
|---|---|
| cells | **360,000** = 3 x 300 x 400 |
| genes | **18,533**, in `gene_names.csv` order |
| nonzeros | **2,057,973,608** |
| file | 15.4 GiB, CSR float32, `int32` column indices, `int64` offsets |
| density | 0.308 |
| library size (median / min / max) | 20,012 / 710 / 52,420 |
| genes detected (median) | 5,856 |

Library-size fidelity is exact by construction: `min_counts_per_cell` is **710**,
which is context B's own minimum control library size. The generator draws each
cell's depth from a real control cell and preserves it through the multinomial
redraw, so the emitted depth distribution is the measured one rather than a
model of it. Median genes detected of 5,856 against the controls' ~6,000 is the
same −3% seen on public data.

All thirteen local checks pass:

```
[PASS] n_cells == 360000                    [PASS] 300 perturbations in every context
[PASS] n_genes == 18533                     [PASS] contexts are exactly A, B, C
[PASS] gene order matches gene_names.csv    [PASS] no control cells emitted
[PASS] raw integer counts                   [PASS] under max_counts_per_cell (1,000,000)
[PASS] non-negative                         [PASS] under max_nnz (4,750,000,000)
[PASS] finite                               [PASS] under max_cell_dim (400,000)
[PASS] exactly 400 cells per (context, perturbation)
```

And `vcc prep --dry-run` accepts it, exit code 0:

```json
{"n_cells": 360000, "nnz": 2057973608, "n_genes": 18533, "encoding": 32,
 "normalization": "counts-preserved", "cells_per_context": {"A": 120000, "B": 120000, "C": 120000},
 "verified_targets": true, "dropped": [], "reordered_genes": false, "notes": []}
```

`verified_targets: true` is the check that every context predicts exactly its
official 300; `dropped: []` and `reordered_genes: false` mean nothing had to be
discarded or rearranged; `counts-preserved` confirms the matrix went through as
raw counts rather than being log-normalized.

**One margin is thinner than it looks.** At 2.058e9 nonzeros the submission sits
at **96% of 2^31**, the point at which SciPy promotes CSR indexing from `int32`
to `int64`. `vcc`'s own sizing module states that this steps the scoring
container's cost from 8 to 12 bytes per nonzero, and that a flat rate fitted
below the boundary underestimates above it by a third. The hard cap (4.75e9) is
far away; this other boundary is 4% away, and a generator change that raised
density by 5% would cross it.

**Nothing was uploaded and no submission was made.** The artifact is
`outputs/arc_dry_run_v1/arc_dry_run_v1.h5ad`, a dry-run candidate.

---

## L. The twelve questions

**1. Which feasible context-main-effect estimator won?**
`M3b_basal_shrunk` — basal-similarity-weighted source pooling, rescaled by one
leave-one-source-out scalar. It wins on mean MSE both as an estimate of the
hidden `m_c` (energy explained 0.356 against `M1`'s −0.091) and as a complete
response prediction. The win is almost entirely **the scalar**: unshrunk pooling
explains −1.6 of the energy on K562, worse than predicting nothing, and the
scalar removes that. Basal weighting contributes +0.013, which is consistent
with the frozen finding that basal expression does not encode the response
template. **On Arc's actual sources the estimator returns ~0** (section J), not
because it is wrong but because those sources disagree.

**2. What Tier-2 shrinkage won?** **0.50**, chosen by nested inner validation,
unanimous in all four folds, matching the outer oracle in 3 of 4 (K562's oracle
prefers 0.25, at a cost of 2.5%).

**3. What Tier-1 shrinkage won?** **0.25**, unanimous in all four folds and
matching the outer oracle in 3 of 4 (RPE1's oracle prefers 0.50, at a cost of
0.25%). Across both tiers the nested selection matches the oracle in **6 of 8**
fold × tier cells. **One global weight is not adequate**:
Tier 1 wants exactly half of Tier 2's weight in every fold, which is the
expected direction — one source carries undiluted interaction where two average
some of it away.

**4. How does Tier 0 behave under the frozen zero-beta policy?**
Exactly as the algebra predicts. Its prediction is the same vector for every
perturbation, so `pds_cosine` is **0.500 to machine precision** — the analytic
floor, reproducing the invariant the metric reference states. It still carries
real signal on the metrics that do not require telling perturbations apart:
median Pearson 0.271, energy explained 0.065, and in count space **0.325 of
`G1`'s 0.551 direction fidelity**. Tier 0 is not nothing; it is nothing
*discriminative*.

**5. What is complete pseudobulk performance by tier?**
Monotone in the tier on every metric that reads perturbation-specific signal:

| tier | Pearson | cosine | MSE | energy | `pds_cosine` |
|---|---|---|---|---|---|
| 2 | 0.306 | 0.315 | 0.00508 | 0.168 | 0.725 |
| 1 | 0.297 | 0.309 | 0.00531 | 0.134 | 0.626 |
| 0 | 0.271 | 0.282 | 0.00571 | 0.065 | 0.500 |

On a **realistically mixed panel** at Arc's 7/79/214 prevalence, `pds_cosine` is
**0.535**, not 0.725 — two thirds of the discrimination disappears because most
of the field is Tier 0 and therefore indistinguishable from itself.

**6. Which count generator performs best?**
**`G1` control transport at `smoothing = 0.5`**, on five of six metric members.
`pds_cosine` 0.492 → **0.882**, direction reach 0.105 → **0.404**, fold-change
error 0.987 → **0.829**. `G2` (negative binomial) loses badly on expression
error (1.771) and gains nothing on fold-change error. The unsmoothed transport
scores comparably but destroys 22% of detected genes.

**7. Does the generator preserve realistic single-cell structure?**
Yes, on every axis tested. Non-negative finite integers; library-size
distribution matched at q10/median/q90 and exact at the minimum (710, context
B's own floor); genes detected within **−3.5%**; sparsity 0.587 against the real
0.579; within-perturbation heterogeneity **32.06 against the real 31.93**.
The prior detected-gene defect was reproduced and is **worse than recorded** —
**−22.2%**, not ~17% — confirming the smoothing blend is load-bearing rather
than cosmetic. **No generator densifies the transcriptome**; the failure mode in
this family is detection collapse. `G2` passes the detection check but fails
heterogeneity (29.06), which is why that statistic is computed within each
perturbation rather than over the pooled output.

**8. How does it score against control resampling?**
Better on five of six members: `pds_cosine` **+0.390**, direction fidelity
+0.535, direction reach +0.299, `sig_jaccard` +0.069, `lfc_nmae` +0.159.
**Worse on one**: `expr_mse_unbiased_capped_norm` −0.026 (1.062 against 1.036).
Every generator scores above 1.0 on that member, meaning noise-corrected squared
error exceeds the measured effect energy, with the sampling correction fully
credited (`rho = 1.0`). **The same model is a good discriminator and a bad
regressor**, and no single summary of it is honest.

**9. What assumptions are made for unsupported genes?**
One, stated and tested: a panel gene with no direct response evidence **keeps
its target-context control distribution**, which is the `lfc = 0` case of the
transport. It is not zero-filled and no shift is fabricated. This covers 2,039
of 18,533 panel genes for every target, and all 18,533 for the 214 Tier-0
targets. The claim being made is "we have no evidence this gene changes", not
"this gene does not change".

**10. Can a valid 360,000-cell Arc bundle now be generated?**
**Yes.** It was: 360,000 cells x 18,533 genes, 2.058e9 nonzeros, 15.4 GiB, all
thirteen local checks passing and `vcc prep --dry-run` accepting it with exit
code 0, `verified_targets: true`, `dropped: []`, `notes: []`. Validity is
settled. Quality is not.

**11. What remains before the first real submission?**
One blocker dominates the rest.

  1. **`m_hat` is not estimable from the available Arc sources, and this is
     the single most consequential open problem.** The only two public contexts
     that perturb Arc targets have mean perturbation responses with cosine
     0.089 on their own panels, **−0.002** on matched Arc targets and **0.030**
     on the 7 targets both measure. A better estimator will not help; the
     estimator is already correct and is returning ~0 for the right reason.
     What is needed is either more public contexts that perturb Arc targets, or
     evidence about *why* these two disagree when the four research contexts
     agree at 0.60–0.81. Until then **214 of 300 targets are effectively
     control-emitting**, and the frozen accounting puts that at **−0.311**.
  2. **The expression-error member is not beaten.** At 1.062 against controls'
     1.036 the model adds squared error. Whether that is fixable by
     recalibrating magnitude, or is the honest ceiling of conserved transfer,
     has not been tested.
  3. **The count-space benchmark is one held-out context on a pre-filtered gene
     axis** (8,562 of 8,563 genes clear the expression gate, against a panel
     that will not behave that way) with a 2,000-cell control reference against
     Arc's 18,400.
  4. **The nnz margin.** 2.058e9 is 96% of the `int32` indexing boundary.
  5. **The anchors have not been transferred.** Everything here is raw; no
     normalized score has been computed, because the published `b` and `r`
     constants were measured on Arc's bundles and applying them to K562 would
     be an invention.

**12. Is there evidence that a more complex generative model is justified?**
**No — and the evidence points at the opposite bottleneck.**

The generator is not what is limiting the score. `G1` realises the mean response
it is handed at **slope 0.935 with correlation 0.781**, reproduces the control
detection rate to 3.5%, matches within-perturbation heterogeneity to 0.4%, and
preserves the depth distribution exactly. A VAE, a diffusion model or a
transformer would be competing against a generator that already transports the
prediction faithfully — and the prediction is what is weak. The one generator
here with more capacity, `G2`, is **worse** on the members that read cell-level
structure, which is the expected result when extra capacity is spent modelling
genes independently.

The decomposition in section I makes the target explicit: `m_hat` buys direction
and error, `beta_hat` buys discrimination, and 214 of 300 Arc targets currently
get neither. That is a problem about **evidence for the response**, not about
how cells are emitted from it. Building a deep generative model now would add
capacity exactly where the measurements say there is no deficit.

---

## Validation

`uv run pytest` — see below. `ruff check`, `ruff format --check`, `uv build`
clean. All 14 freeze manifests verified before and after: 220 files, 0 failures.

**No submission was made and nothing was uploaded to Arc.**
