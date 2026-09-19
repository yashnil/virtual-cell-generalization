# Independent four-context decomposition — predeclared robustness battery

**This is a sensitivity analysis of the independent four-context decomposition.
It is NOT a Molina & Zhang reproduction.** No option was tuned toward any
published percentage; agreement with their paper is not a criterion anywhere in
this document.

Date: 2026-09-19 · Battery runtime 57.8 min · peak RSS 23.2 GB of 68.7 GB
Reproduce: `uv run python scripts/run_four_context_sensitivity.py`
then `uv run python scripts/plot_four_context_sensitivity.py`

## Canonical v1 was frozen first

`data/provenance/scperteval/canonical_v1_freeze.txt` pinned SHA-256 digests of
the v1 report, the frozen design files, the four derived arrays, `summary.json`,
and the two code files that produced them — **before** the battery ran.
Re-verified afterwards: **12/12 match**, and the earlier protocol freeze is
**4/4**. Canonical v1 was not modified, and the decomposition protocol was not
altered in response to any result.

**Independent cross-check.** The battery re-implements the pipeline in a
separate module (`virtual_cell.analysis.robustness`) rather than reusing the
frozen one. At the canonical settings it reproduces canonical v1 to 0.01 pp
(template 20.27, β 30.07, γ 21.04, noise 28.62; reproducibility β 80.8 %,
γ 49.5 %) using 20 resamples instead of 50 — two independent code paths agreeing.

Design unchanged throughout: the frozen **1,264 perturbations × 6,640 genes**.

---

## The three conclusions under test

1. conserved **β** is substantial
2. reproducible context-specific **γ** is substantial
3. raw **γ** is strongly noise-inflated

## Summary across all 21 variant × feature-space combinations

| quantity | min | max | spread |
|---|---:|---:|---:|
| template share | 19.09 % | 25.26 % | 6.2 pp |
| **β share** | **27.98 %** | **30.79 %** | **2.8 pp** |
| **γ share (corrected)** | **20.55 %** | **22.80 %** | **2.2 pp** |
| noise share | 21.72 % | 32.39 % | 10.7 pp |
| β reproducibility | 77.4 % | 85.0 % | 7.6 pp |
| γ reproducibility | 45.6 % | 57.7 % | 12.1 pp |
| **γ share (uncorrected)** | **38.51 %** | **45.05 %** | — |

**β never falls below 28 %. Corrected γ never falls below 20.5 %. Uncorrected γ
is roughly double corrected γ in every single variant.** All three conclusions
hold everywhere.

Mathematical invariants across every variant: reconstruction ≤ 4.44e-16, SS
partition ≤ 5.15e-16, zero-sum ≤ 1.07e-13.

---

## A. Independent control split

Canonical reuses one full-context control mean in both perturbation halves. Here
the control cells are themselves divided into two disjoint equal halves and each
perturbation half is referenced against its own. All other choices identical
(seed 42, 20 resamples, all 6,640 genes).

| | template | β | γ | noise |
|---|---:|---:|---:|---:|
| shared control (canonical) | 20.267 % | 30.068 % | 21.045 % | 28.621 % |
| **independent control split** | **20.021 %** | **30.068 %** | **21.045 %** | **28.867 %** |
| **change** | **−0.246 pp** | **0.000 pp** | **0.000 pp** | **+0.246 pp** |

Component reproducibility:

| | μ | α | β | γ |
|---|---:|---:|---:|---:|
| shared control | 99.961 % | 99.754 % | 80.791 % | 49.525 % |
| independent split | 99.484 % | 97.296 % | 80.791 % | 49.525 % |
| change | −0.48 pp | **−2.46 pp** | **0.000 pp** | **0.000 pp** |

### Which components absorb control-estimation error — and why it must be μ and α

**β and γ are exactly unchanged, to 0.000 pp.** This is not luck; it is
algebraically necessary. The control profile enters as a term that depends on the
context but **not** on the perturbation:

```
delta[c,p] = pert_mean[c,p] - ctrl[c]
```

In the balanced two-way decomposition:

* `mu = mean_{c,p} pert_mean - mean_c ctrl`, so μ carries `-mean_c ctrl`.
* `alpha_c = mean_p delta[c,p] - mu = mean_p pert_mean[c,p] - ctrl[c] - mu`,
  so α carries the context-specific control deviation.
* `beta_p = mean_c delta[c,p] - mu = mean_c pert_mean[c,p] - mean_c ctrl - mu`.
  Substituting μ, the two `mean_c ctrl` terms **cancel exactly** — β is free of
  the control reference.
* `gamma[c,p] = delta[c,p] - mu - alpha_c - beta_p`. The `-ctrl[c]` in `delta`
  and the `-ctrl[c]` inside `alpha_c` **cancel exactly** — γ is likewise free of
  it.

So the control estimate is a per-context constant vector, and a per-context
constant lives entirely in the μ + α subspace, which is orthogonal to both β
(a perturbation contrast averaged over contexts) and γ (a context × perturbation
interaction). **Control-estimation error can only ever inflate the template.**

The empirical numbers confirm this exactly: the whole 0.246 pp moved from
template into noise, α's reproducibility dropped by 2.46 pp (it now carries real
control-sampling error), and β and γ did not move by a single digit.

**Consequence for the canonical result: reusing the control mean cannot have
inflated β or γ.** The concern that motivated this check is real for the
template but provably absent for the two components the science rests on.

---

## B. Controlled cell-depth experiment

Design: the **same 643 (context, perturbation) pairs** — every pair with ≥ 200
cells, so all four depths are feasible — re-estimated at n = 15, 30, 50, 100
cells per half. For each pair and depth, 25 repeats draw 2n cells without
replacement, split into two disjoint halves of n, and correlate the two half
responses. Perturbation identity and context are held fixed, so **depth is the
only thing that varies**. Nothing conditions on the observed ‖δ‖, which is
itself depth-biased.

Pairs per context: K562 309, RPE1 144, HepG2 43, Jurkat 147.

### Median reliability by depth (bootstrap 95 % CI over pairs)

| n per half | median reliability | 95 % CI | mean |
|---:|---:|---|---:|
| 15 | 0.0966 | [0.0812, 0.1106] | 0.170 |
| 30 | 0.1762 | [0.1490, 0.2052] | 0.254 |
| 50 | 0.2626 | [0.2256, 0.2908] | 0.327 |
| 100 | **0.4188** | [0.3733, 0.4551] | 0.432 |

**The confidence intervals do not overlap between adjacent depths.** Reliability
more than quadruples from n = 15 to n = 100.

### Within-pair change

| transition | median Δ | mean Δ | fraction improved | fraction worse |
|---|---:|---:|---:|---:|
| 15 → 30 | +0.0787 | +0.0847 | 94.71 % | 5.29 % |
| 30 → 50 | +0.0776 | +0.0725 | 96.42 % | 3.58 % |
| 50 → 100 | +0.1146 | +0.1049 | 98.76 % | 1.24 % |
| **15 → 100** | **+0.3010** | **+0.2621** | **99.53 %** | **0.47 %** |

**90.8 % of pairs (584/643) increase at every single step**, and 99.53 % are
higher at n = 100 than at n = 15.

### Every context independently

| context | pairs | median Δ (15→100) | fraction improved |
|---|---:|---:|---:|
| HepG2 | 43 | +0.3700 | 100.00 % |
| K562 | 309 | +0.3322 | 100.00 % |
| Jurkat | 147 | +0.2836 | 98.64 % |
| RPE1 | 144 | +0.2402 | 99.31 % |

Within-repeat sd also falls monotonically with depth (0.0367 → 0.0327 → 0.0278 →
0.0183), so deeper estimates are both higher and more stable.

**Conclusion: holding perturbation identity and context fixed, reliability
increases with cell depth.** This settles the question the v1 report raised. The
weak *marginal* correlation reported in v1 §17.1 was indeed an artifact of the
negative depth–magnitude coupling; this design removes the confound entirely
rather than conditioning on the biased ‖δ‖.

---

## C. Feature-space sensitivity

HVGs are ranked by a **single global rule using control cells only**: for each
context compute the per-gene variance of `X` across that context's `control`
cells on the 6,640 shared genes, then score each gene by the **mean of the four
per-context variances** and rank descending. One ranking is applied to all four
contexts. **No perturbation response is consulted at any point**, so the feature
space cannot leak a response.

(shared control, mean(log), seed 42)

| feature set | genes | total SS | template | β | γ | noise | β repro | γ repro |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| all shared | 6,640 | 41.417 | 20.27 % | **30.07 %** | **21.04 %** | 28.62 % | 80.8 % | 49.5 % |
| top HVG | 4,000 | 35.416 | 22.38 % | **29.69 %** | **21.79 %** | 26.13 % | 82.0 % | 52.7 % |
| top HVG | 2,000 | 26.270 | 25.25 % | **30.79 %** | **22.22 %** | 21.74 % | 85.0 % | 57.7 % |

β moves by 1.1 pp and γ by 1.2 pp across a 3.3-fold change in feature-space size.
Noise falls from 28.6 % to 21.7 % and both reproducibilities rise, exactly as
expected when low-variance, noise-dominated genes are dropped — HVG restriction
concentrates signal, it does not create or destroy either component.

---

## D. Aggregation-order sensitivity

Canonical takes `mean_cells(log1p(CP10K))`. The variant recovers CP10K with
`expm1`, averages cells, then applies `log1p`: `log1p(mean_cells(CP10K))`.
Neither is treated as ground truth. By Jensen's inequality the second is ≥ the
first, so they are genuinely different estimators (asserted in the tests).

(shared control, seed 42)

| feature set | | template | β | γ | noise | β repro | γ repro |
|---|---|---:|---:|---:|---:|---:|---:|
| all 6,640 | mean(log) | 20.27 % | **30.07 %** | **21.04 %** | 28.62 % | 80.8 % | 49.5 % |
| all 6,640 | log(mean) | 19.09 % | **27.98 %** | **20.55 %** | 32.39 % | 77.4 % | 45.6 % |
| 4,000 HVG | mean(log) | 22.38 % | 29.69 % | 21.79 % | 26.13 % | 82.0 % | 52.7 % |
| 4,000 HVG | log(mean) | 21.42 % | 28.54 % | 21.77 % | 28.27 % | 79.9 % | 50.3 % |
| 2,000 HVG | mean(log) | 25.25 % | 30.79 % | 22.22 % | 21.74 % | 85.0 % | 57.7 % |
| 2,000 HVG | log(mean) | 24.81 % | 30.08 % | 22.80 % | 22.31 % | 84.1 % | 57.0 % |

`log(mean)` is the noisier estimator — it shifts 2.1 pp of β and 0.5 pp of γ into
noise on the full gene set, because averaging on the linear scale gives high-count
cells more leverage. The gap narrows as the feature set is restricted to
well-expressed genes (β differs by only 0.7 pp at 2,000 HVGs). **Both conclusions
survive under either aggregation order**, and the direction of the difference is
the one Jensen predicts rather than anything anomalous.

---

## E. Seed stability

Five independent seeds, shared control, mean(log), all 6,640 genes, 20 resamples
each.

| seed | template | β | γ | noise |
|---|---:|---:|---:|---:|
| 7 | 20.2647 % | 30.0610 % | 21.0448 % | 28.6294 % |
| 42 | 20.2666 % | 30.0683 % | 21.0445 % | 28.6206 % |
| 101 | 20.2661 % | 30.0682 % | 21.0502 % | 28.6155 % |
| 2024 | 20.2644 % | 30.0652 % | 21.0479 % | 28.6224 % |
| 31337 | 20.2747 % | 30.0668 % | 21.0489 % | 28.6095 % |
| **range** | **0.0104 pp** | **0.0072 pp** | **0.0057 pp** | **0.0199 pp** |

Seed choice is irrelevant at the third decimal place. The split-half estimator is
effectively deterministic at this sample size.

---

## Figures

In `outputs/four_context_sensitivity/` (git-ignored, regenerable):

| figure | content |
|---|---|
| `variant_comparison.png` | stacked component shares for all 7 variants × 3 feature spaces |
| `beta_gamma_stability.png` | every β and γ estimate, showing neither approaches zero |
| `depth_experiment.png` | median reliability vs depth with bootstrap CIs, per-context curves, within-pair trajectories |
| `seed_stability.png` | all five seeds per component with ranges annotated |

Data: `variants.json`, `variant_table.csv`, `depth_experiment.csv`,
`depth_within_pair.csv`, `hvg_ranking.csv`.

---

## Did the canonical conclusions survive?

| conclusion | verdict | evidence |
|---|---|---|
| **1. conserved β is substantial** | **SURVIVED** | 27.98–30.79 % across all 21 combinations; reproducibility 77.4–85.0 %; *exactly* unchanged by the control split |
| **2. reproducible γ is substantial** | **SURVIVED** | 20.55–22.80 % after noise correction in every variant; reproducibility 45.6–57.7 %; *exactly* unchanged by the control split |
| **3. raw γ is strongly noise-inflated** | **SURVIVED** | uncorrected γ 38.51–45.05 % against corrected 20.55–22.80 % — roughly halved in every variant without exception |

No variant reverses, weakens to degeneracy, or qualitatively changes any of the
three. The largest single effect anywhere in the battery is the aggregation-order
shift of 2.1 pp in β, which is an order of magnitude too small to threaten a
conclusion.

## Gate criteria

| criterion | status | evidence |
|---|---|---|
| mathematical invariants remain valid | **PASS** | reconstruction ≤ 4.44e-16, SS partition ≤ 5.15e-16, zero-sum ≤ 1.07e-13, across every variant |
| β clearly non-degenerate across reasonable variants | **PASS** | never below 27.98 % |
| γ clearly non-degenerate after noise correction | **PASS** | never below 20.55 % |
| conclusion does not depend on one preprocessing choice | **PASS** | holds across control scheme, feature space (3.3× size range), aggregation order, and 5 seeds |
| controlled subsampling confirms reliability improves with depth | **PASS** | median 0.097 → 0.419, non-overlapping CIs, 99.53 % of fixed pairs improve, all four contexts |
| no single reasonable variant reverses the interpretation | **PASS** | none did |

**All six gate criteria now pass.** Combined with the v1 criteria (invariants,
β/γ quantified, split-half separates signal from noise), the independent
four-context decomposition gate is **met**.

This remains an independent re-derivation. It does not become a Molina & Zhang
reproduction, and their percentages were not a target at any point.

## Methodological concerns worth recording

1. **Absolute reliability is low for most perturbations.** Even at n = 100 the
   median split-half r is 0.42, and at the observed median depth (55–139 cells)
   most perturbations sit well below that. The decomposition is trustworthy in
   aggregate — the component shares are stable to ~0.01 pp — but an individual
   (context, perturbation) response is often poorly determined. **Any future
   per-perturbation modelling must carry a reliability weight or filter, and must
   not treat all 1,264 × 4 responses as equally informative.**
2. **γ's reproducibility of ~50 % is the binding constraint on anything that
   predicts γ.** Half of the raw interaction signal is measurement noise, so the
   achievable ceiling for a γ predictor evaluated against raw γ is bounded well
   below 1 regardless of model quality. Evaluation must be against the
   noise-corrected or split-half-reliable portion, or the ceiling must be
   reported alongside.
3. **HepG2 remains the weak context**: only 43 pairs had ≥ 200 cells, and its
   median depth is the lowest of the four. Its per-context estimates carry the
   widest uncertainty even though its *relative* depth response is the strongest.
4. **‖δ‖ is depth-biased upward** (v1 §17.1, confirmed here). It must not be used
   as an effect-size filter or ranking key without a depth correction.
5. Not a defect, but worth stating: the battery used 20 resamples where canonical
   used 50. Given the seed ranges of ~0.01 pp this is immaterial, and the 20-
   resample canonical-settings run reproduced the 50-resample canonical result to
   0.01 pp.

## Standing statement

Independent four-context decomposition on standardized public scPertEval data.
**Not** an exact Molina & Zhang reproduction. Arc contexts A/B/C were not used
anywhere: not trained on, not intersected with the response space, not used to
choose any preprocessing, and no identity inference was attempted. No predictive
model — Ridge, MLP, γ predictor, D predictor, transferability classifier, or
generative model — was built.
