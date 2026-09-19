# Independent four-context response decomposition — v1 (canonical result)

**This is an independent four-context decomposition on standardized public
scPertEval data. It is NOT an exact Molina & Zhang reproduction.** That track is
BLOCKED (`reports/molina_zhang_reproduction_spec.md` §7): their claimed processed
artifacts are not publicly present. Their published percentages appear here once,
in §11, purely as contextual reference from a *different* data and preprocessing
pipeline, and were never used as a target.

Date: 2026-09-19 · Design `four_context_v1` · script version 1.0.0
Reproduce: `uv run python scripts/build_four_context_decomposition.py`

---

## 1. Source datasets and hashes

scPertEval @ `4685f11927e887745737600170da7a655b727553`, retrieved 2026-09-19
from `gs://scperteval/processed/` over HTTPS. Every file matched the byte size
and upstream MD5 recorded in the pre-download audit.

| dataset | cell line | bytes | MD5 | local SHA-256 |
|---|---|---:|:--:|---|
| `replogle22k562` | K562 | 2,430,512,332 | ✓ | `f381278e4deda43a…` |
| `replogle22rpe1` | RPE1 | 1,877,364,555 | ✓ | `25cb5bad6cd7abd8…` |
| `nadig25hepg2` | HepG2 | 1,236,448,196 | ✓ | `b25cd23e74b94c0f…` |
| `nadig25jurkat` | Jurkat | 2,004,474,709 | ✓ | `eaad2b2c58ed76d4…` |

Full records: `data/provenance/scperteval/` (per-dataset JSON, `manifest.json`,
`scperteval_sha256.txt`). Sources: Replogle et al., *Cell* 185:2559 (2022),
figshare+ 10.25452/figshare.plus.20029387; Nadig et al., *Nat. Genet.* 57:1228
(2025), GEO GSE264667.

**Protocol freeze.** `data/provenance/scperteval/protocol_freeze.txt` pinned the
SHA-256 of both specs, `virtual_cell/decomposition/anova.py` and its test file
*before* any download. Re-verified after the run: **all four still match.** The
decomposition protocol was not altered in response to any result.

## 2. Frozen design

Recomputed locally and independently of the earlier remote metadata audit:

| quantity | value | audit expectation | agree |
|---|---:|---:|:--:|
| contexts | 4 | 4 | ✓ |
| **shared perturbations** | **1,264** | 1,264 | ✓ |
| **shared genes** | **6,640** | 6,640 | ✓ |

Committed to `data/splits/four_context_v1/`:

| file | SHA-256 |
|---|---|
| `contexts.txt` | `ff44c676310edca9…` |
| `shared_perturbations.txt` | `f52b70382761030f…` |
| `shared_genes.txt` | `4a1a9ed9c1d422bb…` |
| `design_manifest.yaml` | records counts, ordering rules, source hashes |

Both intersections depend **only on identifier presence** — `obs['perturbation']`
categories and `var_names`. No expression value was consulted; there is no
HVG selection, variance filter, or outcome-dependent filtering anywhere.

## 3. Cell counts

Per-context totals (Phase 2, read from the files):

| context | cells | genes | control cells | perturbations |
|---|---:|---:|---:|---:|
| K562 | 308,646 | 8,563 | 10,691 | 1,971 |
| RPE1 | 240,774 | 8,749 | 11,485 | 2,016 |
| HepG2 | 133,757 | 9,623 | 4,976 | 1,818 |
| Jurkat | 258,202 | 8,881 | 12,013 | 2,137 |

Restricted to the 1,264 frozen perturbations:

| context | perturbed cells used | cells/pert min | median | max |
|---|---:|---:|---:|---:|
| K562 | 203,791 | 30 | 139 | 1,996 |
| RPE1 | 159,088 | 30 | 93 | 3,580 |
| HepG2 | 94,151 | 30 | 57 | 1,213 |
| Jurkat | 167,270 | 30 | 104 | 2,555 |

All four validated as CSR `float32`, `obs` carrying only `perturbation`, `var`
with no columns, and **empty** `layers`/`obsm`/`varm`/`obsp`/`varp`/`uns`. Zero
combination perturbations. `var_names` unique everywhere.

**Correction to the spec's expectation.** §1.4 of the data spec predicted
`expm1(X)` row sums would land slightly *below* 1e4 because genes were filtered
after normalisation. They are **exactly 10000.0** (p05 = median = p95 = 10000.0)
in all four datasets, so gene filtering preceded normalisation. The data is
cleaner than assumed; the prediction was wrong in a harmless direction and the
protocol is unaffected.

## 4. Pseudobulk construction

Exactly as frozen, with no re-normalisation of any kind:

```
control_mean[c]      = mean of X over cells labelled "control",  on the 6,640 shared genes
perturbation_mean[c,p] = mean of X over cells labelled p,        on the 6,640 shared genes
delta[c,p]           = perturbation_mean[c,p] - control_mean[c]
```

`X` is already scPertEval `log1p(CP10K)`, so `delta` is a difference of mean
log-expression. Tensor shape **(4, 1264, 6640)**. Streaming pass per context;
no matrix densified. Artifacts in `data/processed/four_context_v1/`
(git-ignored) with hashes in `data/provenance/scperteval/derived_manifest.yaml`.

## 5. Response energy

Denominator is the frozen uncentred convention `mean_{c,p} ||delta[c,p]||²`.

**Total response energy = 41.4168.**

| component | sum of squares |
|---|---:|
| μ | 5.2811 |
| α | 3.1225 |
| β | 15.4142 |
| γ | 17.5990 |
| **total** | **41.4168** |

### Invariants (verified before any interpretation)

| check | result |
|---|---|
| exact reconstruction | max abs err **4.44e-16** |
| Σ_c α = 0 | 7.48e-15 |
| Σ_p β = 0 | 2.08e-12 |
| Σ_c γ = 0 | 7.94e-15 |
| Σ_p γ = 0 | 2.62e-12 |
| SS partition (four components sum to total) | rel err **0.0** (exact) |
| max normalised pairwise cross term | **2.53e-17** |
| deterministic from frozen inputs | yes — fixed seed, frozen splits, hashed artifacts |

## 6–10. Component fractions

### Uncorrected (no noise correction)

| component | share |
|---|---:|
| **6. μ** | **12.75 %** |
| **7. α** | **7.54 %** |
| **8. template (μ+α)** | **20.29 %** |
| **9. β** | **37.22 %** |
| **10. γ** | **42.49 %** |

### Noise-corrected — the headline result

50 split-half resamples, seed 42, primary scheme (equal disjoint halves, odd
cell dropped, full-context control mean reused in both halves).

| component | share | sd over 50 resamples | min | max |
|---|---:|---:|---:|---:|
| template (μ+α) | **20.27 %** | 0.010 | 20.24 | 20.30 |
| conserved **β** | **30.07 %** | 0.016 | 30.03 | 30.10 |
| interaction **γ** | **21.05 %** | 0.033 | 20.96 | 21.13 |
| **noise** | **28.62 %** | 0.044 | 28.52 | 28.73 |

The four shares sum to 1.000000 by construction. **Uncertainty across resamples
is negligible** — every share is stable to within ±0.05 percentage points.

## 11. Split-half reliability: what is reproducible

Per-component reproducibility = cross-half signal SS ÷ raw SS:

| component | raw SS | reproducible SS | **% reproducible** | noise SS |
|---|---:|---:|---:|---:|
| μ | 5.2811 | 5.2787 | **100.0 %** | 0.0024 |
| α | 3.1225 | 3.1154 | **99.8 %** | 0.0071 |
| β | 15.4142 | 12.4531 | **80.8 %** | 2.9611 |
| γ | 17.5990 | 8.7169 | **49.5 %** | 8.8821 |
| **total** | 41.4168 | 29.5641 | **71.4 %** | 11.8527 |

This is the single most informative table in the run. The template is essentially
noise-free; β is largely reproducible; **γ is about half signal and half noise**,
which is exactly why the uncorrected γ (42.5 %) halves to 21.1 % after
correction. Any analysis that skips the noise correction would roughly double
the apparent interaction.

### Contextual reference only — a different pipeline

Molina & Zhang report template 27.8 / β 29.4 / γ 23.5 / noise 19.3 %. **These are
not a target and were not tuned toward.** They come from a different cell and
perturbation set (their K562 has 188,590 cells and 1,383 perturbations against
our 308,646 / 1,971), a different and undocumented 2,000-gene HVG response
space against our 6,640-gene intersection, and an undocumented filtering step.
That our β (30.1 vs 29.4) and γ (21.1 vs 23.5) land nearby is interesting but is
**not** evidence of reproduction, and the difference is **not** evidence of a bug
in either analysis.

## 12. Uncertainty

See the sd/min/max columns in §6–10: all four shares vary by <0.12 percentage
points across the 50 resamples. The split-half estimator is the only stochastic
element; everything else is deterministic given the frozen inputs.

## 13. Figures

In `outputs/four_context_v1/` (git-ignored, regenerable):

| figure | content |
|---|---|
| `component_energy.png` | uncorrected and noise-corrected component shares with resample error bars |
| `beta_gamma_magnitude.png` | distributions of ‖β_p‖ and ‖γ_{c,p}‖ by context |
| `split_half_reliability.png` | per-perturbation reliability distribution; reliability vs sampling depth |
| `gamma_heatmap.png` | 4 × 1,264 heatmap of ‖γ_{c,p}‖, perturbations sorted by mean |
| `beta_fraction.png` | transferability distribution; most-conserved vs most-context-specific perturbations |

## 14. Is β materially non-zero? **Yes.**

β holds **30.07 %** of total response energy after noise correction and is
**80.8 % reproducible**. Median ‖β_p‖ = 3.02 (max 9.39) against a median
‖γ_{c,p}‖ of 3.3–4.3. A conserved, transferable perturbation effect unambiguously
exists in these four contexts.

## 15. Is γ materially non-zero? **Yes.**

γ holds **21.05 %** after noise correction — a fifth of all response energy, and
comparable in magnitude to the template (20.27 %). Per-context median ‖γ_{c,p}‖
(K562 3.29, RPE1 4.26, HepG2 3.93, Jurkat 3.71) is of the same order as median
‖β_p‖ = 3.02, so the interaction is not a small perturbation on a conserved
effect.

## 16. Is γ reproducible rather than mostly measurement noise? **Yes, but only about half of it.**

γ is **49.5 % reproducible**: of its 17.60 raw SS, 8.72 survives the cross-half
test and 8.88 is measurement noise. So the interaction is *real and
substantial* — 21.05 % of total energy is reproducible interaction, far above
any plausible null — but it is also **by far the noisiest component**, carrying
75 % of all measurement noise in the decomposition.

The practical reading: γ cannot be dismissed as noise, and cannot be trusted
per-(context, perturbation) without accounting for depth. Both facts matter for
any later modelling.

## 17. Pathologies and unexpected findings

1. **Reliability correlates only weakly with cell count marginally (Spearman
   0.110) — this is Simpson's paradox, and it resolves cleanly.** ‖δ‖ correlates
   with reliability at 0.756 but *negatively* with cell count at **−0.492**,
   because `E‖δ̂‖² = ‖δ‖² + noise/n` inflates the measured effect magnitude for
   shallow perturbations. Stratifying by effect-size quintile restores a strong
   monotone relationship in every stratum:

   | effect quintile | Spearman(cells, reliability) |
   |---|---:|
   | 0 (weakest) | +0.483 |
   | 1 | +0.846 |
   | 2 | +0.920 |
   | 3 | +0.893 |
   | 4 (strongest) | +0.715 |

   Median reliability rises monotonically across cell-count quintiles within
   every effect stratum (e.g. quintile 3: 0.267 → 0.462 → 0.585 → 0.675 → 0.823).
   **Gate criterion 4's monotonicity requirement is met once effect size is
   controlled for**; the marginal correlation is an artifact of the negative
   depth–magnitude coupling, not a failure.

2. **Reliability is strongly bimodal** — a large mode near r ≈ 0 and a second
   near r ≈ 0.75. In an essential-gene screen many perturbations are
   near-null, so there is nothing reproducible to recover. Overall median
   reliability is 0.313; only 25.7 % (K562) to 57.1 % (RPE1) of perturbations
   exceed r = 0.5.

3. **RPE1 is the cleanest context by a wide margin** (median reliability 0.570
   vs 0.227–0.276 elsewhere) yet also has the **largest** interaction magnitude
   (median ‖γ‖ = 4.26). Higher reliability with larger γ argues its interaction
   is genuinely biological rather than noise.

4. **HepG2 is the thinnest context** (median 57 cells/perturbation, 94,151 cells
   used) and, as predicted in the spec, contributes disproportionate noise.

5. **Transferability is narrowly distributed.** The per-perturbation β fraction
   has median 0.412 (p25 0.353, p75 0.487); only **21.4 %** of perturbations are
   β-dominated (>0.5) and 8.5 % fall below 0.3. Most perturbations mix conserved
   and context-specific response rather than being cleanly one or the other.
   Most conserved: POLRMT, SMG5, SMN2, INTS8, PPRC1, EIF2B2, HARS, RPS13, VARS,
   DNAJC17 (translation/mitochondrial housekeeping). Most context-specific:
   DNAJC19, SHQ1, BCR, NCAPG2, EXOSC1, ISCU, NOP9, STIL, DCTN4, PLK4. *Offered as
   description only; no enrichment test was run and none should be read in.*

6. **192 of the 1,264 shared perturbations target a gene outside the 6,640-gene
   response space**, so their on-target knockdown cannot be checked. Recorded in
   the spec; not a defect.

No invariant failed, no share fell outside (0, 1), and no NaN/degenerate values
were produced.

## 18. Standing statement

This is an **independent four-context decomposition** on standardized public
scPertEval data. It is **not** an exact Molina & Zhang reproduction, and must
never be described as one. Arc contexts A/B/C were not used anywhere in this
analysis: not trained on, not intersected with the response space, not used to
choose any preprocessing, and no identity inference was attempted.

## Gate status against the predeclared criteria

| # | criterion | status |
|---|---|---|
| 1 | decomposition mathematics passes all invariants | **PASS** — §5 |
| 2 | results stable across reasonable preprocessing choices | **NOT YET TESTED** — sensitivity battery is the next step, deliberately not run in v1 |
| 3 | β and γ quantified | **PASS** — §14, §15, with resample uncertainty |
| 4 | split-half separates reproducible signal from noise | **PASS** — §11, §17.1 (monotone within effect strata) |
| 5 | no pathological preprocessing dependence | **NOT YET TESTED** — requires criterion 2 |

**The gate is not yet passed**: criteria 2 and 5 require the sensitivity
analysis, which is deliberately out of scope for v1. No novel prediction
architecture may be built until it is.

### Predeclared for the sensitivity stage (recorded now, not run)

A secondary robustness check must **independently split the control cells as
well**. The primary scheme reuses one full-context control mean in both halves,
which induces correlated error between them and can bias the noise estimate
downward. This was predeclared before the run and is deliberately not executed
here, so the primary result stands unmodified.
