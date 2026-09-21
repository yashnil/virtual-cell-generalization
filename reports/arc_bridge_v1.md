# Arc Challenge Bridge — v1

**The output space is not the problem. The input space is.**

Mapping our validated findings onto Arc's 18,533-gene raw-count output turns out
to be straightforward: every gene the challenge's differential-expression
metrics actually test is measured in public data, and three count generators
reproduce the official controls' depth and sparsity. The blocker is upstream.
**Zero of the 300 Arc perturbation targets appear in any of the four contexts
the entire frozen research programme was built on.** Across all seven public
scPertEval datasets, 86 of 300 targets have any perturbation data at all, and
only 7 have it in more than one context — which is the minimum a conserved
perturbation effect can be estimated from.

A second finding compounds the first. The challenge's score is anchored so that
**0 is the context's mean perturbation response and 1 is a replicate of the
experiment.** Predicting no effect — which is what our conserved-transfer
predictor degenerates to when a perturbation has no source data — does not score
0. It scores **−0.31**, because a submission that calls nothing significant
scores 0 on direction fidelity against a baseline of 0.51.

So the honest position after this phase is: we can build a complete, valid,
submittable object today, and for 214 of 300 perturbations we have no evidence
that would let it beat the trivial baseline.

**No model was fitted. No submission was made. No conclusion from any prior
phase was modified.** Gamma modelling, pathway residual modelling and confidence
modelling all remain closed.

Date: 2026-09-20
Reproduce: `uv run python scripts/run_arc_bridge.py`
Requirements audit: `reports/arc2026_submission_requirements.md`
Outputs: `outputs/arc_bridge_v1/`

All 13 freeze manifests verified before and after this phase (140 files OK).
New code is confined to a new package, `virtual_cell.arc`; nothing frozen was
touched.

---

## 1. What exactly does an Arc submission have to be?

Audited in full in `reports/arc2026_submission_requirements.md` from the
installed `vcc 0.2.0` CLI and the official controls bundle, not from memory.
The short form: **360,000 cells × 18,533 genes of raw integer counts**, covering
all 300 targets in each of contexts A, B and C at exactly 400 cells each, with
no control cells, in the official gene order.

Validated end-to-end: `vcc sample --full` → 360,000 × 18,533 CSR float32,
integer-valued, then `vcc prep --dry-run` → *"targets: verified against the
official list"*, *"normalization: counts-preserved"*. The pipeline works. Only
the contents are missing.

The scoring metric was **not** available from the installed tooling, so it was
read from the published `cell-eval2` sources: `configs/vcc2026.yaml` and
`docs/vcc2026_metrics/vcc2026-metrics-brief.md` (2026/08/19, `rule_version` 3).
Six metrics, equally weighted per context:

| member | what it reads |
|---|---|
| `pds_cosine` | is a predicted effect closest to its own measured effect? |
| `expr_mse_unbiased_capped_norm` | squared profile error, sampling noise subtracted |
| `de_wilcoxon_direction_fidelity_yield_raw` | are the called genes' directions right? |
| `de_wilcoxon_direction_reach_raw` | how deep do correct directions persist? |
| `de_wilcoxon_sig_jaccard` | same responding-gene set? |
| `de_wilcoxon_lfc_nmae` | right fold-change magnitudes? |

Each is rescaled `s = (u − b) / (r − b)`.

---

## 2. How many of the 300 Arc perturbations are represented in public sources?

**86 of 300 (28.7%). In the four contexts this project's science was built on,
zero.**

| dataset | Arc genes measured | % of panel | Arc targets perturbed | % of 300 | a research context? |
|---|---|---|---|---|---|
| `arch1` | 18,017 | 97.2% | **13** | 4.3% | no |
| `kaden25rpe1` | 16,828 | 90.8% | **80** | 26.7% | no |
| `wessels23` | 15,439 | 83.3% | 0 | 0.0% | no |
| `replogle22k562` | 7,938 | 42.8% | **0** | 0.0% | **yes** |
| `replogle22rpe1` | 8,259 | 44.6% | **0** | 0.0% | **yes** |
| `nadig25hepg2` | 9,024 | 48.7% | **0** | 0.0% | **yes** |
| `nadig25jurkat` | 8,284 | 44.7% | **0** | 0.0% | **yes** |
| **union** | **18,367** | **99.1%** | **86** | **28.7%** | |

- Arc targets in **≥2** datasets: **7** (`HSBP1`, `MTA1`, `RNF2`, `SMARCA5`,
  `STAT6`, `TARBP2`, `ZNF714`) — all of them the `arch1` ∩ `kaden25rpe1`
  intersection.
- Arc targets in **0** datasets: **214**.

### Why the four research contexts contribute nothing

This is not a naming artefact; it was checked. Source-to-source perturbation
intersections are large and correct (K562 × RPE1 = 2,055 in the raw deposits),
NTC counts match their published values exactly, and all 300 Arc targets *are*
present as measured output genes in those files — they are simply never
perturbed. The reason is structural: Replogle and Nadig screen **essential-gene**
libraries, and Arc's 300 constructs are drawn from a largely non-essential
regulatory panel. The two libraries are close to disjoint by construction.

The consequence is blunt. Our decomposition writes a response as
`delta[c,p] = mu + alpha_c + beta_p + gamma[c,p]`, and our validated
contribution is that **`beta_p` (the conserved perturbation effect) transfers
across contexts while `gamma` does not**. `beta_p` is indexed by perturbation.
For 214 of Arc's 300 perturbations we cannot estimate `beta_p` from public data
at all, and for 207 more we have it in only one context, where `beta_p` and
`gamma[c,p]` are not separable.

**The axis of generalization we spent eleven phases characterising — unseen
*context* — is not the axis Arc's difficulty actually sits on for most of its
panel. For 71% of the targets, the problem is an unseen *perturbation*.**

---

## 3. What genes can we predict, and what happens to the rest?

`virtual_cell.arc.panel` assigns every panel gene an explicit support category
rather than zero-filling silently. For a target that *does* have public response
data:

| category | genes | meaning |
|---|---|---|
| `predicted` | 18,367 | measured in ≥1 public dataset |
| `measured_no_response` | 0 | |
| `unmeasured` | **166** | absent from every public gene space |

**The 166 unmeasured genes do not matter, and this is measurable rather than
assumed.** They are olfactory receptors, defensins, KIR cluster members,
testis- and keratin-specific genes — `GUCA1A`, `DEFB127`, `OR1E2`, `KIR2DL1`,
`PRAMEF8` and similar. Two facts settle them:

1. The four DE metrics test only genes above **5 CPM in Arc's own control
   cells**. Of the 166, **zero** clear that gate in any of A, B or C. They are
   never tested by four of the six members.
2. For the two members that read the whole axis (`pds_cosine`,
   `expr_mse_unbiased_capped_norm`), the 166 hold **0.0000%** of the control
   profile's squared norm at the metric's own normalisation.

Restricted to the 11,957 genes that clear the 5 CPM gate in at least one
context, public coverage is **11,957 / 11,957 = 100.00%**.

So: zero-filling the 166 is defensible, and it is recorded as a decision with
its justification rather than applied by default. `project_response` takes
`fill` as an explicit argument for exactly this reason.

---

## 4. What do the trivial submissions score?

Using the published anchors (section 8 of the metric reference, midpoints of the
A/B/C ranges) and the values the reference states analytically for a submission
that emits the control unchanged:

| member | baseline `b` | replicate `r` | control-emitting `u` | scaled score |
|---|---|---|---|---|
| `pds_cosine` | 0.500 | 0.956 | 0.500 | **0.000** |
| `expr_mse_unbiased_capped_norm` | 0.989 | 0.037 | 1.0032 | 0.000 (clamped) |
| `de_wilcoxon_direction_fidelity_yield_raw` | 0.514 | 0.814 | 0.000 | **−1.712** |
| `de_wilcoxon_direction_reach_raw` | 0.072 | 0.968 | 0.000 | −0.080 |
| `de_wilcoxon_sig_jaccard` | 0.029 | 0.399 | 0.000 | −0.078 |
| `de_wilcoxon_lfc_nmae` | 1.0013 | 0.400 | 1.000 | +0.002 |
| **mean of six** | | | | **−0.311** |

Against a mean-response baseline of **0.000** by construction.

Two things follow, and both are non-obvious:

**(a) "Predict no effect" is not the safe fallback.** It is a third of a point
worse than the panel mean, almost entirely through direction fidelity: a
submission that calls nothing significant scores a flat 0 on a member whose
baseline is 0.51. For the 214 targets with no evidence, the correct fallback is
the **context's mean perturbation response**, not the control.

**(b) The mean response is exactly the part of our decomposition we can
estimate.** Averaging `delta[c,p] = mu + alpha_c + beta_p + gamma[c,p]` over `p`
kills both `beta_p` and `gamma[c,p]` by the decomposition's own centring
constraints, leaving `mu + alpha_c` — the context main effect. Recovering a
*context* main effect from basal controls is precisely what the transferability
foundations phase studied. The scale's zero is therefore not an arbitrary
reference point for us: **it is a quantity our existing work is about.**

The achievable structure of an Arc submission is then:

```
prediction[c,p]  =  (mu + alpha_c)          <- scores 0 if estimated well;  all 300 targets
                 +  beta_p                  <- scores > 0;  available for 86, separable for 7
                 +  gamma[c,p]              <- not zero-shot predictable (frozen, phases 6-10)
```

A per-perturbation gain on 86 of 300 lifts the four per-perturbation means by
roughly `86/300 ≈ 29%` of the per-target gain.

---

## 5. Can we actually emit raw counts?

Yes. `virtual_cell.arc.generate` provides three generators, validated against
3,000 real control cells from each of A, B and C:

| generator | median library (control) | density (control) |
|---|---|---|
| `G0` resample controls | 20,602 (20,340) | 0.327 (0.323) |
| `G1` transport controls | 20,145 (20,340) | 0.320 (0.323) |
| `G2` negative-binomial count model | 19,952 (20,340) | 0.331 (0.323) |

All three emit non-negative integers, draw their per-cell depths from the
context's own controls, and sit far below the 10⁶ per-cell cap. Context B and C
numbers are in `outputs/arc_bridge_v1/generator_validation.csv`.

**One real defect was found and fixed.** `G1` transports a control cell by
multiplying its composition by the predicted fold change and re-drawing counts
at the same depth. Drawn from the cell's *own* empirical composition, that is a
second round of sampling loss on top of the one the measurement already made —
a gene the cell happened not to catch has probability zero of reappearing — and
it cost **17% of detected genes** (density 0.269 against the controls' 0.323).
That would have biased every Wilcoxon-based member of the score. Blending the
cell's composition with the pooled control composition (`smoothing=0.5`)
restores the detection rate to 0.320. The default is a modelling choice, stated
as one, and the regression test uses a fixture at realistic counts-per-gene
because the effect does not appear at unrealistic depth.

---

## 6. Can we measure a candidate before submitting it?

Yes. `virtual_cell.arc.metrics` reimplements all six scored members from the
published specification. It is not the official scorer and any disagreement is
this module's bug; what makes it trustworthy is that the reference states
several values analytically and each is asserted in the test-suite:

- a zero-effect prediction scores **exactly 0.5** on `pds_cosine`;
- a single shared profile scores **exactly 0.5** on the panel mean;
- a prediction carrying signal only in the excluded target-gene coordinates
  scores **exactly 0.5**;
- a zero-fold-change prediction scores **exactly 1** on `de_wilcoxon_lfc_nmae`;
- a submission calling nothing scores **exactly 0** on direction fidelity;
- emitting the control makes `expr_mse`'s numerator and denominator identical.

42 tests; full suite 356 passing.

---

## 7. What remains unsupported?

| | count | status |
|---|---|---|
| Arc genes with no public measurement | 166 | **resolved** — all below the 5 CPM gate, 0.0000% of profile energy |
| Arc genes in the DE-tested set with no public measurement | **0** | resolved |
| Arc targets with no public perturbation data | **214** | **open — the blocker** |
| Arc targets with data in exactly one context | 79 | `beta_p` not separable from `gamma[c,p]` |
| Arc targets with data in ≥2 contexts | 7 | usable |
| The 2026 metric definition | — | resolved, read from `cell-eval2` |
| Raw-count generation | — | resolved, three generators validated |
| Local scoring | — | resolved, six members reimplemented |

---

## 8. Should we acquire more raw data?

**Yes, and the decision is now well-posed rather than speculative.** Two
datasets carry every Arc target we have any evidence for:

| dataset | size | Arc targets | Arc genes |
|---|---|---|---|
| `kaden25rpe1` | 5.64 GB | 80 | 90.8% |
| `arch1` | 4.38 GB | 13 | 97.2% |
| total | **10.0 GB** | **86** (7 shared) | 99.1% union |

1.2 TB of disk is free. The four datasets already held contribute **nothing** to
Arc target coverage and should be kept only for the frozen science, not for the
bridge.

`arch1` is additionally the natural Arc-like public benchmark: 221,273 cells,
18,020 genes (97.2% of the Arc panel), 150 perturbations, and a context distinct
from all four research contexts — so it supports a genuine held-out-context
evaluation on almost exactly Arc's gene axis.

**This download has not been started.** It is a new acquisition beyond the four
approved in phase 4, so it waits for approval.

---

## 9. A warning that changes how our numbers should be read

Basal `log1p(CPM)` profile correlations between Arc's own contexts:

| | A | B | C |
|---|---|---|---|
| A | 1.000 | 0.780 | 0.749 |
| B | 0.780 | 1.000 | 0.842 |
| C | 0.749 | 0.842 | 1.000 |

The four public contexts sit at **0.89–0.93** (`zero_shot_recoverability_v1.md`).
A, B and C are **further apart from each other than any pair we have measured
transfer across**. Every zero-shot number in this project was therefore obtained
on easier context pairs than the ones Arc will score. Our results are an
optimistic bound on Arc performance, not a neutral estimate.

This is a statement about basal expression distance only. No identity inference
was attempted, and none of the above depends on knowing what A, B and C are.

---

## 10. What is the major gap before submission?

**Supervision for 214 perturbations, not machinery.**

Everything on the output side is done: the format is audited, the metric is
transcribed and tested, the generators reproduce the controls' depth and
sparsity, and the panel mapping accounts for every gene. A complete, valid
`.vcc` can be produced today.

What cannot be produced today is a *reason* for that file to beat −0.31 on most
of the panel. The frozen science gives us `beta_p` transfer; Arc gives us 300
perturbations of which 86 have any `beta_p` evidence and 7 have it twice. The
remaining 214 fall back to the context mean response, which scores 0 at best.

Three directions follow, in the order their premises are supported:

1. **Estimate `mu + alpha_c` for A, B and C from basal controls.** This sets the
   floor at 0 instead of −0.31 for all 300 targets, uses only permitted inputs,
   and is the one component the decomposition says is estimable. It is the
   cheapest real gain available and nothing about it reopens a closed question.
2. **Acquire `arch1` + `kaden25rpe1` (10 GB)** to lift `beta_p` coverage from 0
   to 86 targets, and to obtain an Arc-like held-out-context benchmark on a
   97%-overlapping gene axis.
3. **Recognise that the remaining 214 need a different kind of model** — one
   that predicts a perturbation's effect from the *identity* of the gene
   knocked down rather than from measurements of that same knockdown elsewhere.
   That is a perturbation-generalization model, not a context-generalization
   model, and this project has produced no evidence about it either way. It
   would be new work, and it is out of scope for this phase.

**STOP.** No submission has been made and none should be until direction 1 is
done and direction 2 is approved.
