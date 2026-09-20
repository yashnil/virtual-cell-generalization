# Pathway gamma — representation falsification battery v1

**Diagnostic only. No pathway predictor, D predictor, neural network, context
encoder, GNN, foundation-model predictor or Arc model was trained.** Independent
four-context study on public scPertEval data — **not** a Molina & Zhang
reproduction.

Date: 2026-09-20 · runtime 5.0 min · peak RSS 12.9 GB
Reproduce: `uv run python scripts/run_pathway_falsification.py`
then `uv run python scripts/plot_pathway_falsification.py`

**The question:** is the pathway-level gamma gain genuine biological pathway
organization, or just noise reduction from averaging genes?

All six freeze records maintained and re-verified (protocol 4/4, canonical 12/12,
decomposition 21/21, zero-shot 14/14, foundations 12/12, MSigDB 2/2). Nothing
earlier was modified.

---

## Design: every representation is the same arithmetic

Biological pathways, matched random gene sets and random projections are all
applied as a single linear map `W` of shape `(k, G)`, `R @ W.T`. Output
dimensionality and arithmetic are therefore **exactly matched**, and any
remaining difference is attributable to *which genes are grouped together*.

| null | preserves | destroys |
|---|---|---|
| **permuted** (primary) | number of sets, **every set size**, **every pairwise set-set overlap**, the multiset of per-gene membership counts — because `(MP)(MP)ᵀ = MMᵀ` | only which genes are grouped |
| **resampled** | number of sets, size distribution | overlap structure, per-gene degree distribution |
| **gaussian** | dimensionality only | everything else |

The permutation null is the strict one: it preserves the entire aggregation
geometry. Each invariance is asserted as a test in `tests/test_falsification.py`.
100 replicates per null type, seeds fixed in advance, **construction never
adjusted after seeing a result**.

---

## 1. Hallmark performance by fold

MSigDB 2024.1.Hs Hallmark, 45 of 50 sets with ≥10 genes in the frozen
6,640-gene space. Predeclared unweighted-mean aggregation. Baseline
`basal_affine`, as in the foundations study.

| held out | r_gamma | rho_full | ceiling | **normalised** |
|---|---:|---:|---:|---:|
| K562 | +0.570 | 0.764 | 0.874 | **0.652** |
| RPE1 | +0.182 | 0.878 | 0.937 | **0.194** |
| HepG2 | +0.026 | 0.642 | 0.801 | **0.032** |
| Jurkat | +0.510 | 0.727 | 0.852 | **0.599** |

Gene level, for reference: normalised 0.305 / 0.002 / 0.062 / 0.372.

## 2. Matched-random null distribution by fold

Reliability-normalised recovery, 100 replicates each:

| held out | Hallmark | permuted | resampled | gaussian | p (permuted) | z |
|---|---:|---|---|---|---:|---:|
| K562 | **0.652** | 0.305 ± 0.057 | 0.307 ± 0.058 | 0.300 ± 0.040 | **0.010** | **+6.1** |
| RPE1 | **0.194** | 0.015 ± 0.040 | 0.016 ± 0.042 | 0.015 ± 0.027 | **0.010** | **+4.5** |
| HepG2 | 0.032 | 0.065 ± 0.037 | 0.057 ± 0.040 | 0.065 ± 0.029 | 0.842 | **−0.9** |
| Jurkat | **0.599** | 0.378 ± 0.032 | 0.376 ± 0.034 | 0.363 ± 0.024 | **0.010** | **+6.8** |

`p = 0.010` is the floor attainable with 100 replicates under the
`(count+1)/(n+1)` convention — no null replicate ever reached the observed value
in those three folds.

## 3. Random-projection distribution by fold

Given in the table above (`gaussian`). It is statistically indistinguishable
from the two gene-set nulls, which is itself informative: **the aggregation
geometry — overlaps, gene degrees, sparsity — contributes nothing.** Only which
genes are grouped matters.

## 4. Reactome performance by fold

C2:CP:REACTOME 2024.1.Hs, 874 of 1,736 sets qualify, same aggregation rule,
against its own permuted null (20 replicates):

| held out | Reactome (normalised) | permuted null | p | z |
|---|---:|---|---:|---:|
| K562 | **0.536** | 0.305 ± 0.029 | **0.048** | +7.9 |
| RPE1 | **0.195** | 0.008 ± 0.024 | **0.048** | +7.9 |
| HepG2 | 0.091 | 0.063 ± 0.018 | 0.095 | +1.5 |
| Jurkat | **0.468** | 0.376 ± 0.018 | **0.048** | +5.1 |

`p = 0.048` is the floor for 20 replicates. **An independently constructed
ontology reproduces the pattern exactly**: significant for K562, RPE1 and
Jurkat, not for HepG2. Hallmark (45 coarse sets) beats Reactome (874 finer sets)
in the two strong folds, consistent with coarser aggregation buying more noise
suppression on top of the biology.

## 5. Reliability-normalised comparison — the decisive panel

| representation | K562 | RPE1 | HepG2 | Jurkat |
|---|---:|---:|---:|---:|
| genes (identity) | 0.305 | 0.002 | 0.062 | 0.372 |
| permuted null (mean) | 0.305 | 0.015 | 0.065 | 0.378 |
| gaussian null (mean) | 0.300 | 0.015 | 0.065 | 0.363 |
| Reactome | 0.536 | 0.195 | 0.091 | 0.468 |
| **Hallmark** | **0.652** | **0.194** | 0.032 | **0.599** |

**Random 45-dimensional aggregation reproduces gene-level performance almost
exactly** (K562 0.305 vs 0.305; Jurkat 0.378 vs 0.372). Raw `r_gamma` does rise
under random aggregation — from 0.185 to ~0.194 for K562 — purely because
averaging raises reliability, and the normalisation removes precisely that.

**So the entire benefit of dimensionality reduction per se is a reliability
artefact, and the reliability correction catches it.** What remains — Hallmark
and Reactome roughly doubling the normalised value — is attributable to
biological gene grouping.

## 6/7. K562–Jurkat dependence

Sensitivity diagnostic, not a benchmark. Hallmark, every two-source subset, so
that "lost a partner" is separable from "has fewer sources":

| target | all 3 | drop K562 | drop RPE1 | drop HepG2 | drop Jurkat |
|---|---:|---:|---:|---:|---:|
| K562 | +0.570 | — | +0.341 | +0.567 | **+0.337** |
| RPE1 | +0.182 | +0.171 | — | +0.011 | +0.182 |
| HepG2 | +0.026 | +0.003 | +0.116 | — | −0.078 |
| Jurkat | +0.510 | **−0.195** | +0.424 | +0.533 | — |

**6. K562 without Jurkat: +0.337** — down from 0.570 (−41 %), but still far
above the null mean of 0.194 raw. Losing RPE1 costs K562 about the same
(+0.341), so K562's signal is *not* specifically Jurkat-dependent.

**7. Jurkat without K562: −0.195** — a complete collapse, and the only negative
value in the table. Dropping HepG2 instead leaves Jurkat at +0.533, essentially
unchanged. **Jurkat's pathway gamma recoverability is entirely carried by
K562.** Reactome reproduces this (Jurkat 0.369 → −0.125; K562 0.437 → 0.297).

The two-source reduction is not the cause: RPE1 is unchanged at two sources
(0.182 → 0.182), and both K562 and Jurkat are unchanged when the *non*-partner
is dropped. The spread across dropped sources is 0.23 for K562 and **0.73** for
Jurkat.

Also visible: RPE1's modest signal depends on **HepG2** (drops to 0.011 without
it), a third and quite different pairing.

### A caveat I have been carrying since the zero-shot phase is now refuted

I repeatedly flagged that K562–Jurkat might be confounded by shared dataset
ancestry. **It is not.** K562 and RPE1 come from Replogle; HepG2 and Jurkat from
Nadig. So:

| pair | same dataset? | gamma excess over null |
|---|---|---:|
| K562–Jurkat | **no** (Replogle–Nadig) | **+0.196** |
| K562–HepG2 | no | +0.033 |
| RPE1–HepG2 | **no** (Replogle–Nadig) | −0.011 |
| HepG2–Jurkat | **yes** | −0.074 |
| K562–RPE1 | **yes** | −0.077 |
| RPE1–Jurkat | no | −0.081 |

The two **same-dataset** pairs rank 4th and 5th of six; mean excess is
**−0.076 same-dataset vs +0.034 cross-dataset**. Both informative pairings
(K562–Jurkat, RPE1–HepG2) cross the dataset boundary. Dataset ancestry does not
explain gamma sharing — if anything it is mildly anti-correlated. With only two
same-dataset pairs this is weak evidence, but it points firmly away from the
confound. The **lineage** confound for K562–Jurkat (both suspension leukaemia)
is untouched by this and remains.

## 8. Does pathway biology add information beyond dimensionality reduction?

**Yes, in three of four contexts, decisively.**

* Hallmark exceeds every one of 100 matched-random replicates in K562, RPE1 and
  Jurkat (z = +4.5 to +6.8 on the normalised metric, +5.2 to +16.5 on raw).
* The strict permutation null — which preserves set sizes, all pairwise
  overlaps and all gene degrees — is no easier to beat than the crude ones, so
  the gain is not aggregation geometry.
* Random aggregation is **identical to gene level** after reliability
  normalisation, so the gain is not dimensionality.
* An independent ontology reproduces the pattern.

**And no, not in HepG2**, where Hallmark (0.032) sits *below* the null mean
(0.065) at z = −0.9. HepG2 has no recoverable interaction at any resolution
tested, and pathway aggregation does not manufacture one. That is a reassuring
property of the analysis rather than a defect: the method is not producing
signal where none exists.

## 9. Is pathway-level gamma modelling scientifically justified?

Against the four predeclared criteria:

| # | criterion | verdict |
|---|---|---|
| 1 | biological representations beat matched random controls by a meaningful amount | **PASS** in 3/4 contexts (normalised roughly doubles; p = 0.010, the 100-replicate floor). **FAIL in HepG2.** |
| 2 | not entirely a reliability artefact | **PASS** — random aggregation reproduces gene-level performance exactly after normalisation; only real gene grouping moves it |
| 3 | some signal survives outside the K562–Jurkat pairing, or the limit is quantified | **PARTIAL PASS, quantified.** K562 retains +0.337 without Jurkat and RPE1 (+0.194, driven by HepG2) is independent of the pairing entirely. But **Jurkat is wholly K562-dependent (−0.195 without it)** and HepG2 has nothing. |
| 4 | an independent ontology behaves compatibly | **PASS** — Reactome reproduces the pattern fold-for-fold |

**Justified, with the scope stated explicitly.** Three of four criteria pass
cleanly and the third passes in quantified, partial form. The honest framing is
not "pathway gamma is predictable" but **"pathway-level gamma is recoverable for
some target contexts, and which ones is itself context-dependent and currently
unpredictable from four contexts."** A model must therefore report per-context
applicability rather than a single headline number, and must be evaluated
per-fold — a pooled average would hide that one of four folds is at zero and
another depends entirely on a single source.

## 10. Limitations from n = 4 contexts

* **Two of four folds carry the result.** K562 and Jurkat provide the large
  effects; RPE1 is modest; HepG2 is null. Any claim rests on effectively two
  strong observations.
* **Jurkat's result is one pairing.** Remove K562 and it is negative. It should
  never be quoted as independent evidence alongside K562.
* **The K562–Jurkat lineage confound is unresolved and unresolvable here.** Both
  are suspension leukaemia lines. Dataset ancestry is refuted (§6/7) but lineage
  is not, and with four contexts it cannot be.
* **Six context pairs.** The basal-similarity/gamma-sharing association and the
  dataset-ancestry refutation both rest on n = 6 with two same-dataset pairs.
* **Ceilings are themselves estimates** from 10 split-half repeats; pathway
  reliabilities (0.53–0.88) are high and stable, but the normalised numbers
  inherit their error.
* **Hallmark's 45 sets were not selected on performance**, but the choice of
  Hallmark over Reactome as the headline collection was made *after* seeing
  foundations results. Reactome is reported in full alongside for that reason.

## 11. Tests / lint / build

`pytest` **233 passed** (28 new, each null construction tested for exactly the
invariances it claims); `ruff check` passed; `ruff format --check` clean;
`uv build` succeeded.

---

## Recommendation

**Proceed to a pathway-level gamma model, scoped as follows**, or reconsider if
this scope is unacceptable:

* Evaluate **per held-out context**, never pooled.
* Report the **matched-random null** alongside every number, and the
  reliability-normalised value as primary.
* Treat **HepG2 as a known-negative control** — a model that "succeeds" there
  is suspect.
* Do **not** claim independence between the K562 and Jurkat results.
* Carry scale calibration (foundations §A) and source-agreement confidence
  (foundations §B) regardless.

## Standing statement

Independent four-context study on public scPertEval data. **Not** a Molina &
Zhang reproduction. Arc contexts A/B/C were not used anywhere: not trained on,
not intersected with the response space, not used to choose preprocessing, and
no identity inference was attempted. **No model was trained.**
