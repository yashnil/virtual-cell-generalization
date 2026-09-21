# External unseen-perturbation validation — v1

**Kaden cannot answer the question, and finding that out is the result.**

`kaden25rpe1` was chosen to separate *study/target-regime shift* from *true
cellular-context shift*: it is an independent screen in RPE1, and Replogle RPE1
is already a source context. Under the frozen estimators it fails — unseen
`beta_p` scores `r = +0.016`, worse than predicting zero on energy. But the
**positive control fails too**: directly measured perturbations, the same 113
genes in the same cell line, transfer at `r = −0.033`.

The reason is not biology. **Kaden's own responses barely reproduce**:
split-half reliability 0.115–0.128, capping any method at `r ≈ 0.34`. `arch1`,
by contrast, reaches 0.84–0.95. When a benchmark's positive control is dead, its
negative result on the experimental arm carries no information.

**So this is Interpretation C, quantified: Kaden is uninformative, and the
study-shift-versus-context-shift question remains open.** It provides no
evidence for Interpretation B and none against the arch1 conclusion — which is,
if anything, strengthened, because arch1's reliability turns out to be
excellent.

For Stage B, the same reliability criterion decides the candidates, and it
overturns the obvious choice. **TeloHAEC has by far the better Arc overlap and
gene axis but weak per-perturbation signal — only 4% of its targets have ≥10
significantly changed genes. Feng's targeted screen has 100%.** The
recommendation is Feng, with its limitations stated plainly.

**No model was built or tuned. The frozen conclusion is unchanged:**
*STRING-based unseen-perturbation prediction succeeds internally but fails on
arch1.* Tier-0 policy is unchanged.

Date: 2026-09-21
Reproduce: `uv run python scripts/run_external_validation.py`
Outputs: `outputs/external_validation_v1/`

---

## 0. The protocol is provably the frozen one

`scripts/run_arch1_external.py` is frozen, so the shared protocol was
transcribed into `virtual_cell.modelling.external_benchmark` and the
transcription is **checked, not asserted**: the script re-scores `arch1` first
and refuses to touch Kaden unless every number matches the frozen
`arch1_external.csv`.

```
max |difference| against frozen arch1_external.csv          9.021e-17
max |difference| against frozen seen-perturbation reference 1.110e-16
REPRODUCED. The protocol is the frozen one; Kaden may be scored.
```

Estimators, `k = 25`, ridge `alpha = 10.0`, 64 PCA components, STRING v12.0,
the preprocessing and both confidence statistics are constants. Nothing was
fitted to Kaden.

## 1. Kaden unseen-perturbation result

| quantity | value |
|---|---|
| eligible unseen perturbations | **1,601** |
| shared gene axis | 6,640 (all of ours) |
| STRING coverage, training | 0.882 |
| STRING coverage, unseen set | 0.674 |

| estimator | beta `r` | beta unexplained | response `r` | response unexplained |
|---|---|---|---|---|
| U0 zero | — | 1.000 | +0.040 | 1.627 |
| U1 nearest | −0.000 | 4.061 | +0.007 | 3.701 |
| U2 k-NN | +0.016 | 1.583 | +0.012 | 1.056 |
| U3 ridge | +0.018 | **1.105** | +0.040 | 1.225 |

**No estimator beats predicting zero**, matching arch1 (+0.057) and, if
anything, weaker.

The feasible context main effect fails as it did on arch1: E2 leaves 14.34
unexplained, and **with a perfect oracle scalar 79.4% is still unexplained**
(`s* = 0.079`, norm ratio 4.14).

## 2. Kaden directly-measured transfer result

This is the control that should have worked.

| | `r` | unexplained |
|---|---|---|
| raw conserved transfer | −0.033 | 3.982 |
| scale-calibrated (s = 0.550) | **−0.033** | 1.893 |

On arch1 the same control gave `r = +0.302`, unexplained 0.966.

**Directly measured transfer fails on Kaden and works on arch1** — the opposite
of what a context-shift explanation predicts, since Kaden holds the cell line
fixed and arch1 does not.

Before interpreting that, the pipeline was checked:

- Kaden's basal control profile correlates with Replogle RPE1 at **r = 0.974**,
  the highest pair in the entire project. The cell line matches and the gene
  axis is aligned.
- For the *same 113 perturbations*, Kaden versus each source gives median
  `r` of −0.013 to −0.016, while Replogle RPE1 versus the other three cell
  lines gives **+0.157 to +0.205**.

So two independent screens of the same perturbations in the *same* cell line
agree *less* than screens of those perturbations in *different* cell lines.
That is not a plausible biological result, which is what prompted the
reliability check.

## 3. What Kaden says about study shift versus context shift

**Nothing, and the reason is measurable.**

Split-half reliability of each target context's own responses (Spearman-Brown;
the ceiling on any correlation is `sqrt(rho)`):

| dataset | subset | median cells/pert | Spearman-Brown | ceiling | frac > 0.2 | median ‖delta‖ |
|---|---|---|---|---|---|---|
| **arch1** | measured (17) | 790 | **0.952** | 0.976 | 1.000 | 3.02 |
| **arch1** | unseen (100) | 1,182 | **0.838** | 0.916 | 0.930 | 1.29 |
| `replogle22rpe1` | — | — | 0.570 | 0.755 | — | 6.07 |
| `nadig25jurkat` | — | — | 0.276 | 0.525 | — | 4.81 |
| `replogle22k562` | — | — | 0.243 | 0.493 | — | 3.60 |
| `nadig25hepg2` | — | — | 0.227 | 0.476 | — | 4.80 |
| **kaden25rpe1** | measured (113) | 379 | **0.128** | **0.358** | 0.204 | 1.97 |
| **kaden25rpe1** | unseen TF (300 sampled) | 400 | **0.115** | **0.339** | 0.257 | 1.92 |

Kaden's responses are the least reproducible in the study, by a factor of two
against the weakest source and by a factor of seven against arch1. Its effects
are also the smallest (‖delta‖ ≈ 1.9 against 3.6–6.1 in the sources) despite a
healthy 400 cells per perturbation — so this is weak knockdown or weak
transcriptional consequence, not thin sampling.

Disattenuating changes nothing: Kaden's measured transfer at `r = −0.033`
against a ceiling of 0.358 is still approximately zero, and arch1's `r = +0.057`
against a ceiling of 0.916 was never attenuation-limited in the first place.

**Verdict — Interpretation C.** Kaden is not evidence for B. It is an
uninformative benchmark whose target-side signal is too weak to discriminate
any hypothesis, and the study-versus-context question is still open. The one
thing it does establish is that **an external benchmark must be qualified on
target reliability before it is used at all** — a criterion this project did
not previously apply, and which now governs Stage B.

The confidence statistics degrade in the same way, as expected when the target
is noise:

| dataset | neighbour agreement | support distance |
|---|---|---|
| arch1 | **+0.314** | +0.195 |
| kaden25rpe1 | +0.027 | +0.025 |

## 4. TeloHAEC audit

**First correction: the accession in the plan is the pilot study.** GSE212396
is *"[Pilot scRNA-seq]"* — 50-gene and 200-gene libraries in TeloHAEC and
Eahy926, *"used to optimize conditions for and cross validate the larger
Perturb-seq study reported in GSE210681"*. The ~2,285-gene screen is
**GSE210681**, which is what was audited.

| property | value |
|---|---|
| accession | **GSE210681** (pilot: GSE212396) |
| context | CRISPRi TeloHAEC, telomerase-immortalized human aortic endothelial |
| design | 36,880 guides targeting gene promoters, CROP-seq, 20 10x lanes |
| perturbation targets | **2,345** columns in the released log2FC matrix |
| target selection | all expressed genes within 500 kb of CAD GWAS loci |
| measured genes | 25,054 (p-value matrix) / 17,472 (NMF spectra) |
| cells | 32,991 in one audited lane × 20 lanes ≈ **660,000** (~280/target) |
| control | non-targeting guides; log2FC released already referenced to control |
| representation | released as **log2 fold change**, not counts; raw counts in `RAW.tar` |
| publication | PMID 38326615 |

Assets: `ALL_log2fcs_dup4_s4n3.99x.txt.gz` 475 MB, `ALL_Pvalues...` 415 MB,
`RAW.tar` 6.4 GB, `aggregated...RDS.gz` 1.4 GB, MAST results 5.4 MB.

**Coverage is excellent:**

| | value |
|---|---|
| Arc 300 targets perturbed | **66** |
| Arc panel genes measured | 15,321 (82.7%) |
| shared with our 6,640-gene axis | **6,613** |
| in the 1,264-perturbation core (measured-transfer control) | 250 |
| unseen (absent from all four sources) | **1,925** |

**Signal is not.** Two independent statistics from the authors' own released
tables agree:

| statistic | result |
|---|---|
| program-level MAST, experiment-wide FDR < 0.05 | 304 / 2,357 targets (**12.9%**) with ≥1 significant program |
| gene-level, BH per perturbation, FDR < 0.05 | **93 / 2,345 targets (4.0%)** with ≥10 significant genes; **median 0** |

Restricting to targets with ≥10 significantly changed genes leaves **48**
unseen perturbations, **20** measured-transfer controls and **3** Arc targets.

This is consistent with the paper's own design: the authors aggregate into 60
NMF programs precisely because individual CAD-locus perturbations are weak in
endothelial cells. For a program-discovery study that is fine. For a
per-perturbation transfer benchmark it is the Kaden failure mode.

## 5. Feng 19-line audit

Published as Cell Genomics 2026 (PMC12903452; bioRxiv 2024.11.28.625833).

| property | value |
|---|---|
| context | human iPSC, CRISPRi |
| targeted screen | **444 genes**, 1,355 guides (3/gene), 20 non-targeting |
| cell lines | **20 lines from 10 donors** (genome-scale: 34 lines / 26 donors) |
| cells per target per line | median **74** (pooled over 20 lines ≈ 1,480) |
| measured genes | **6,520** (targeted screen) |
| raw data | ENA **ERP165335** |
| processed | Figshare `10.6084/m9.figshare.26819743` (9.12 GB, 11 files) |
| counts | Figshare `10.6084/m9.figshare.27989294` (7.93 GB, 9 files) |
| licence | MIT |

Key assets: `TargetedScreen_LFC_byGene.tsv.gz` **77.6 MB** (pooled over lines —
downloaded during this audit as an incidental file, 2,894,880 rows = 444 × 6,520),
`TargetedScreen_LFC_byGene-perLine.tsv.gz` **1.92 GB**,
`TargetedScreen_RNA-UMI-Counts.csv.gz` 4.48 GB, `TargetedScreen_Cell-Metadata.tsv.gz` 13.5 MB.

| | value |
|---|---|
| Arc 300 targets perturbed | **5** |
| Arc panel genes measured | 6,194 (33.4%) |
| shared with our 6,640-gene axis | 5,213 |
| in the 1,264-perturbation core (measured-transfer control) | **170** |
| unseen (absent from all four sources) | **141** |

**Signal is excellent**, and this is the decisive contrast:

| statistic | Feng targeted | TeloHAEC |
|---|---|---|
| targets with ≥10 significant genes (FDR < 0.05) | **444 / 444 (100%)** | 93 / 2,345 (4.0%) |
| median significant genes per target | **189** of 6,520 (2.9%) | 0 of 25,054 |

## 6. Arc-target overlap for each candidate

| candidate | Arc targets | with real signal | Arc panel genes |
|---|---|---|---|
| TeloHAEC GSE210681 | **66** | 3 | 15,321 (82.7%) |
| Feng targeted | 5 | **5** | 6,194 (33.4%) |
| *(kaden25rpe1, scored)* | 80 | unmeasurable | 16,828 (90.8%) |
| *(arch1, scored)* | 13 | 13 | 18,017 (97.2%) |

**Neither candidate adds meaningful Arc-target coverage in usable form.** That
is worth stating separately from the scientific question: these datasets are
for testing *whether unseen-perturbation prediction generalizes*, not for
improving Arc Tier assignments.

## 7. Eligible unseen test-perturbation counts

| candidate | nominal unseen | unseen with ≥10 significant genes |
|---|---|---|
| TeloHAEC GSE210681 | 1,925 | **48** |
| Feng targeted | 141 | **141** |
| *(kaden25rpe1)* | 1,601 | unmeasurable (reliability 0.115) |
| *(arch1)* | 100 | 93 (reliability > 0.2) |

## 8. Proposed second external benchmark

**Feng et al. targeted screen (444 genes, 20 iPSC lines).**

It is the only audited candidate whose **positive control can work**. Kaden
proved that a benchmark with a dead positive control cannot distinguish a
failed method from a failed measurement, and TeloHAEC would reproduce that
failure: 20 usable measured-transfer controls and 48 usable unseen
perturbations, most with a median of zero significantly changed genes.

Feng gives 170 measured-transfer controls and 141 unseen perturbations, all
with strong per-perturbation signal.

**Its limitations, stated plainly rather than discovered later:**

1. **Only 5 Arc targets.** It tests the scientific question, not Arc coverage.
2. **iPSC is pluripotent, as arch1 (H1 hESC) is.** It is a different lab,
   protocol and cell source, but it is not a maximally distant context, so a
   failure there is weaker evidence about *context* distance than a failure in,
   say, endothelium would have been.
3. **5,213 of our 6,640 genes**, against TeloHAEC's 6,613.
4. **Median 74 cells per target per line.** Per-line reliability will be poor;
   the pooled estimate (~1,480 cells) is the usable one, and per-line analysis
   must be treated as a separate, reliability-qualified question.

Its distinctive advantage: **20 cell lines** allow the degradation-across-
contexts test directly, with the study and protocol held fixed — which is
exactly the variable Kaden was meant to isolate and could not.

## 9. Decision rule used

Predeclared, applied before any model was run on either candidate, and with one
criterion added by Stage A's finding:

| # | criterion | TeloHAEC | Feng |
|---|---|---|---|
| 1 | single-gene CRISPRi | yes | yes |
| 2 | genuinely held-out cellular context | yes (endothelial) | partly (pluripotent, as arch1) |
| 3 | perturbation outcomes not previously used | yes | yes |
| 4 | sufficient globally-unseen perturbations | 1,925 nominal | 141 |
| 5 | ≥30–50 *useful* test perturbations | **48** | **141** |
| 6 | meaningful Arc-target overlap | 66 nominal / 3 usable | 5 |
| 7 | compatible transcriptomic output | 6,613 / 6,640 | 5,213 / 6,640 |
| 8 | accessible provenance | yes | yes |
| **9** | **target reliability / signal** — added from Stage A | **4% of targets** | **100% of targets** |

Criterion 9 is decisive and was **not** chosen to favour a result: it is the
direct lesson of Kaden, where a benchmark passing criteria 1–8 produced an
uninterpretable answer. TeloHAEC wins criteria 2, 6 and 7; it fails 5 and 9,
and 9 determines whether any of the others matter.

**No candidate was selected because the frozen model performs well on it.**
The frozen model has not been run on either. Both were assessed only on
published statistics and structural metadata.

## 10. Is a download warranted?

**Yes, for one 1.92 GB file, pending approval.**

`TargetedScreen_LFC_byGene-perLine.tsv.gz` (1.92 GB) is preferred over the
4.48 GB count matrix because it supplies everything needed:

- per-line deltas for 444 targets across 20 lines → a **cross-line reliability
  estimate**, which is a better control than split-half and is the criterion
  Stage A established;
- the pooled estimate by averaging lines;
- the degradation-across-contexts test.

Already obtained during this audit as incidental files (no further approval
needed): Feng `TargetedScreen_LFC_byGene.tsv.gz` (77.6 MB, pooled),
TeloHAEC MAST results (5.4 MB) and p-value matrix (415 MB), NMF spectra
(10 MB). TeloHAEC's `RAW.tar` (6.4 GB) and log2FC matrix (475 MB) were **not**
downloaded and are not recommended.

**Nothing further has been downloaded. This is the stop point.**

## 11. Does the Tier-0 zero-effect policy change?

**No.**

The policy changes only on new external evidence. Stage A produced none:
Kaden's reliability of 0.115 makes its result uninterpretable in both
directions — it neither supports nor refutes prior-based Tier-0 prediction.
The only valid external evidence remains arch1, which is negative, and whose
standing is *strengthened* by this phase because its reliability (0.84–0.95)
turns out to be the highest of any target context measured.

Arc policy is therefore unchanged:

| tier | n | prediction | confidence |
|---|---|---|---|
| **Tier 2** | 7 | direct conserved transfer, scale-calibrated | raw source agreement |
| **Tier 1** | 79 | direct single-context evidence; no prior-based extrapolation | raw source agreement (weak) |
| **Tier 0** | 214 | **zero perturbation-specific beta** | neighbour agreement, reported only |

No Arc submission has been created.

## Validation

`uv run pytest` — 411 passed. `ruff check`, `ruff format --check`, `uv build`
clean. All 19 freeze manifests verified.
