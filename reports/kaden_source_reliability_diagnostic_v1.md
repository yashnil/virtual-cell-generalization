# Kaden source-reliability diagnostic — v1

**The arch1 / Kaden disagreement is not measurement noise.** Kaden's individual
responses are weak, with a median split-half reliability of 0.15–0.17 against
arch1's 0.90. Averaged over its 1,836 perturbations, though, its main effect is
moderately reliable (0.76). So the two sources' main effects *could* agree up
to a noise ceiling of 0.86, yet they correlate at 0.095. Kaden also disagrees
with the same-cell-line Replogle RPE1 screen: per perturbation (median r ≈ 0
against a ceiling of 0.30), and in main effect by more than any two
*different* cell lines disagree among the four research screens. Under the
predeclared rules the evidence is **CASE E: mixed / inconclusive**. The frozen
Arc model is **not changed**.

Date: 2026-09-25 ·
Reproduce: `uv run python scripts/run_kaden_source_reliability.py` (27.6 min, peak RSS 14.6 GB; re-run once after whitespace-only formatting so the predeclared code hashes match the code that ran. Every result table was byte-identical) ·
Outputs: `outputs/kaden_source_reliability_v1/` ·
Figure: `reports/figures/fig7_source_reliability.png` ·
Tests: `tests/test_source_reliability.py`, `tests/test_kaden_diagnostic.py`

---

## A. What was and was not done

**Done:**
- A predeclared measurement-quality diagnostic of six public sources.
- Constants, gene axes, perturbation panels, quality bands and CASE rules were
  written to `predeclaration.json`, with the SHA-256 of the script and module,
  before any expression value was read. A test checks that the code that ran
  is the code that was declared.

**Not done:**
- No Arc model was evaluated, refitted or modified.
- No alternative `m_hat`, tier weight, magnitude calibration or G1 setting was
  tried.
- No Arc hidden outcome or leaderboard value was read. The script touches only
  the Arc gene list and the Arc target list, and a test pins that.
- The dry-run bundle was not regenerated.
- All 15 freeze manifests (258 digests) and all 19 raw-data checksums verified
  before the phase.

**Consistency check:** the diagnostic's full-data main effects reproduce the
frozen dry-run cosines to within 3e-17:

| cosine | recomputed | frozen |
|---|---|---|
| natural panels | 0.089149 | 0.089149 |
| Arc-target subsets | −0.002293 | −0.002293 |
| matched targets | 0.029519 | 0.029519 |

## B. Method

- **Split halves.** Each perturbation's cells were dealt at random into 20
  blocks, in one streaming pass per dataset. Each of **50 repeats** (seed
  20260925) takes a random 10-of-20 block half. Controls were blocked the same
  way and **split independently in every repeat**, following the protocol of
  `external_benchmark.target_reliability`. This keeps a 2.45e9-nonzero
  dataset out of memory; expectations match a cell-level split. The block
  construction is tested against direct cell means and against a planted
  latent signal.
- **Reliability.** Per-repeat Pearson between the two half responses, averaged
  over repeats, then projected to the full estimate with Spearman–Brown. The
  same was done for the perturbation-centred response, which is the quantity
  `beta_hat` uses.
- **Main effect.** The mean half-response over the *same* perturbation panel,
  giving `main_effect_half1` and `main_effect_half2`. Their Pearson, cosine and
  norm ratio were recorded per repeat.
- **Null.** 40 control pseudo-perturbations per dataset, built from control
  cells and referenced to the remaining controls. A perturbation's signal
  counts as *detected* when its reliability exceeds the null's 95th
  percentile.
- **Noise ceiling.** Two independent measurements of one latent response
  cannot correlate above `sqrt(R_x R_y)`. Disattenuation divides by that and is
  reported only where both reliabilities are ≥ 0.1. It is a model-based
  correction, not ground truth.
- **Gene axes** (identifier presence only):

  | axis | definition | genes |
  |---|---|---|
  | `G_ARC` | the frozen Arc response space | 16,494 |
  | `G_RPE` | Kaden ∩ Replogle RPE1 | 8,747 |
  | `G_3` | arch1 ∩ Kaden ∩ Replogle RPE1 | 8,166 |
  | `G_4CTX` | the frozen four-context axis | 6,640 |

- **Perturbation intersections**, recomputed locally:

  | intersection | perturbations |
  |---|---|
  | arch1 ∩ Kaden | **39** |
  | Kaden ∩ Replogle RPE1 | **205** |
  | Arc targets in both arch1 and Kaden | **7** (exactly the frozen Tier 2) |

## C. Per-perturbation measurement quality

Spearman–Brown reliability of the full-data response, all perturbations with
≥ 20 cells:

| dataset | axis | perts | median cells | q25 | **median** | q75 | detected above null | signal energy (median) |
|---|---|---|---|---|---|---|---|---|
| arch1 | `G_ARC` | 150 | 1,045 | 0.776 | **0.906** | 0.956 | 100% | 5.33 |
| arch1 | `G_3` | 150 | 1,045 | 0.773 | **0.900** | 0.958 | 99.3% | 4.07 |
| Replogle RPE1 | `G_3` | 2,016 | 82 | 0.371 | **0.677** | 0.833 | 91.6% | 29.3 |
| Kaden RPE1 | `G_ARC` | 1,836 | 400 | 0.096 | **0.167** | 0.274 | 81.3% | 0.80 |
| Kaden RPE1 | `G_3` | 1,836 | 400 | 0.089 | **0.155** | 0.258 | 84.2% | 0.77 |
| Kaden RPE1 | `G_4CTX` | 1,836 | 400 | 0.065 | **0.123** | 0.212 | 73.1% | 0.54 |
| Replogle K562 | `G_4CTX` | 1,971 | 125 | 0.122 | 0.396 | 0.682 | 85.7% | 5.62 |
| Nadig HepG2 | `G_4CTX` | 1,818 | 55 | 0.129 | 0.381 | 0.746 | 86.1% | 9.12 |
| Nadig Jurkat | `G_4CTX` | 2,137 | 89 | 0.170 | 0.423 | 0.650 | 87.7% | 10.17 |

The null's 95th percentile is 0.05–0.10 everywhere, and its median is ≈ 0.

- **The earlier single-split estimate replicates.** On the 6,640-gene axis
  Kaden's median is 0.123, against the frozen 0.115–0.128.
- **Kaden is weak, not empty.** 81% of its perturbations carry signal above
  the control null. But the median reliable signal energy is 0.80, against
  5.3 for arch1 and 29 for Replogle RPE1, despite 400 cells per perturbation
  (5× Replogle RPE1's 82). The deficit is effect size, not sampling depth.
- **Kaden is lower than every research context**, including the thinnest
  (HepG2, 55 cells per perturbation, median 0.38).

## D. Main-effect reliability: does averaging rescue it?

| dataset | axis | panel | perts | half-level Pearson (5–95%) | **Spearman–Brown** |
|---|---|---|---|---|---|
| arch1 | `G_ARC` | all | 150 | 0.958 (0.955–0.961) | **0.979** |
| arch1 | `G_ARC` | Arc targets | 13 | 0.916 | **0.956** |
| arch1 | `G_ARC` | matched Arc | 7 | 0.901 | **0.948** |
| Kaden | `G_ARC` | all | 1,836 | 0.617 (0.596–0.633) | **0.763** |
| Kaden | `G_ARC` | Arc targets | 80 | 0.253 (0.219–0.277) | **0.404** |
| Kaden | `G_ARC` | matched Arc | 7 | 0.176 | **0.300** |
| Kaden | `G_RPE` | shared with Replogle RPE1 | 205 | 0.390 | **0.561** |
| Replogle RPE1 | `G_RPE` | all | 2,016 | 0.989 | **0.995** |
| Replogle RPE1 | `G_RPE` | shared with Kaden | 205 | 0.978 | **0.989** |
| research contexts | `G_4CTX` | 1,264 shared | 1,264 | 0.897–0.992 | **0.946–0.996** |

Half-level norm ratios all sit near 1.0 (Kaden, natural panel: 5–95% range
0.88–1.16).

**Averaging partly rescues Kaden.** Its main-effect reliability rises from
0.17 per perturbation to 0.76 over 1,836 perturbations, which is the
averaging effect the phase anticipated. It is still the least reliable main
effect measured. Over the panels actually relevant to Arc it is poor: 0.40 on
its 80 Arc targets and 0.30 on the 7 it shares with arch1. The frozen `m_hat`
uses Kaden's full natural panel, i.e. the 0.763 value.

## E. Kaden vs Replogle RPE1, same cell line

Basal control profiles correlate at 0.974 (`external_validation_v1`), so the
cell identity and gene axis are right.

**Per perturbation**, 205 shared perturbations on `G_RPE`:

| response | median Pearson (95% CI) | median noise ceiling | median disattenuated (CI) | below half the ceiling | norm ratio (Kaden/Replogle) | sign agreement (top-100 genes) |
|---|---|---|---|---|---|---|
| raw | **−0.010** (−0.019, 0.005) | 0.298 | −0.065 (−0.104, 0.020), n = 133 | 90.7% | 0.33 | 0.50 |
| centred (β-like) | **−0.029** (−0.042, −0.016) | 0.304 | −0.102 (−0.128, −0.058), n = 145 | 94.6% | 0.33 | 0.47 |

**Main effect:**

| panel | Pearson | cosine | noise ceiling | fraction of ceiling |
|---|---|---|---|---|
| natural panels | 0.349 | 0.388 | 0.840 | 0.42 |
| 205 shared | 0.279 | 0.320 | 0.745 | 0.37 |
| *four research contexts, 6 pairs (different cell lines)* | *0.58–0.79* | *0.60–0.81* | *0.96–0.98* | *0.59–0.81* |

Kaden's responses are a third the size of Replogle RPE1's.

**Reading.** The per-perturbation ceiling is modest (0.30), so weak agreement
was expected. But the observed agreement is *zero*, with sign agreement at
chance (0.50). At the main-effect level Kaden agrees with its **own cell line**
at 42% of the ceiling, while four *different* cell lines agree with each other
at 59–81% of theirs. What signal Kaden has is not the signal a same-cell-line
screen measures. This is a **study- or protocol-level difference, not
noise**. The data cannot say whether the cause is knockdown efficiency,
library, preprocessing or biology.

## F. arch1 vs Kaden: the Arc-relevant comparison

On `G_ARC`:

| comparison | n (arch1 / Kaden) | cosine | Pearson | noise ceiling | Pearson / ceiling | verdict |
|---|---|---|---|---|---|---|
| main effect, natural panels (the frozen `m_hat` inputs) | 150 / 1,836 | 0.089 | 0.095 | **0.864** | 0.11 | not explained by noise |
| main effect, Arc-target subsets | 13 / 80 | −0.002 | −0.006 | 0.621 | −0.01 | not explained by noise |
| main effect, all 39 shared | 39 / 39 | 0.015 | 0.018 | 0.712 | 0.03 | not explained by noise |
| main effect, 7 matched Arc targets | 7 / 7 | 0.030 | 0.029 | 0.533 | 0.05 | not explained by noise |

A half-level check agrees. The cross-dataset half correlation on natural
panels is 0.085, against a half-level ceiling of 0.769.

Per perturbation, arch1 and Kaden agree at median r = −0.025 on all 39 shared
perturbations (CI −0.047 to −0.013) and on the 7 Arc targets. Noise ceilings
are 0.35–0.42, sign agreement is 0.49, and 97–100% of pairs fall below half the
ceiling. The raw arch1/Kaden norm ratio is 1.24–1.49, but raw norms are
depth-biased upward, so the signal-energy table in §C is the fair magnitude
comparison.

**The disagreement cannot be attributed to measurement noise.** Both main
effects are reproducible enough that, if they measured the same thing, they
would correlate near 0.86. They correlate at 0.1. So "Kaden is the problem"
is not what the uncertainty analysis shows: the two sources measure different
mean responses. Kaden also disagrees with its own cell line (§E), which makes
it the more anomalous of the two. But arch1 is a different cell type (hESC),
so its disagreement with Kaden is not by itself evidence against either.

## G. Arc-supported targets

Full table: `outputs/kaden_source_reliability_v1/arc_supported_target_reliability.csv`.
It has one row per (target, source): 93 rows covering all 86 supported
targets. **Nothing is removed and nothing is weighted.**

Quality bands were fixed before target identities were read. They apply to the
centred (β-like) Spearman–Brown reliability: **high ≥ 0.5, moderate 0.2–0.5,
low / uninformative < 0.2**.

| source | targets | high | moderate | low / uninformative | median reliability (centred) | detected above null (centred) | median cells |
|---|---|---|---|---|---|---|---|
| arch1 | 13 | **13** | 0 | 0 | 0.925 | 100% | 938 |
| Kaden | 80 | **5** | **26** | **49** | 0.158 | 85.0% | 430 |
| Kaden, Tier-1 only (73) | 73 | 5 | 24 | 44 | — | — | — |
| Kaden, Tier-2 (7) | 7 | 0 | 2 | 5 | — | — | — |

- **Measurable signal.** 85% of Kaden-supported targets (68 of 80) exceed the
  control null, but 61% are in the low band. Their signal is detectable and
  small.
- **Tier 2 averages unequal evidence.** Kaden's centred reliability on the 7
  Tier-2 targets is 0.02–0.48 (median 0.12); arch1's is 0.78–0.99. The frozen
  `beta_hat` averages the two with equal weight, so for Tier 2 roughly half of
  each prediction comes from the low-reliability source.
- Kaden's Arc targets are no worse than its other perturbations (median 0.153
  vs 0.168). Reliability rises only weakly with cell count (Spearman 0.19).

## H. Predeclared classification

| criterion (fixed in `predeclaration.json`) | value | verdict |
|---|---|---|
| **β source**: median centred reliability of Kaden Arc targets; fraction detected | 0.158; 0.85 | **partial** (median < 0.2, but detection ≥ 0.5, so not "unsuitable") |
| **`m_hat` source**: Kaden main-effect reliability, natural panel, `G_ARC` | 0.763 | **moderate** (between 0.5 and 0.8) |
| arch1–Kaden main-effect disagreement vs noise | ceiling 0.864, observed 0.095 | **not explained by noise** |
| **CASE** | — | **E: mixed / inconclusive** |

Two borderline facts, stated rather than argued away:
- Kaden's main-effect reliability (0.763) is close to the 0.8 "reliable" line.
- Its β evidence falls short of "unsuitable" only because most of its weak
  responses still clear a strict null.

Neither threshold was moved after the results.

## I. The twelve questions

1. **Kaden's perturbation-level reliability:** median Spearman–Brown **0.167**
   on the Arc response axis (0.155 on the shared three-source axis, 0.123 on
   the research axis). 59% of perturbations are below 0.2.
2. **Compared with arch1 and Replogle RPE1:** far lower. arch1 is **0.906**
   and Replogle RPE1 **0.677**, on the identical axis. Kaden is below all four
   research contexts too (0.38–0.42), despite 400 cells per perturbation.
3. **Is Kaden's main effect reliable?** Moderately: **0.763** over its full
   panel, the lowest of any source (the others are 0.95–1.00). On its Arc
   targets it is 0.40, and on the 7 matched targets 0.30.
4. **Does averaging rescue it?** Partly. Averaging 1,836 perturbations lifts
   reliability from 0.17 to 0.76, but not to the level of the other sources,
   and not on Arc-relevant subsets.
5. **Kaden vs Replogle RPE1 on matched perturbations:** no per-perturbation
   agreement (median r −0.010, ceiling 0.30; sign agreement 0.50). Main
   effects agree at 0.35, which is 42% of the ceiling.
6. **Disagreement left after accounting for noise:** most of it. The
   disattenuated per-perturbation agreement is about −0.07 (CI −0.10 to 0.02).
   Main-effect agreement reaches 0.37–0.42 of the ceiling, below the 0.59–0.81
   seen between *different* research cell lines.
7. **Kaden vs arch1 on shared Arc targets:** none. Per perturbation r −0.025;
   main effect 0.029 on the 7 matched targets and 0.095 on natural panels,
   against ceilings of 0.53–0.86.
8. **Kaden-supported Arc targets with measurable signal:** **85%** exceed the
   control null (68 of 80). Only **39%** (31 of 80) reach moderate or high
   reliability, and **6%** (5) are high.
9. **Is Kaden suitable as a β source?** Partial, under the predeclared rule.
   Its per-target responses are mostly detectable but low-reliability (median
   0.16), and they do not agree with either independent screen that shares
   perturbations with it.
10. **Is Kaden suitable as an `m_hat` source?** Moderate reliability (0.76).
    But reliability is not the obstacle: its main effect is reproducible
    enough to have agreed with arch1's, and it does not.
11. **Which CASE?** **E: mixed / inconclusive.**
12. **Is a reliability-qualified Arc v2 justified?** **Not by this
    diagnostic.** For `m_hat`, the disagreement is not a reliability problem,
    so excluding a source for unreliability would not address it. For β, the
    evidence is concerning (Tier-2 equal averaging with a source at 0.02–0.48
    reliability; 44 of 73 Kaden-only Tier-1 targets in the low band), but it
    does not reach the predeclared "unsuitable" bar. Any v2 must be motivated
    and predeclared on its own terms, and selected on public held-out
    contexts. It is **not implemented here**.

## J. Limitations

- Block split-halves give slightly smaller repeat-to-repeat variation than
  cell-level splits. Point estimates are unaffected, and the 5–95% bands are
  narrow anyway.
- Reliability is computed on the log1p(CP10K) representation that the frozen
  model uses. A count-space reliability could differ.
- The null uses pseudo-perturbations capped at half the controls. For arch1
  that means 477 cells rather than its median 1,045, which makes the null
  slightly wider and the detection threshold conservative.
- "Not explained by noise" rests on the independent-noise model behind
  `sqrt(R_x R_y)`. Shared artefacts within a dataset would *raise*
  within-dataset reliability without raising cross-dataset agreement. That is
  one possible form of the study-specific difference found here, not an
  alternative to it.
