# Paper figure design spec

Written 2026-10-07, before any paper-figure plotting code. Every number below is quoted from a frozen
report or read from a frozen output table. Nothing here is a new experiment. Where a panel needs a quantity
that is not stored verbatim, it is a **deterministic transformation of frozen per-draw or per-perturbation
values** (ratio, difference, median, percentile). Each such transformation is named in the "derived" column
and implemented in `scripts/paper_figures/build_fig*_sources.py`.

## 0. Paper-level thesis and novelty boundary

> Aggregate adaptation to a new cellular context can look mature before genuinely context-specific
> perturbation information has been learned. Response components differ sharply in how much target data they
> need. Extra source contexts mainly improve the zero-shot starting point rather than removing target-data
> requirements. Source compatibility can strongly influence transfer.

| claim | status in this paper | figure |
|---|---|---|
| δ = μ + α + β + γ decomposition; β more reproducible than γ | **background / replication** of Molina & Zhang 2026 (independent re-derivation, not a reproduction of their data) | Fig 1, ED Fig 1 |
| γ is hard to predict zero-shot | **not new** (Molina & Zhang; State) | Fig 1 caption only |
| target perturbations help γ | **not new** at ~30 % of the panel (Molina & Zhang; State); our contribution is the **budget curve** | Fig 2 |
| component-resolved target-budget learning curve across six contexts (template/scale vs γ⊥) | **primary contribution** (N1 confirmatory, N3-A preregistered replication) | Fig 2, ED Fig 2 |
| more sources shift the zero-shot start, not the γ learning rate | **exploratory** decomposition of a preregistered experiment (N3-B); the preregistered Outcome-B criterion fired but is judged confounded | Fig 3, ED Fig 5 |
| within-study source compatibility (GWPS ≫ RPE1) survives reliability/depth control | **confirmatory** (N5 Outcome A); magnitude largely expected from Replogle 2022 / Nadig 2025 | Fig 4 |
| cell identity vs study/lab | **unresolved** (N6 Outcome D) | Fig 4D, ED Fig 4 |
| source agreement as trust score | **demoted to a baseline** (N4 WEAK; v1 inflated) | ED Fig 3 |
| low-rank structure of learnable γ⊥ | **not novel** (Replogle 2022 programs, Systema, COMPASS) | not shown |

## 1. Shared definitions used in every caption

* **Context**: a cell line × screen. Six contexts: K562, RPE1 (Replogle 2022), HepG2, Jurkat (Nadig 2025),
  HCT116, HEK293T (X-Atlas/Orion, research use under CC BY-NC-SA 4.0).
* **k**: number of measured target perturbations ("anchors") available to the predictor; drawn at random from
  the anchor pool. Tested grid: k ∈ {0, 1, 2, 5, 10, 20, 50, 100, 200, k_ref}. k_ref = full pool: 885 (N1),
  743 (N3).
* **m**: number of source contexts (N3-B, m = 1…5; all 31 subsets enumerated).
* **γ⊥ recovery (M3)**: unbiased fraction of the *reliable* held-out γ⊥ energy explained, where γ⊥ is the
  target response component orthogonal to the source consensus for each perturbation. Template, global-scale
  and per-perturbation-scale predictors score exactly 0 by construction.
* **Template/scale recovery**: gain of estimator E2 (target template + scale from anchors) over the zero-shot
  E0s in M0 (full response, template included), expressed as a fraction of its k_ref gain. This is the
  preregistered k_T50 construct (N3 protocol §1.2).
* **E4**: the preregistered γ-capable estimator (E2 + kernel ridge of anchor residuals over source-defined
  perturbation similarity). All γ panels use E4.
* **Uncertainty**: (i) *draw interval* = 2.5–97.5 percentiles over 200 anchor draws × cell-split repeats
  (N1/N3-A) or 40 draws (N3-B); (ii) *perturbation bootstrap* = 2,000 resamples of the held-out test
  perturbations (N1/N3-A; N5/N6 paired bootstrap over 1,054 / 637 perturbations). No interval pooled across
  contexts exists in the frozen analysis; cross-context summaries are labelled **descriptive medians** and
  never carry a CI.

## 2. Visual system (implemented in `scripts/paper_figures/style.py`)

Font Arial (Helvetica fallback; STIX only for the ⊥ glyph); 7 pt body, 6 pt minimum for every text element (ticks, legends, annotations), 8.5 pt bold panel letters. Line 0.6 pt axes, 1.0 pt
data, 1.8 pt summaries. Widths: 183 mm (double column), 89 mm (single). White background, no top/right
spines, no gridlines except a faint reference line where a value is scientifically meaningful (0, a gate).

**Semantic palette — fixed across all figures** (categorical triple validated with the dataviz
`validate_palette.js`: lightness, chroma, CVD ΔE ≥ 9 all-pairs, normal-vision ΔE ≥ 18.8, contrast ≥ 3:1):

| meaning | colour | secondary channel |
|---|---|---|
| template / scale (global adaptation) | teal `#2A9D78` | square marker, dashed line where shown with γ |
| conserved β | blue `#1F6FB4` | circle marker, "β" label |
| context-specific γ / γ⊥ | vermilion `#D55E00` | diamond marker, "γ" label |
| measurement noise | light grey `#D0D0D0` + hatch | hatch, label |
| target-calibrated information / emphasised source | near-black `#1A1A1A` | filled marker, bold label |
| zero-shot / source-only / non-emphasised source | mid grey `#8F8F8F` | open marker, label |
| individual contexts in an aggregate panel | 30 % tint of the component colour | thin line; identity given in small multiples |
| failed gate / non-interpretable / off-scale | light grey `#BDBDBD` | italic label |
| positive control | dark neutral `#3A3A3A` | open square |
| uncertainty | 20–25 % alpha fill of the line colour, or 0.8 pt error bar | — |
| sequential magnitude (γ heatmap) | white → vermilion → dark brown (single hue, monotone L) | colourbar |

Six contexts are identified by **facet titles / direct labels**, never by six hues.

## 3. Main figures

### Figure 1 — Conserved and context-specific perturbation effects differ in reproducibility

*Claim:* a substantial context-specific component exists (21.0 % of response energy after noise correction)
but it is far less reproducible than the conserved component (49.5 % vs 80.8 %). **Background/replication.**

| panel | content | data | frozen upstream | type |
|---|---|---|---|---|
| A | schematic: sources (measured) → held-out target (controls only); zero-shot vs target calibration with k anchors | none | N1/N3 protocol §2.2 | schematic |
| B | δ(c,p) = μ + α(c) + β(p) + γ(c,p) with plain-language annotations; "following Molina & Zhang (2026)" | none | `four_context_decomposition_v1.md` | schematic |
| C | noise-corrected energy share: template (μ+α) 20.3 %, β 30.1 %, γ 21.0 %, noise 28.6 % (dots; sd over 50 split-half resamples ≤ 0.044 pp, drawn but smaller than the marker); uncorrected γ 42.5 % as an open marker | `paper/fig1_components.csv` ← `figure_sources/fig1_decomposition.csv` | `outputs/four_context_v1/summary.json` | confirmatory (decomposition v1, predeclared gates) |
| D | split-half reproducibility per component, same rows as C: μ 100.0, α 99.8, β 80.8, γ 49.5 %; plus the range over 21 preprocessing variants (β 77.4–85.0, γ 45.6–57.7) as a bar-free interval | `paper/fig1_components.csv` ← `fig1_decomposition.csv` + `outputs/four_context_sensitivity/variant_table.csv` | decomposition v1 + sensitivity battery | confirmatory |

*Constraints:* four contexts, two labs; the interval in D is a **range over preprocessing variants**, not a
CI. Do not show Molina & Zhang's percentages as a target. No pie chart.

### Figure 2 — Context-specific perturbation effects require substantially more target data to learn (flagship)

*Claim:* template/scale adaptation approaches its plateau within tens of anchors (k_T50 = 10–100), while γ⊥
recovery rises roughly log-linearly, reaches only 9–36 % of its k_ref value at k = 20 and is still rising at
k_ref = 743 (0.26–0.41 of reliable γ⊥ energy). **Primary contribution.**

Data: **N3-A** (six contexts, preregistered replication of N1 with an unchanged protocol; the four original
contexts re-run with 5 sources). N1 (four contexts, 3 sources, confirmatory FAIL-A) is shown in ED Fig 2D for
agreement.

| panel | content | data | type |
|---|---|---|---|
| A | template/scale recovery = [M0_E2(k) − M0_E0s] / [M0_E2(k_ref) − M0_E0s], median over draws; thin tint line per context, thick descriptive median over the six contexts; k ≤ 5 values below −1 (harmful: a single anchor's own response is mistaken for the template) drawn as off-scale markers with the true minimum in text | `paper/fig2_curves.csv` ← `figure_sources/n3/n3a_summary.csv` | confirmatory (preregistered descriptive k_T50 construct) |
| B | γ⊥ recovery M3(E4), median over draws, same x-axis and same layout as A; y from 0 (no normalisation — 1.0 would be all reliable γ⊥ energy) | same | confirmatory |
| C | 2 × 3 small multiples, one held-out context each; M3(E4) point estimate with **perturbation-bootstrap 95 % CI** at every tested k; identical axes; k_T50 and k_γ50 marked | `paper/fig2_curves.csv` ← `outputs/n3/n3a/n3a_boot.csv` | confirmatory |
| D | fraction of attainable gain G(k) = M3(k)/M3(k_ref) at k = 20, 50, 100, 200 per context (vermilion), next to the same ratio for template/scale (teal, open); point = ratio of medians (the preregistered C3 quantity at k = 20), interval = 2.5–97.5 % over draws of the per-draw ratio with k_ref paired by cell-split repeat; dashed line at the preregistered C3 threshold 0.25 (k = 20 only) | `paper/fig2_gain_fraction.csv` ← `outputs/n3/n3a/n3a_draws.csv` (**derived**) | ratio point confirmatory; its draw interval is a descriptive transformation |

*x-axis treatment (A–C):* log10 for k ≥ 1 with k = 0 placed one decade-half to the left behind an axis
break; tick at every tested budget (0, 1, 2, 5, 10, 20, 50, 100, 200, 743); faint verticals at k = 20, 50 and
k_ref. Markers at every tested budget; connecting lines only guide the eye.

*Constraints:* RPE1 is the exception (k = 20 reaches 36 %; PASS under the N1 rule). HEK293T has a negative
zero-shot start, so its template/scale gain is the largest. The template claim was **weakened** by N3 to
"10–100 anchors" — the figure must not imply ≤ 20. E4 is a deliberately simple linear smoother on half-depth
anchors: the curve is a property of this estimator class, not a lower bound for all models. No CI on the
six-context median.

### Figure 3 — Additional source contexts improve zero-shot transfer but do not eliminate target-data requirements

*Claim:* going from 2 to 5 sources raises the zero-shot template-removed full-response score by +0.10 to +0.15
in five of six targets (HEK293T falls by 0.16), while the 20-anchor gain in the same metric is constant or smaller
in five targets, and own-frame γ⊥ recovered from anchors changes by only −0.013 to +0.019 (k = 20) and −0.007 to
+0.027 (k = 50) from 2 to 5 sources (up to +0.045 from 1 to 5). **Exploratory** decomposition of the
preregistered N3-B ablation. The preregistered fixed-reference R_γ is *not* used as the main readout because
it conflates conserved-estimate error with γ learning (N3 §3.3); it is shown transparently in ED Fig 5.

| panel | content | data | type |
|---|---|---|---|
| A | schematic: held-out target, m of 5 sources (all 31 subsets), k anchors; two axes varied | none | schematic |
| B | 2 × 3 small multiples (one per target): template-removed full response M1 vs m, zero-shot (E0s, k = 0; bars = range over 5 repeats) **and** after k = 20 anchors (E4; draw interval) in the same frame-free metric, so the starting-point offset and the calibration gain (vertical gap) are read together; m = 1 zero-shot off-scale (single unshrunk source, fit scale fixed at 1.0: −0.19 to −5.17) | `paper/fig3_source_count.csv` ← `figure_sources/n3/n3b_R_gamma.csv`, `outputs/n3/n3b/n3b_draws.csv` (repeat range **derived**) | exploratory reading of preregistered quantities |
| C | 2 × 3 small multiples aligned with B: own-frame γ⊥ calibration gain R_γ,own(k, m) − R_γ,own(0, m) (= R_γ,own(k, m), since own-frame zero-shot ≡ 0) vs m at k = 20 (open) and k = 50 (filled), draw interval (2.5–97.5 % over 40 draws of the subset-averaged value) | `paper/fig3_source_count.csv` ← `outputs/n3/n3b/n3b_draws.csv` (**derived** intervals; medians asserted equal to the frozen `R_gamma_own`) | exploratory |
| D | heatmap small multiples (one per target, identical colour limits 0–0.42): rows m = 1…5, columns k = 5…743, colour = own-frame γ⊥ recovery; values < 0 (min −0.053, k = 5) drawn in the "under" colour | `paper/fig3_heatmap.csv` ← `figure_sources/n3/n3b_R_gamma.csv` | exploratory |

*Constraints:* do not say "more data does not help"; the large-budget own-frame plateau does rise from m = 1
to 3. The own-frame target itself depends on the subset's consensus direction (a reviewer-relevant caveat), which
is why B uses the frame-free M1. HEK293T is the exception in B. Source *reliability* and a basally similar *partner* (Q3, partial
ρ 0.22 and 0.62) — not diversity — are the useful properties; mentioned in caption only.

### Figure 4 — Within-study source compatibility strongly influences perturbation transfer

*Claim:* for the K562-essential target, an independent K562 screen from the same study (GWPS) has a
noise-corrected latent cosine of 0.822 [0.808, 0.837], against 0.38–0.56 for five other lines; the GWPS − RPE1
advantage survives reliability adjustment, depth matching, reliability matching and subsampling. **Confirmatory
(N5 Outcome A).** Cell identity is **not** separated from same-study compatibility.

| panel | content | data | type |
|---|---|---|---|
| A | forest plot: C_S with 95 % paired perturbation bootstrap for GWPS, depth-matched GWPS (open), Jurkat, HepG2, HCT116, HEK293T, RPE1, ordered by C_S; GWPS emphasised | `paper/fig4_compatibility.csv` ← `outputs/n5/n5_decision.json` | confirmatory primary / preregistered secondary |
| B | forest plot of the GWPS − RPE1 advantage under every preregistered control, in two blocks by estimand: **pooled latent cosine difference** (primary Δ_BC; depth-matched Δ_BC) and **per-perturbation template-removed Pearson difference** (raw median; Adj-2 fixed effects + reliability + magnitude, full depth; Adj-2 depth-matched; Adj-3 reliability-matched pairs, n = 108); zero line | `paper/fig4_robustness.csv` ← `outputs/n5/n5_decision.json`, `n5_summary.json` | confirmatory |
| C | ECDF of per-perturbation r_GWPS − r_RPE1 (n = 1,054): full depth (solid) and depth-matched (dashed); zero line; fraction favouring GWPS 0.858 / 0.705 and median annotated | `paper/fig4_per_perturbation.csv` ← `outputs/n5/per_perturbation_terms.csv` | preregistered descriptive |
| D | study-design matrix: rows GWPS, RPE1, VIPerturb-seq; columns same cell line / same study & lab / usable (reliability gate); VIPerturb row greyed "gate failed (N6, Outcome D)"; one-line statement of the unresolved confound | none (facts from N5/N6 reports) | schematic |

*Constraints:* never "same cell lines transfer better across labs". The ladder failure (RPE1 last despite
sharing study and library) is noted in the caption, not promoted. N6 C_S values are not shown anywhere in
Fig 4.

## 4. Extended Data figures

| fig | content | data | type |
|---|---|---|---|
| **ED 1** decomposition robustness | A: energy share per component across 21 preprocessing variants (strip, canonical highlighted); B: β vs γ reproducibility across the same variants; C: uncorrected vs noise-corrected γ share paired per variant; D: split-half reliability vs cells per half in the controlled depth experiment (median, bootstrap CI over 643 pairs) | `fig2_decomposition_robustness.csv`, `outputs/four_context_sensitivity/{variant_table,depth_experiment}.csv` | confirmatory (predeclared battery) |
| **ED 2** N1/N3 validity controls | A: anchor-permutation null (F3a), permuted vs unpermuted E4 M3 at k = 20, 100, six contexts; B: shared-control diagnostic (F3b), shared − split difference for M0 (template-level, inflated) and M3 (identical) at k = 20; C: a single anchor is harmful: E1 M0 gain over zero-shot at k = 1, 2, 5, 10, 20; D: N1 (4 contexts, k_ref 885) vs N3-A (k_ref 743) k = 20 / k_ref ratio with the 0.25 threshold | `outputs/n3/n3a/{n3a_f3a,n3a_f3b_summary,n3a_summary}.csv`, `n1_summary.csv` | confirmatory failure checks |
| **ED 3** agreement correction | A: paired before/after ρ(agreement, transfer quality): v1 averaged-halves 0.55–0.79 → per-repeat 0.29–0.59, with the frozen v1 report values; B: partial Spearman (T2) with 95 % CI, N1/N4 (4 folds) and N3 (6 folds); C: T1 difference ρ(agreement) − \|ρ(b)\| against magnitude, reliability, reliable energy and noise ceiling, N3 six folds | `n1_n4/{exploratory_halves_bias,n4_t1,n4_t2}.csv`, `outputs/n3/n4/{n4_t1,n4_t2}.csv` | A exploratory (labelled); B, C confirmatory |
| **ED 4** N6 reliability gate | A: pooled split-half reliability vs the frozen gate 0.10: VIPerturb 0.0866 (gate), GWPS at VIPerturb depth 0.109 and full-depth GWPS 0.341 (both exploratory, from `n6_results.md` §4); B: per-perturbation reliability ECDFs (VIPerturb, GWPS at matched depth, GWPS full); C: positive control C_GWPS full vs at VIPerturb depth with the d3 tolerance ±0.10; D: gate checklist d1–d4 | `outputs/n6/{source_gate.json,n6_decision.json,per_perturbation_terms.csv}` | gate = confirmatory; A's GWPS values exploratory |
| **ED 5** N3 reference-frame artefact | A: preregistered fixed-reference R_γ(k, m) for one target vs own-frame; B: zero-shot offset R_γ,ref(0, m) by target (negative for m < 5 by construction) | `figure_sources/n3/n3b_R_gamma.csv` | transparency for the B-criterion judgement |

ED 4 **omits every VIPerturb compatibility value**. They are non-interpretable after the gate failure and
showing them, even greyed, invites the reading "independent-lab K562 transfers worse", which N6 does not
support.

## 5. Interpretation constraints that every caption repeats where relevant

1. β/γ decomposition is Molina & Zhang's; ours is an independent re-derivation on public data.
2. No number is pooled across contexts for inference; cross-context medians are descriptive.
3. N1/N3 curves describe a simple linear estimator on half-depth anchors; better estimators could learn γ⊥
   faster. Random anchors only; selection was not tested (N2 not justified).
4. N3-B's preregistered Outcome-B criterion fired; the paper reports it and explains why the own-frame
   decomposition is the substantive reading.
5. N5 is within one study and lab; the K562 target and GWPS share protocols and cell stock.
6. N6 is Outcome D: the same-cell / different-study question is unanswered, not answered negatively.

## 6. Changes made during implementation (after this spec was first written)

* Fig 3 B/C became aligned 2 × 3 small multiples. The first draft of C (six targets overplotted with intervals) was
  unreadable; B gained the k = 20 curve so that the parallel-offset statement is made in one frame-free metric.
* Fig 2A annotation states k_T50 = 10–100 explicitly, so the panel cannot be read as "template solved by k = 20"
  (the N1-era claim that N3 weakened).
* Fig 4D became a 2 × 2 design grid (same cell line × same study) instead of a three-row table: the empty
  identification cell is the point of the panel.
* ED 1 omits a template reproducibility range (values 98.7–102.3 %; > 100 % arises from the split-half estimator
  with shared control-estimation error, explained in the sensitivity report) to avoid a misleading > 100 % mark.

## 7. Revision 2 (final structure)

* **Typography:** 6 pt floor everywhere (checked by `qc_figures.py`).
* **Schematics (1A, 1B, 3A, 4D)** are drawn on millimetre-scaled, equal-aspect canvases (`style.schematic_axes`) with
  one corner radius, one stroke weight and one arrow geometry (`style.box`, `style.arrow`, `style.cell_grid`). Fig 1A
  carries an explicit information-boundary key (available / anchor / held out).
* **Fig 2:** panel subtitles state the normalisation ("normalised to each context's own k_ref gain" vs "absolute
  scale, not normalised"); the caption opens with the normalisation sentence and the single estimator-scope
  statement; guides at k = 20, 50, 100, k_ref; D gains a descriptive cross-context median row.
* **Fig 3** now follows the requested sequence: A design; B zero-shot M1 vs m (preregistered quantity, per-target
  trajectories, descriptive median, HEK293T shown worsening); C own-frame calibration gain at k = 20 / 50 as paired
  m = 1 → 5 dumbbells with the exact change printed (−0.005 to +0.045; not described as flat); D heatmaps. Title
  changed to "Source breadth improves zero-shot transfer but does not replace target-context measurements".
* **ED 5** holds the full N3-B analysis: fixed vs own frame curves, the frame-free zero-shot vs 20-anchor view that
  was previously Fig 3B, and the preregistered Q1 Δ(20) beside the same pairing in the own frame
  (`ext5_delta20.csv`, derived from `n3b_draws.csv`; the fixed-frame values are copied from `n3_decision.json`).
* **Attributions** (Replogle 2022 CC BY 4.0; X-Atlas/Orion CC BY-NC-SA 4.0 research use; VIPerturb-seq CC BY 4.0)
  are in the captions that use those data.

## 8. Final freeze (2026-10-08)

* **Titles:**
  * Fig 1: Conserved and context-specific perturbation effects differ in reproducibility.
  * Fig 2: Context-specific perturbation effects require substantially more target data to *estimate*. The word
    "learnability" was removed because the result is estimator-specific.
  * Fig 3: Broader source panels improve zero-shot transfer across most target contexts.
  * Fig 4: Source choice strongly influences perturbation transfer. C_S is described as a noise-corrected
    similarity, not a biological compatibility, because cell identity is confounded with study/lab.
* **Fig 3 main-text claim** depends only on the frame-free zero-shot quantity R_full(0, m) (preregistered N3-B
  quantity, read descriptively):
  * A: design schematic.
  * B: zero-shot M1 vs m.
  * C: per-target m = 2 → 5 change, re-displaying B's values. No new estimand.
  * Supported claim: breadth improves zero-shot transfer in 5/6 contexts; HEK293T is the exception. No claim about
    substituting for target data.
* **ED 5 absorbs the former Fig 3C–D:**
  * A: preregistered fixed-frame Δ(20) vs own frame. This is the primary panel, and the disagreement is shown.
  * B: own-frame calibration gain vs m.
  * C: fixed vs own frame curves.
  * D: frame-free full response with 20 anchors.
  * E: own-frame k × m heatmaps.
* Validation gained a text-on-data check (`qc_figures.py`) and a full frozen-manifest / provenance / manifest
  re-hash (`validate_figures.py`).
