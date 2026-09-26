# C1: license-clean competitive atlas baseline (v1)

**Competition track, 2026-09-26. Nothing was submitted.** Submission #2 is still unspent.
Branch `competition/c1-license-clean`, starting HEAD `d375f935`. The start state is in
`c1_start_state.md`.

## Where things are

| artifact | path |
|---|---|
| license register | `reports/competition_v2/data_license_register.md` (+ `src/virtual_cell/competition_v2/licensing.py`) |
| X-Atlas permission status | `reports/competition_v2/xatlas_permission_status.md` (**PENDING**) |
| predeclaration, frozen before scoring | `reports/competition_v2/c1_predeclaration.md`, SHA-256 `fd4efa6f…eda32957` (`data/provenance/competition_v2/c1_predeclaration_digest.txt`), unchanged at the end |
| license-clean coverage | `data/splits/arc_target_support_license_clean_v1.csv`, `outputs/competition_v2/c1_license_clean/coverage_summary.json` |
| our backbone | `src/virtual_cell/competition_v2/{sources,atlas,fusion,generator,evaluation,licensing}.py` |
| equality with C0 | `tests/test_competition_v2_c1.py`, `outputs/…/c1_license_clean/sources/equality_*.json`, `…/equivalence_c0_bundle.json` |
| public folds | `outputs/competition_v2/c1_license_clean/folds/{H1,K562,CD4}/` |
| decision, tables | `…/c1_license_clean/{decision.json, comparison_public.csv, xatlas_ablation.csv, promoter_ablation.csv, agreement_transfer.csv}` |
| candidate (**not submitted**) | `…/c1_license_clean/{prediction.h5ad, prediction_compact.h5ad, c1_license_clean_val.vcc, c1_manifest.json}` |
| figures | `reports/competition_v2/figures/c1_{A..E}_*.{png,svg}` with `figures/sources/*.csv` |

Scripts, in run order (all under `scripts/competition_v2/`):

1. `verify_c1_state.py`
2. `prepare_c1_sources.py`
3. `check_c1_equivalence.py`
4. `build_license_clean_coverage.py`
5. `run_c1_public_folds.py`
6. `analyse_c1_folds.py`
7. `build_c1_candidate.py --stage emit|package|validate`
8. `c1_mechanistic.py`
9. `plot_c1_figures.py`

Attribution: the backbone concepts are AtlasShift's (kaipengm2/Virtual-Cell-Challenge-2026 @
`d24ce4f`, MIT). They are reimplemented in our own code, and every module docstring
names the upstream function it reproduces. The vendored copy was not modified. Only its
packaging utilities (`compact.py`, `pack.py`, which do lossless reordering and call
`vcc.prep`) were run, to package the C1 bundle.

---

## 1. Which data sources are GREEN?

| source | license | evidence |
|---|---|---|
| **VCC 2025 H1** (train / validation / test) | CC0 1.0 | Arc Virtual Cell Atlas page ("openly licensed via CC0 1.0"). The 2026 datasets page says "You are free to train on the H1 hESC data released for the 2025 Challenge" |
| **K562 GWPS** (Replogle 2022) | CC BY 4.0 | figshare API `license` field, re-fetched this phase |
| **CD4 GWCD4i DE** (Marson lab) | MIT License | CZI Virtual Cells Platform dataset listing. Caveat: the MIT text is written for "Software", and the authors' own bucket readme names no license |
| **GENCODE v47** | no named license; "no restrictions" | EMBL-EBI Terms of Use and the Ensembl disclaimer |
| Kaden 2025 RPE1 | CC BY 4.0 (Zenodo) | GREEN, but **excluded on scientific grounds**: not in the backbone; median split-half reliability 0.17; deposit names a CRISPRa screen |

## 2. Which are BLOCKED?

**X-Atlas/Orion HCT116 and HEK293T**: CC BY-NC-SA 4.0 (HF tag re-checked, plus the local
`LICENSE.md`). Status is BLOCKED_PENDING_PERMISSION. It is enforced in code:
`licensing.assert_sources_allowed` refuses `HCT116` / `HEK293T`, and a test asserts
that the candidate builder never references X-Atlas files.

## 3. Which are UNKNOWN?

**None.** Every source AtlasShift uses resolved to GREEN or BLOCKED.

## 4. How many of the 300 targets have direct evidence without X-Atlas?

**287 / 300 (95.7 %).**

| contexts with usable direct evidence | V1 | C0 (with X-Atlas) | **C1 (license-clean)** |
|---|---|---|---|
| 0 | 214 | 0 | **13** |
| 1 | 79 | 1 | **53** |
| 2 | 7 | 18 | **210** |
| 3+ | 0 | 281 | **24** |
| any | 86 | 300 | **287** |

* Usable targets per source: K562 269, CD4 251, H1 25.
* A target loses 2.21 contexts on average relative to C0.
* 296 targets have *some* GREEN measurement if Kaden is counted.

## 5. Which targets lose evidence?

These 13 targets end up with no GREEN direct evidence (all were V1 Tier 0, and X-Atlas
covered them all):

**ABCD1, ANKRD52, CAPRIN2, EPHB2, KIF21B, NICN1, PARP3, PBLD, PHF19, SEMA4F, SLC44A1, SNN,
TAF4.**

In the C1 bundle these targets are the context control plus the promoter cap. Their
pseudobulk effect norm is 0.13, against about 1.8 for supported targets. They pass the
cross-context reproducibility test only because their deviation equals minus the panel
mean, which is shared across contexts. That is **not** a target signature.

## Implementation correctness (§5–6)

* **Source statistics recomputed from raw by our code are bit-identical to upstream.**
  This covers K562 (437 × 8,246), H1 (300 × 18,080) and CD4 (3 × 437 × 10,282), every
  array, `max_rel = 0`. The CD4 `adjusted_p` array matches with NaN-equality.
* **Function-level equality** (`tests/test_competition_v2_c1.py`, 17 tests):
  * exact equality, dtype included, with the vendored code for source loading and
    shrinkage;
  * centred source response in both spaces, with and without centering;
  * fused effect, CD4 family and mapped fold change (`desired_mean`);
  * promoter pairs and cap, control template, and integer emission (same RNG);
  * the same checks repeated on a deterministic subset of the real C0 inputs, including
    X-Atlas, read only to prove the implementation.
* **End-to-end:** our backbone, fed C0's five sources and weights, regenerates
  **9/9 (context, target) blocks of the frozen C0 `prediction.h5ad` count-for-count**.
* **Harness cross-check:** C0 with X-Atlas on the H1 fold scores PDS **0.8270**, FID
  **0.5246**. These are the identical numbers from the previous phase's independent
  harness.

---

## Public folds (leave-one-atlas-out)

The held-out source never enters prediction, agreement, centering or selection.

| fold | truth | panel | predictors, C1 |
|---|---|---|---|
| H1 | H1 2025 train cells | 150 targets | K562 + CD4 |
| K562 | K562 GWPS cells | 390 targets | H1 + CD4 |
| CD4 | CD4 centred log2FC, effect-level PDS only | 346 targets | K562 + H1 |

Scores use local anchors (0 = the oracle mean response through the same emitter,
1 = split half). MSE is clamped.

### 6–7. How much of C0's advantage survives? (same folds, same harness)

| fold | model | coverage | sources/target | PDS raw | scaled PDS | MSE | NMAE | FID | REACH | JAC | **local Overall** | pred. DE (median) | dir. precision | DE yield |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| H1 | C0 (X-Atlas, **BLOCKED**) | 149/150 | 3.53 | 0.8270 | 0.655 | 0 | −0.070 | −0.143 | +0.009 | +0.035 | **+0.081** | 5,450 (real 1,296) | 0.541 | 1.0 |
| H1 | C0 w/o X-Atlas | 141/150 | 1.61 | 0.8277 | 0.656 | 0 | −0.146 | −0.125 | −0.010 | +0.061 | +0.073 | 6,051 | 0.541 | 1.0 |
| H1 | **C1a** | 141/150 | 1.61 | 0.8123 | 0.626 | 0 | −0.120 | −0.122 | −0.012 | +0.048 | **+0.070** | 5,712 | 0.545 | 1.0 |
| H1 | C1b (S2, f\* = 0.5) | 141/150 | 1.61 | 0.7777 | 0.556 | 0 | −0.089 | −0.117 | −0.018 | +0.030 | +0.060 | 5,336 | 0.549 | 1.0 |
| K562 | C0 (X-Atlas, **BLOCKED**) | 389/390 | 3.13 | 0.6915 | 0.507 | 0 | +0.128 | −0.149 | +0.231 | −0.000 | **+0.119** | 1,387 (real 3) | 0.548 | 1.0 |
| K562 | C0 w/o X-Atlas | 355/390 | 1.20 | 0.6659 | 0.440 | 0 | +0.124 | −0.158 | +0.140 | +0.002 | +0.091 | 1,377 | 0.547 | 1.0 |
| K562 | **C1a** | 355/390 | 1.20 | 0.6648 | 0.437 | 0 | +0.117 | −0.158 | +0.136 | +0.001 | **+0.089** | 1,378 | 0.548 | 1.0 |
| K562 | C1b (S2, f\* = 0.5) | 355/390 | 1.20 | 0.6468 | 0.389 | 0 | +0.089 | −0.150 | +0.113 | +0.001 | +0.074 | 1,279 | 0.549 | 1.0 |

The G0 control (zero effect) scores −0.062 (H1) and −0.052 (K562).

Survival, C1a relative to C0 on the same fold:

| measure | H1 | K562 | mean |
|---|---|---|---|
| scaled PDS (C0's PDS advantage) | 95.5 % | 86 % | ≈ 91 % |
| local Overall | 86 % | 74 % | **79 %** (+0.079 vs +0.100) |

Caveats:

* **The K562 fold's DE members are weakly informative.** K562 GWPS gives ~190 cells per
  target, and the truth calls a median of **3** DE genes. There FID/REACH/NMAE mostly
  measure over-calling against an almost empty reference. PDS is the reliable K562
  number.
* **Every arm over-calls DE.** This comes from the emission: 4-cell-averaged template
  cells, as in C0 (median 5.4–6.1k called against 1.3k real on H1). Directional
  precision is 0.54–0.55 everywhere, below the oracle mean-response arm. This is
  inherited from the C0 emitter and was not changed (§15).
* **Equal weights cost a little.** C1a (equal weights) is 0.002–0.003 Overall below
  C0-without-X-Atlas (upstream weights K562 2 / H1 2 / CD4 0.5) in both folds. It was not
  tuned, as predeclared.

**V1 is not on these folds.** Its only comparable numbers are:

* official hidden score **−0.062** (PDS +0.022, MSE 0, NMAE −0.006, FID −0.349, REACH
  −0.003, JAC −0.036);
* frozen public scPertEval LOCO (a different panel and folds): local Overall −0.090;
* 86 / 300 targets covered, median cross-context signature reproducibility 0.003.

It is not placed in the table above, because that would be a false apples-to-apples
comparison.

### 8. Does source-agreement shrinkage beat equal averaging?

**No. The pass rule did not fire, and C1a is selected.** Nested selection picked the
floor f\* = 0.5 in both outer folds, which is the least shrinkage in the grid.

| criterion (predeclared) | result | pass |
|---|---|---|
| 1. mean local Overall, C1b > C1a | 0.0671 vs 0.0793 | ✗ |
| 2. PDS better in ≥ 3/3 folds | 0/3 (H1 0.778 < 0.812; K562 0.647 < 0.665; CD4 tie 0.699) | ✗ |
| 3. no material DE worsening (≥ −0.02) | NMAE +0.002, FID +0.007, REACH −0.015, JAC −0.009 | ✓ |
| 4. gain in both folds | H1 −0.009, K562 −0.015 | ✗ |
| 5. no BLOCKED data | K562, H1, CD4, GENCODE | ✓ |

The global-shrinkage diagnostic (S1) is monotone the same way: less shrinkage is always
better, and λ = 0.75 is still below C1a.

**Correction to my own predeclared expectation.** I wrote that shrinkage would be "nearly
PDS-neutral in counts". It is not. At count level, PDS is measured against a
*reference* control that differs from the template pool. That pool-vs-reference offset
is shared by every prediction, so shrinking the effect lets the offset dominate. PDS
therefore falls monotonically with shrinkage (H1: 0.812 → 0.737 at f = 0). At the effect
level (CD4) the invariance held exactly, as predicted.

### 9. Does agreement generalise across held-out atlases?

**As a trust score, yes. As a shrinkage signal, no.** Raw source agreement predicts
held-out transfer quality in all three held-out atlases:

| fold | targets with agreement | Spearman ρ, agreement vs transfer cosine | ρ, agreement vs per-target PDS | median transfer cosine, low / mid / high tercile |
|---|---|---|---|---|
| H1 | 101 | **0.50** (p = 1e−7) | 0.24 | 0.026 / 0.038 / 0.068 |
| K562 | 112 | **0.49** (p = 3e−8) | 0.40 | 0.041 / 0.048 / 0.096 |
| CD4 | 112 | **0.34** (p = 3e−4) | 0.27 | 0.024 / 0.031 / 0.048 |

This replicates the frozen research finding ("confidence = raw source agreement") on new
atlases. Low-agreement targets still carry direction, though. Shrinking them costs PDS
and gains nothing on the DE members, so agreement belongs in a reported trust score, not
in the amplitude. On Arc, agreement is defined for 234 of 300 targets (median 0.017).

### 10. Is the promoter correction usable?

**Yes. GENCODE v47 is GREEN**, so the correction is retained exactly as frozen. It was
not retuned. The C1 Arc pairs are byte-identical to C0's `official_pairs.csv`:
75 pairs, 74 targets.

| fold | pairs on panel | PDS raw on / off | Δ scaled PDS | Δ local Overall | other members |
|---|---|---|---|---|---|
| H1 | 47 | 0.8123 / 0.8020 | +0.021 | +0.002 | NMAE −0.008; others ≈ 0 |
| K562 | 67 | 0.6648 / 0.6545 | +0.027 | +0.007 | REACH +0.010, NMAE +0.003 |

### 11. What is X-Atlas's marginal public value?

C0 backbone and weights, with vs without X-Atlas, on the same folds, untuned. The
ablation was not used for any C1 decision.

| fold | coverage | sources / target | PDS raw | local Overall |
|---|---|---|---|---|
| H1 | 141 → 149 | 1.61 → 3.53 | 0.8277 → 0.8270 (−0.001) | +0.073 → +0.081 (**+0.008**) |
| K562 | 355 → 389 | 1.20 → 3.13 | 0.6659 → 0.6915 (+0.026) | +0.091 → +0.119 (**+0.028**) |
| CD4 (effect) | 328 → 346 | 1.27 → 3.21 | 0.6989 → 0.6726 (**−0.026**) | — |

**+0.018 mean local Overall. PDS goes flat / up / down across the three folds, so it is
not uniformly positive.** On the Arc panel, permission is worth two further things:

* evidence for the 13 unsupported targets;
* genome-wide depth for the unseen D/E/F panel (18.3k targets per X-Atlas line, against
  K562 9.9k and CD4 11.5k).

### 12. Is C1 strong enough to justify submission #2?

**Yes, on public evidence.**

* C1a beats the local mean-response baseline in both folds: +0.070 and +0.089 local
  Overall, PDS 0.81 and 0.66. The zero-effect control scores −0.06 in both.
* It keeps ≈ 79 % of C0's public Overall and ≈ 91 % of its PDS advantage.
* It is the **only** candidate we may legally submit today, and it is the pipeline we
  would actually use for D/E/F.

The one real risk is the hidden panel. There, 13 targets are control-only and 53 rest on
a single context, so the public-to-hidden discount will be larger than on the folds.
Submitting it calibrates that discount on the pipeline that matters.

### 13. If X-Atlas permission arrives, should we add it back?

**Yes**, as a separate, predeclared change: C1 + X-Atlas with C0's frozen weights,
compared against C1 on these same folds. It is worth +0.018 public Overall and restores
the 13 unsupported targets and genome-wide D/E/F depth. It lowered CD4 effect-level PDS,
so check that fold explicitly. Do not add it silently, and record the permission first in
`xatlas_permission_status.md`.

### 14. What exact candidate should submission #2 use?

**`C1_LICENSE_CLEAN`, variant C1a** (equal-weight, centred fusion of K562 GWPS + VCC 2025 H1
+ CD4 DE, with the frozen amplitude 0.6 / 0.3, clip 3, promoter cap and C0-style
dual-moment emitter). No shrinkage, no X-Atlas, no Kaden.

| item | value |
|---|---|
| `prediction.h5ad` | `5072e45c23173d25fd51443a4debcb5f2d1c7d6dd100e8dcf537616b18f77b50` (360,000 × 18,533, nnz 2,437,123,261) |
| `prediction_compact.h5ad` | `cd560f5194a24b8243a4c95fc0e7b1b49fd474a24a55b0194cebe325a5985a12` |
| **`c1_license_clean_val.vcc`** | **`fcc4e2508798805d93b9f1bb3a6957ca7cc8fbaa9f8318296161e4184bc43ea7`** (3.36 GB) |
| git HEAD at build | `d375f935…` (C1 code uncommitted at build time) |
| sources / status | K562 GREEN, H1 GREEN, CD4 GREEN, GENCODE GREEN; excluded: HCT116 / HEK293T BLOCKED, Kaden GREEN (scientific) |
| checks | all 16 local checks PASS. `vcc prep --dry-run` exit 0: `verified_targets: true`, `normalization: counts-preserved`, nothing dropped, genes not reordered |
| submitted | **false** |

---

## 15. Generator (unchanged C0-style emission)

Emission fidelity is excellent: median L1 between the emitted and expected compositions
is 2.3e-4 (bulk) / 1.6e-3 (mean CPM) on H1. Cells are less realistic than real ones:

| | H1 fold: C1a / real perturbed | K562 fold: C1a / real | Arc bundle C1, A / B / C (C0 is identical to ±1) |
|---|---|---|---|
| library median | 54,985 / 53,631 | 7,790 / 7,712 | 20,037 / 20,118 / 20,199 (controls 20,109 / 19,946 / 20,034) |
| genes detected | 9,689 / 8,752 (+11 %) | 3,531 / 2,945 (+20 %) | 7,240 / 6,753 / 7,043 (controls 6,147 / 5,756 / 6,006, about +17 %) |
| sparsity | 0.468 / 0.517 | 0.550 / 0.623 | density 0.365 |
| within-group heterogeneity | 13.0 / 21.8 | 32.3 / 39.8 | 22.0 / 21.2 / 22.0 (V1 ≈ 33) |

Output validity holds: integers, no stored zeros, 400 cells per group, nnz 2.44e9
(under the 4.75e9 cap).

The 4-cell-averaged template (too many genes detected, too little heterogeneity) is the
likely driver of the DE over-calling in §6–7. This is the main biological-plausibility
debt, and a generator question for a later phase.

## Perturbation-signature diversity (Arc bundles)

| | V1 | C0 | C1 |
|---|---|---|---|
| deviation effective rank (A/B/C) | 288–292 (noise) | 235–237 | **201–208** |
| mean pairwise cosine of effects | 0.002 | 0.003–0.007 | 0.007–0.010 |
| median cross-context reproducibility | 0.003 | 0.91 | **0.78** |
| effect norm (median) | 3.9 (noise) | 1.41–1.46 | 1.79–1.91 |

C1 signatures are less diverse and less reproducible than C0's. There are fewer sources
per target, and 13 targets are control-only. Per-target effects are larger, because
there are fewer sources to average down (Fig D).

## 16. STATE

Not started (out of scope).

## 17. Figures

`reports/competition_v2/figures/`:

* **A** coverage;
* **B** public scores;
* **C** agreement vs transfer;
* **D** signature diversity;
* **E** X-Atlas ablation.

Each has PNG + SVG and a source CSV.

## 18. Permission

X-Atlas status: **PENDING**. No permission arrived during the run, and C1 was not
altered.

## 21. Validation

* `uv run pytest -q`: **587 passed, 1 skipped**.
* `ruff check`: clean. `ruff format --check`: clean.
* `uv build`: sdist + wheel OK.
* End-of-phase `verify_c1_state.py --full`
  (`outputs/competition_v2/c1_license_clean/state_full_end.json`): **0 failures**:
  * 16 freeze manifests + V1 digests (287 digests);
  * 7 raw checksum files (20 digests);
  * 8 competition downloads;
  * 8,689 X-Atlas files;
  * C0 code, statistics, prediction / compact / `.vcc` unchanged;
  * V1 bundle and `.vcc` unchanged;
  * V1 tier file unchanged.
* Predeclaration hash unchanged.

## Next step

**SUBMIT C1** (C1a, `c1_license_clean_val.vcc`, SHA-256 `fcc4e250…`). This is a
recommendation for the human, and the submission is theirs to make; nothing was
submitted here. Before submitting, confirm the CD4 MIT caveat is acceptable. After
scoring, compare the hidden result to the public expectation above, as was done for V1.
