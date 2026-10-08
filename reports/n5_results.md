# N5 results: within-study source compatibility for a K562 target

Date: 2026-10-07. Protocol: `reports/n5_protocol.md` (SHA-256 `43161f6e…`). Audit:
`reports/n5_coverage_audit.md`.

**Integrity, verified before writing:**

* Every file in `data/provenance/research_v3/n5_protocol_digest.txt` (protocol, N5 code, frozen axes,
  source-metadata manifest) re-hashes OK.
* Every file in `outputs/n5/manifest_sha256.txt` re-hashes OK.
* Nothing was re-run, tuned or modified for this report.
* The tables come from `outputs/n5/{n5_summary.json, n5_decision.json, matched_comparator.json,
  bootstrap_C.csv, per_perturbation_terms.csv, r4_calibration.csv}`.

**GWPS gates** (`outputs/n5/gwps/build_gates.json`):

* 343,127 cells kept of 343,127 (no cell below 200 genes), including 75,328 controls.
* G2 reconstruction error 9.1e-16; G3 median on-target delta −0.276 (n = 887); G4 asserted.
* Depth-matched variant: median 89 cells. 144 perturbations had fewer GWPS cells than RPE1, so all GWPS cells
  were used for them, as preregistered.

**Matched comparator** (source-only rule, written before the target was read): **Jurkat**. Its median
source reliability is 0.255, against GWPS's 0.283; gap 0.029.

**Attribution.**

* Replogle et al. 2022, *Cell* (figshare+ 10.25452/figshare.plus.20029387), CC BY 4.0; changes: pseudobulk
  statistics derived.
* X-Atlas/Orion (Xaira Therapeutics; Huang et al. 2025), CC BY-NC-SA 4.0, research-only
  (`reports/xatlas_license_memo.md`).

---

## 1. Primary result

The primary endpoint is the pooled noise-corrected (both-sided) latent cosine between template-removed source
and target response fields, over 1,054 perturbations × 6,408 genes. Intervals are 95 % paired perturbation
bootstrap (2,000).

| source | relation to K562-essential target | **C_S** | 95 % interval |
|---|---|---|---|
| **K562 GWPS** | same cell, same study, different screen | **0.822** | [0.808, 0.837] |
| **RPE1** | different cell, same study, same library | **0.379** | [0.358, 0.398] |

**Primary contrast.** Δ_BC = C_GWPS − C_RPE1 = **+0.444 [+0.420, +0.469]**.

**Preregistered verdict on the primary contrast: PASS.** The interval is far from 0.

## 2. Reliability-controlled result

**Adj-2: perturbation fixed effects, plus source split-half reliability (linear and quadratic) and log source
magnitude.** Outcome: per-perturbation template-removed Pearson.

| contrast vs RPE1 | full-depth GWPS | depth-matched GWPS |
|---|---|---|
| **K562 GWPS** | **+0.221 [+0.207, +0.235]** | **+0.181 [+0.168, +0.194]** |
| HepG2 | +0.069 [+0.061, +0.077] | +0.066 [+0.057, +0.074] |
| Jurkat | +0.098 [+0.089, +0.107] | +0.093 [+0.083, +0.102] |
| HCT116 | +0.063 [+0.051, +0.075] | +0.058 [+0.045, +0.071] |
| HEK293T | +0.039 [+0.024, +0.053] | +0.033 [+0.018, +0.048] |

Intervals are 95 % perturbation-cluster bootstrap.

**The same-cell advantage survives adjustment.** It shrinks by about a fifth when GWPS is depth-matched
(0.221 → 0.181) but stays large and precisely estimated.

**Reliability is not what separates the sources:**

| source | median source reliability | C_S |
|---|---|---|
| RPE1 | 0.531 (the highest) | 0.379 (the lowest) |
| full GWPS | 0.283 | 0.822 |
| depth-matched GWPS | **0.146** | **0.818** |

The depth-matched GWPS is *less* reliable than RPE1 and leaves C essentially unchanged. The noise-corrected
endpoint is behaving as designed: a 2.5-fold depth reduction moved it by only 0.005.

## 3. Matched analyses

| analysis | result |
|---|---|
| **Adj-3: reliability-matched perturbations** (\|ρ_GWPS − ρ_RPE1\| ≤ 0.05; n = 108) | median r_GWPS − r_RPE1 = **+0.233 [+0.200, +0.273]** |
| **Adj-4: GWPS subsampled to RPE1's cell count per perturbation** | C = 0.818 [0.802, 0.833]; Δ_BC = **+0.439 [+0.414, +0.464]** |
| raw (unadjusted) median per-perturbation r, GWPS − RPE1 | +0.138 [+0.125, +0.149] |
| raw, depth-matched GWPS − RPE1 | median r 0.246 vs 0.181 |
| predeclared matched comparator: C_GWPS − C_Jurkat | **+0.264 [+0.244, +0.287]** |
| C_GWPS − mean of the 5 non-K562 sources | +0.368 [+0.351, +0.387] |

**No conclusion changes under any matched analysis.**

## 4. Per-perturbation distribution (descriptive; no subgroups)

Per-perturbation difference of template-removed Pearson, r_GWPS − r_RPE1, over n = 1,054:

| | fraction favouring GWPS | fraction favouring RPE1 | 5 % | 10 % | 25 % | median | 75 % | 90 % | 95 % |
|---|---|---|---|---|---|---|---|---|---|
| full GWPS | **0.858** | 0.142 | −0.064 | −0.023 | +0.042 | +0.138 | +0.260 | +0.388 | +0.455 |
| depth-matched GWPS | 0.705 | 0.295 | −0.111 | −0.076 | −0.011 | +0.066 | +0.171 | +0.282 | +0.391 |

* **The advantage is broad, not driven by a few perturbations.** At full depth GWPS is better for 86 % of
  perturbations; depth-matched, for 71 %.
* The pooled endpoint is energy-weighted. The top 10 % of perturbations by contribution account for 58 % of
  the aggregate difference, close to their 48 % share of target reliable energy.
* Removing the strongest contributors leaves the contrast large:
  * top 1 % removed: Δ_BC 0.433;
  * top 5 %: 0.411;
  * top 10 %: 0.392;
  * top 25 %: 0.336.
* Against the other sources, the fraction of perturbations favouring GWPS is: HepG2 0.90, Jurkat 0.85,
  HCT116 0.94, HEK293T 0.96.

## 5. Other sources (preregistered secondary / descriptive rules only)

**C_S (95 % interval):**

| source | C_S |
|---|---|
| Jurkat | 0.558 [0.540, 0.575] |
| HepG2 | 0.490 [0.472, 0.507] |
| HCT116 | 0.456 [0.433, 0.478] |
| HEK293T | 0.387 [0.369, 0.406] |
| RPE1 | 0.379 [0.358, 0.398] |

**Ladder (N5-D, preregistered ordered differences):**

| step | median difference [95 %] | holds |
|---|---|---|
| GWPS > RPE1 | +0.444 [+0.420, +0.469] | ✓ |
| RPE1 > mean(HepG2, Jurkat) | **−0.145 [−0.163, −0.129]** | **✗ (significantly reversed)** |
| mean(HepG2, Jurkat) > mean(HCT116, HEK293T) | +0.102 [+0.084, +0.119] | ✓ |

**The ladder does not hold.** The same-study, same-library different cell line (RPE1) is the *least*
compatible non-K562 source. It is below both companion-study lines from the same lab, and tied with the
other-lab HEK293T.

**Secondary endpoints:**

| source | R2 raw cosine | R3 directional accuracy |
|---|---|---|
| GWPS | 0.463 | 0.662 |
| Jurkat | 0.306 | 0.600 |
| HepG2 | 0.271 | 0.587 |
| RPE1 | 0.241 | 0.586 |
| HCT116 | 0.222 | 0.567 |
| HEK293T | 0.182 | 0.566 |

**R4 (single-source calibration, E4; medians over 40 draws):**

| source | M1 (full response, template removed) at k = 20 | M1 at k = 50 | own-frame γ⊥ at k = 20 | at k = 50 |
|---|---|---|---|---|
| GWPS | 0.371 | 0.384 | 0.016 | 0.041 |
| GWPS (depth-matched) | 0.250 | 0.270 | 0.020 | 0.058 |
| Jurkat | 0.163 | 0.183 | 0.025 | 0.044 |
| HepG2 | 0.128 | 0.175 | 0.025 | 0.085 |
| RPE1 | 0.105 | 0.155 | 0.008 | 0.050 |
| HCT116 | 0.094 | 0.143 | 0.019 | 0.068 |
| HEK293T | 0.054 | 0.094 | 0.005 | 0.042 |

* A compatible source gives a far better *starting point*.
* With every source, target-specific (own-frame) γ⊥ learned from 20–50 anchors stays small (0.005–0.085).
  This agrees with N1/N3.

These are descriptive and are not promoted to the headline.

## 6. Outcome classification

**Mechanical preregistered verdict: Outcome A ("same-cell advantage survives reliability control").** All
three A conditions hold:

* Δ_BC interval > 0;
* the Adj-2 GWPS − RPE1 interval > 0;
* the depth-matched Δ_BC interval > 0.

The other outcomes:

* B is false (A holds).
* C′ is false (Δ_BC > 0).
* D′ is false: the C_GWPS upper bound is 0.837, not < 0.70.
* E is false (spread 0.44).
* Kill rule: not triggered.

**Scientific interpretation.** Within one study and lab, an independent K562 screen predicts held-out K562
responses with a noise-corrected latent cosine of 0.82. Other cell lines reach 0.38–0.56. This advantage is
not explained by:

* source reliability or depth: it survives fixed-effect adjustment, matched pairs and depth matching, and
  persists when the K562 source is the *less* reliable one;
* perturbation signal strength: perturbation fixed effects absorb it;
* a handful of perturbations: 71–86 % favour GWPS;
* library or study sharing: the source that shares the target's study **and** library (RPE1) is the least
  compatible non-K562 source.

**Where the verdict and the interpretation agree, and where they don't.** They agree that the effect is real
and robust. They differ in what it means:

* "Outcome A" is worded as biological compatibility beyond reliability.
* The design cannot separate cell identity from *same-cell-line, same-lab* compatibility, because the only
  K562 source is from the target's own study.
* The defensible statement is therefore **within-study source compatibility**: a same-line screen from the
  same study transfers far better than any other line, beyond reliability.

The ladder failure adds a second, descriptive observation. Among non-K562 sources, sharing the target's
study and library confers no advantage: RPE1 is last. This argues against a pure study/library-compatibility
explanation *for the non-K562 ordering*. It does not identify cell identity as the cause of the GWPS
advantage.

**The headline magnitude has a ceiling reading.** Two K562 screens from one lab share only about two-thirds
of their latent response variance (0.82² ≈ 0.68). About a third of the latent response field differs even
within cell line, study and lab, once timepoint and library scale change. That is not the D′ regime (< 0.70
cosine), but it is a substantial screen-level floor on "context-specific" error.

## 7. Confounding

1. **Same study/lab.**
   * GWPS and the target are both Replogle 2022 / Weissman lab: shared protocols, cell stock, reagents,
     processing and possibly shared technical response signatures.
   * Any K562-specific technical artefact common to both screens would inflate C_GWPS.
   * This is the largest unresolved confound.
2. **GWPS differs in library scale and timepoint** (~day 8, genome-wide library vs ~day 6, essential
   library).
   * This should *penalise* GWPS, so it makes the advantage conservative.
   * But a longer timepoint can also amplify K562-specific secondary programs. The direction of that bias on
     the cosine is not identifiable.
3. **RPE1 shares the target's library and study.** If shared library or protocol drove transfer, RPE1 should
   lead the non-K562 sources. It is last. That argues against library/protocol sharing as the dominant
   driver, but it is one line pair. RPE1's biology (p53-proficient, epithelial) may simply be unusually far
   from K562.
4. **There is no same-cell, different-study source.** "Cell identity" therefore cannot be separated from
   "same line in the same lab". VIPerturb-seq K562 (different lab) was not usable without conversion
   (audit §1.2).
5. **There is no technical replicate.** The platform/sequencing ceiling is unknown. The 0.82 includes
   biological (timepoint) and technical (library, batch) differences inseparably.
6. **Processing asymmetry.** GWPS was processed by us from raw counts with the scPertEval recipe. The target
   and the other scPertEval sources were processed by scPertEval, with upstream QC as deposited. This would
   be expected to *lower* C_GWPS, not raise it.

No claim is made that biological similarity causes transfer.

## 8. Paper decision

**REVISE AGAIN.**

**Why not CONTINUE.**

* The robust result, "an independent same-cell-line screen from the same lab transfers much better than other
  cell lines, beyond reliability", is largely expected. Nadig et al. 2025 and Replogle et al. 2022 already
  describe stronger within-line than cross-line agreement.
* As it stands, N5 confirms and quantifies that expectation rigorously. That is not yet a novel lead.
* The one thing that would make it novel is separating cell identity from study/lab compatibility. The
  current data cannot.

**Why not KILL.**

* The preregistered rule says the direction survives (Outcome A), and the effect is large and robust.
* Two observations are informative and not obviously in the literature in this noise-corrected form:
  1. the reliability-independent compatibility ordering (Jurkat > HepG2 > HCT116 > HEK293T ≈ RPE1 for a K562
     target, with the same-study, same-library RPE1 last);
  2. the about-one-third latent mismatch between two same-lab K562 screens.

**The revision:** reframe from "biological partner selection" to "**how much of cross-context transfer error
is cell line vs screen?**" with noise-corrected compatibility as the instrument. That question is only
answerable with the missing ladder level.

**Single most informative next experiment.** A preregistered **same-cell, different-lab** source test:
VIPerturb-seq K562 (CC BY 4.0) → K562 essential.

* Use the frozen N5 endpoint, with GWPS and Jurkat as fixed comparators.
* Prerequisites:
  * a Seurat-to-AnnData conversion (R environment) with provenance;
  * a predeclared minimum reliable-energy criterion, because VIPerturb is shallow (median 35 cells on the
    panel). The depth-matched GWPS result shows the estimator tolerates about 2.5× depth loss, but VIPerturb
    is shallower still.
* Decision logic, fixed in advance:
  * C_VIPerturb ≈ C_GWPS (about 0.8): cell identity dominates, and the direction is promoted.
  * C_VIPerturb at non-K562 levels (≤ 0.56): the GWPS advantage was largely study/lab compatibility. Pivot to
    cross-screen harmonisation.
  * In between: a quantified split.

## 9. Critical review: five strongest objections

1. **"This is expected, not discovered."** Two K562 screens agreeing better than K562 and RPE1 is the null
   expectation of any biologist.
   * *Needed:* show the comparison that is *not* expected, namely same cell / different lab against different
     cell / same lab (VIPerturb).
   * Also show that the noise-corrected magnitudes add something beyond published cross-line correlations.
2. **Same-lab technical sharing.** C_GWPS may be inflated by K562-specific technical signatures shared within
   one lab: reagent lots, cell stock drift, a 10x chemistry version, processing.
   * *Needed:* a different-lab K562 source; ideally a technical re-measurement pair to bound platform effects.
3. **One target, one same-cell source.** n = 1 cell line pair, with no replication on a second target.
   * *Needed:* the mirrored design (e.g. an RPE1 target with an independent RPE1 source of adequate
     reliability). Kaden is too shallow on this panel, so a new or alternative dataset is required.
4. **Validity of the noise correction.** The latent cosine assumes independent noise between split halves.
   Shared within-screen structure (gem-group batches, shared controls across a screen's halves, cell-cycle or
   UMI-depth programs) could bias the reliable-energy denominators.
   * The 0.822 vs 0.818 full/depth agreement is reassuring but internal.
   * *Needed:* batch-blocked split halves (by gem group), a check on synthetic data with injected shared batch
     effects, and agreement of raw and corrected rankings. The rankings agree here.
5. **The essential-gene panel and energy weighting.** The pooled metric is energy-weighted: the top 10 % of
   perturbations give 58 % of the difference. Essential-gene knockdowns produce strong, partly generic stress
   and lineage programs (e.g. erythroid or proliferation programs in K562).
   * The "same-cell advantage" may therefore be about generic K562 response programs, not perturbation-specific
     transfer. The shared on-target knockdown contributes to every source.
   * *Needed:* a predeclared analysis of the component orthogonal to each line's generic response axis
     (N1/N3 machinery); exclusion of the on-target gene; a non-essential, genome-wide panel where GWPS coverage
     allows.
