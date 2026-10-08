# N5 preregistered protocol: source compatibility within and across studies (K562 target)

Written 2026-10-07, after the read-only coverage audit (`reports/n5_coverage_audit.md`, verdict **PARTIAL**).

**Before any N5 transfer quantity existed:**

* No GWPS response vector had been built.
* No source–target similarity had been computed.

The SHA-256 of this file, and of the N5 code it names, is recorded in
`data/provenance/research_v3/n5_protocol_digest.txt` before anything is run.

After that, nothing below may change. Implementation bugs found while running are fixed and reported as
deviations, with their consequence. Outcome criteria are never changed.

**Hard constraints:**

* No new ML model; no embeddings, neural networks, GNNs, foundation models or learned source selectors.
* No perturbation, gene or source chosen after seeing transfer outcomes.
* No frozen N1/N3 result modified.

---

## 0. Question and what this design can and cannot identify

**Primary question.** For the fixed K562 essential-screen target: does an independent **K562** screen (GWPS)
predict held-out K562 perturbation responses better than a **non-K562** source, after controlling for source
measurement reliability and signal strength?

**Identifiable.** A same-cell vs different-cell contrast **within one study and lab** (Replogle 2022):

* B = K562 GWPS (same cell; different screen, library scale and timepoint ~d8);
* C = RPE1 essential (different cell; same essential library as the target, ~d7);
* target = K562 essential (~d6).

The library and timepoint differences **penalise B**. An advantage of B over C therefore cannot be
attributed to library or timepoint compatibility.

**Not identifiable** (audit §3):

* a technical-platform ceiling (no Illumina/Ultima pair);
* same cell / different study vs different cell / same study (no usable K562 source from another lab).

Results are therefore described as **within-study source compatibility**, never as "cell identity causes
transfer" in general.

## 1. Data object

**Axes.** `data/splits/n5_k562/` (frozen, hashed): **1,054 perturbations × 6,408 genes**. These are the N3
panel ∩ GWPS ≥ 30 cells, and the N3 genes ∩ unique GWPS symbols. Order follows the N3 files.

**Seven contexts:**

| index | context |
|---|---|
| T | K562 essential (target) |
| 1 | K562 GWPS |
| 2 | RPE1 |
| 3 | HepG2 |
| 4 | Jurkat |
| 5 | HCT116 |
| 6 | HEK293T |

**Sources other than GWPS, and the target.**

* Subset exactly from the N3 object `outputs/n3/data`: canonical deltas and the disjoint F/E1/E2 part means
  (perturbed and control cells), 5 repeats.
* Subsetting a mean of per-gene values is exact.

**GWPS** (`scripts/research_v3/build_n5_gwps.py`):

* **Cells:** all cells whose `obs.gene` is an N5 perturbation (all transcripts and guide pairs pooled, matching
  scPertEval's suffix-stripped labels), plus all `non-targeting` cells. Drop cells with < 200 non-zero genes
  on the 8,248-gene axis (scPertEval `filter_cells(min_genes=200)`).
* **Value:** `log1p(1e4 · count / Σ_{8,248 genes} count)`, then restricted to the 6,408 genes. This is the
  scPertEval recipe.
* **Canonical delta:** mean over all of a perturbation's cells − mean over all controls (pooled over gem
  groups, as scPertEval).
* **Disjoint parts:** the N1 three-way algorithm.
  * Within each perturbation, cells are sorted by row index, then shuffled.
  * Parts are F = ⌊n/2⌋, E1 = E2 = ⌊n/4⌋.
  * Controls are split the same way, before any delta exists.
  * 5 repeats, seed `[20261006, 6, r]`.
* **Depth-matched variant** (Adj-4). Per perturbation, a random subset of exactly `n_RPE1,p` cells is drawn
  first (seed `[20261006, 7, p_index]`), where `n_RPE1,p` is the RPE1 canonical cell count. Its canonical delta
  and F/E1/E2 parts (5 repeats, seed `[20261006, 8, r]`) are then built the same way, against the full control
  set split as above. If `n_RPE1,p` > the GWPS count, all GWPS cells are used, and this is reported.

**GWPS gates (stop on failure):**

* G2: the canonical delta equals the count-weighted part combination where no cell is dropped (1e-6 relative).
* G3: the median on-target canonical delta is < 0 among perturbations on the gene axis.
* G4: part sizes and disjointness are asserted.

## 2. Primary endpoint (one)

**Pooled reliability-corrected latent cosine.** For each source S, define

```
C_S = Σ_p ⟨S̃_p, T̃_p⟩  /  sqrt( Σ_p E_r⟨S̃a_p, S̃b_p⟩ · Σ_p E_r⟨T̃a_p, T̃b_p⟩ )
```

* `S̃`, `T̃` are the canonical deltas, each **centred over the 1,054 perturbations** (template removed).
* `a`/`b` are the F part and the (E1+E2)/2 part. These are disjoint in perturbed and control cells and centred
  the same way.
* `E_r` is the mean over the 5 repeats.
* The numerator is unbiased for the latent cross product, because S and T are independent experiments. The
  denominator terms are the project's unbiased reliable energies (`transferability.signal_energy` logic).
* C_S estimates the cosine between the noise-free source and target response fields. It is interpretable (1 =
  identical latent responses), robust to the template, comparable across screens, and corrected for
  measurement reliability on **both** sides.

**Uncertainty.**

* Paired perturbation bootstrap: 2,000 resamples of the 1,054 perturbations, the **same resample for every
  source**, seed `[20261006, 9]`.
* Each resample recomputes the three sums from the per-perturbation terms.
* Report the 95 % percentile interval for every C_S and for every contrast.

## 3. Secondary endpoints

| id | endpoint |
|---|---|
| R1 (raw) | median over perturbations of Pearson `r(S̃_p, T̃_p)` across genes |
| R2 (raw) | pooled raw cosine `Σ⟨S̃,T̃⟩ / sqrt(Σ‖S̃‖² Σ‖T̃‖²)` (no reliability correction) |
| R3 | directional accuracy. Per perturbation, the 200 genes with the largest \|T̃a\| (target F part); the fraction where `sign(S̃) = sign(T̃b)` (target E part); median over perturbations |
| R4 | calibration (N3 machinery, single source). New test split: 30 % of 1,054 = 316, seed `[20261006, 21]`; pool 738. E4 with source set {S} only, at k ∈ {20, 50}, 40 draws (8 per repeat × 5), anchors seeded `[20261006, 22, k, r, d]`. Report own-frame M3 (γ⊥ relative to S's own consensus) and M1. Anchors use the target F part; truth uses the target E parts |

## 4. Reliability and signal-strength control (fixed now)

**Adj-1.** Built into C_S (both-sided disattenuation).

**Adj-2 (perturbation-level regression).**

* Model: `y_{p,S} = α_p + β_S + γ1·ρ_{S,p} + γ2·ρ_{S,p}² + δ·log‖S̃_p‖ + ε`, by OLS over all (p, S), where:
  * `y_{p,S}` is R1's per-perturbation Pearson;
  * `ρ_{S,p}` = E_r Pearson(S̃a_p, S̃b_p), the source split-half reliability;
  * perturbation fixed effects α_p absorb target reliability and target signal strength.
* Report the contrasts β_S − β_RPE1, with a perturbation-cluster bootstrap (2,000, seed `[20261006, 10]`).

**Adj-3 (reliability-matched perturbations, B vs C).**

* Perturbations with `|ρ_GWPS,p − ρ_RPE1,p| ≤ 0.05`.
* Statistic: the median of `r_GWPS,p − r_RPE1,p`, with a bootstrap 95 % interval (seed `[20261006, 11]`).
* Report n.

**Adj-4 (design-based depth matching).** C_GWPS, R1 and Adj-2 recomputed with the depth-matched GWPS
(n_GWPS,p = n_RPE1,p).

**Predeclared matched non-K562 comparator** (rule computed from source data only, before any target
quantity):

* The non-K562 source whose median over perturbations of `ρ_{S,p}` is closest to GWPS's.
* Ties are broken by the closest median cell count.
* The rule's output is written to `outputs/n5/matched_comparator.json` before the target is read.

## 5. Contrasts, ladder and outcome criteria

**Primary contrast.** Δ_BC = C_GWPS − C_RPE1, paired bootstrap.

**Secondary contrasts.** C_GWPS − C_matched; C_GWPS − mean of the 5 non-K562 C_S.

**Ladder (N5-D).** Tested as ordered adjacent differences, each with its paired 95 % interval:

```
C_GWPS  >  C_RPE1  >  mean(C_HepG2, C_Jurkat)  >  mean(C_HCT116, C_HEK293T)
```

"Ladder holds" iff all three adjacent intervals are > 0. Otherwise report which steps hold.

**Outcomes.** These are fixed now and may co-occur, except that A, B and C′ are mutually exclusive by
construction.

| outcome | criterion | interpretation |
|---|---|---|
| **A: same-cell advantage survives reliability control** | Δ_BC interval > 0, **and** the Adj-2 contrast β_GWPS − β_RPE1 interval > 0, **and** the Adj-4 (depth-matched) Δ_BC interval > 0 | Within one study, a same-cell source carries transfer information beyond source reliability, depth and signal strength, despite a library and timepoint mismatch that favours the comparator. Supports source-context selection as a research direction. Not a general causal claim |
| **B: advantage disappears after control** | the raw R1 difference (median r_GWPS − r_RPE1, bootstrap seed `[20261006, 12]`) has interval > 0, **but** A fails | The basal-similarity / partner signal mainly reflected source quality and depth. Drop biological-partner selection as the main direction |
| **C′: screen/library compatibility dominates** | Δ_BC interval < 0 | RPE1, a different cell line sharing the target's library and study, beats a same-cell screen. Library/screen compatibility outweighs cell identity within this study. Pivot toward harmonisation and assay effects. (This replaces N5-C's same-study-vs-same-cell-other-study test, which is not identifiable) |
| **D′: same-cell cross-screen mismatch is already large** | C_GWPS upper interval < 0.70 | Even two K562 screens from the same lab share less than about half of their latent response variance (0.70² ≈ 0.49). A large share of apparent "context-specific" error may be screen-level. Quantify before biological claims |
| **E: no source structure** | max_S C_S − min_S C_S < 0.05, **and** the Δ_BC interval includes 0 | Neither cell match nor study explains transfer. Kill this direction |

If none of A, B, C′ or E holds: "inconclusive". In particular, Δ_BC whose interval includes 0 and a raw
difference that is not significant is inconclusive.

**Kill rule for the direction.** E holds, or B holds, or C′ holds. Each points away from biological source
selection.

## 6. Leakage controls

* The primary analysis fits nothing: every perturbation is evaluated, and no target value enters any source
  quantity.
* In R4, anchors enter only through the target F part, and truth only through the target E parts. This is the
  N1/N3 machinery, with its tests.
* GWPS and the target are separate experiments (separate cells and gem groups). No target cell can enter
  GWPS.
* The matched comparator, axes and depth-matching counts use only source data (RPE1 cell counts) and
  identifiers.
* **Test:** replacing the target with noise leaves every source-side quantity and the comparator choice
  bit-identical (`tests/test_n5.py`).

## 7. Outputs

* `outputs/n5/` (git-ignored): GWPS data object and gates, per-source and per-perturbation terms, bootstrap
  draws, contrasts, `n5_decision.json`, a SHA-256 manifest.
* `reports/n5_results.md`.
* Stop after N5. No acquisition model.

## 8. Execution

```
bash scripts/research_v3/run_n5_pipeline.sh
```

The script runs, in order:

1. `build_n5_gwps.py` (GWPS object + gates);
2. `run_n5.py` (endpoints, adjustments, R4);
3. `analyse_n5.py` (contrasts, outcomes).

It refuses to run if this protocol's hash differs from the recorded digest.
