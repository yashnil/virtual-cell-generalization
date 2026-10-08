# N6 preregistered protocol: does an independent-lab K562 screen retain the transfer advantage?

Written 2026-10-07, after the coverage audit (`reports/n6_coverage_audit.md`, **PARTIAL (strong)**) and the
conversion audit (`reports/viperturb_conversion_audit.md`).

**State before any N6 outcome existed:**

* No VIPerturb response vector had been built.
* No N6 compatibility quantity had been computed.
* `outputs/n6` did not exist.

The SHA-256 of this file, the N6 code, the frozen axes and the converted artefact are recorded in
`data/provenance/research_v3/n6_protocol_digest.txt`. The pipeline refuses to run if any of them change.

**Frozen inputs:**

* **N5** is frozen and read only (`n5_freeze_sha256.txt`; `outputs/n5` is read-only).
* **N3** objects are read only.

**Hard constraints:**

* No new ML model, embeddings, foundation-model features, active selection or hyperparameter optimisation.
* No comparator, perturbation or gene chosen after outcomes.

> **Pre-run correction (2026-10-07, before any N6 outcome existed).** The first axis freeze applied the alias
> map to the measured-gene axis as well as to perturbation labels. That wrongly dropped GET1, GET3 and MICOS10,
> which VIPerturb's probe axis lists under their *current* symbols; only its perturbation labels use the older
> ones. The rule was corrected to "current symbol if present, else alias", and the gene axis went from 6,080 to
> 6,083. The panel (637) is unchanged. The protocol was re-hashed after this correction, with no outcome
> computed.

## 0. Question and identification

**Primary question.** Does VIPerturb-seq K562 (Satija lab; dCas9-KRAB-MeCP2; 3 single-sgRNA library; 10x Flex
fixed cells; Ultima) transfer into the Replogle K562 essential target substantially better than non-K562
sources?

**Identifies.** Whether the N5 same-cell advantage survives a screen that is independent in lab, library,
effector, chemistry and sequencer **together**.

**Does not identify.** Which of those technical axes causes any loss. They change jointly. Results are stated
as "same-cell compatibility across an independent screen". A loss is attributed to the joint "lab/platform
bundle", never to a single factor.

## 1. Data and frozen inclusion rule

**Axes** (`data/splits/n6_k562/`, frozen; built by `make_n6_axes.py` from identifiers and QC counts only):

* **Perturbations.** The N5 panel (1,054) with ≥ **30** VIPerturb cells passing the scPertEval cell filter
  (≥ 200 detected genes), pooled over guides. The VIPerturb label comes through the frozen alias map {GET1 →
  WRB, GET3 → ASNA1, MICOS10 → MINOS1}. Result: **637 perturbations**.
* **Genes.** N5 genes (6,408) present on the VIPerturb probe axis: **6,083 genes**. The current symbol is
  used if present; otherwise the frozen alias.
* **Rationale for the thresholds.** The 30-cell floor and the 200-gene cell filter are the values used
  unchanged throughout N1–N5. They were not chosen for N6.
* **Reliability enters through the preregistered gate (§4) and a secondary panel (§3), not the primary
  inclusion.** A per-perturbation reliability filter selects on the source's own noise. It is reserved for the
  secondary analysis.

**Sources and target** (`build_n6.py`):

| role | name | source |
|---|---|---|
| target | K562 essential | N3 object, subset |
| primary source | **VIPerturb_K562** | converted h5ad. Cells with ≥ 200 genes; `log1p(1e4·count / Σ_{19,068 genes} count)`; restricted to the axis; canonical delta against all 7,949 `NO-TARGET` cells; disjoint F/E1/E2 parts with the N1 algorithm, 5 repeats, seed `[20261006, 30, r]` |
| fixed comparator (positive reference) | **K562_GWPS** | frozen N5 GWPS object, subset |
| fixed comparator | **RPE1** | N3 object, subset |
| fixed comparator (primary, non-K562; N5's predeclared matched comparator) | **Jurkat** | N3 object, subset |
| positive control at VIPerturb depth | **K562_GWPS_vipdepth** | GWPS cells subsampled per perturbation to min(n_GWPS, n_VIP), seed `[20261006, 31, p]`; parts seed `[20261006, 32, r]`; all controls; N5 recipe |

**Build gates.** G3: median on-target canonical delta < 0, for VIPerturb and for the positive control. G4: part
sizes asserted.

## 2. Primary endpoint (frozen N5 endpoint, unchanged)

The **pooled both-sided noise-corrected latent cosine**, exactly as N5 §2 (`source_compatibility`):

```
C_S = Σ_p⟨S̃_p, T̃_p⟩ / sqrt(Σ_p E_r⟨S̃a_p, S̃b_p⟩ · Σ_p E_r⟨T̃a_p, T̃b_p⟩)
```

* Centring is over the **637 N6 perturbations**.
* Uncertainty: paired perturbation bootstrap, 2,000 resamples, seed `[20261006, 41]`.

## 3. Contrasts and secondary analyses

**Primary contrast.** C_VIP − C_Jurkat.

**Position statistic.** `f = (C_VIP − C_Jurkat) / (C_GWPS − C_Jurkat)`, computed per bootstrap draw. 0 means
VIPerturb behaves like the non-K562 comparator; 1 means it matches the same-study K562 screen.

**Also reported:** C_VIP − C_GWPS, C_VIP − C_RPE1.

**Reliability controls** (N5 logic):

* **Raw:** R1 median per-perturbation Pearson; R2 pooled raw cosine; R3 directional accuracy.
* **Adj-2.** Perturbation fixed-effect regression of the per-perturbation Pearson on source indicators
  {VIP, GWPS, RPE1, Jurkat}, with source split-half reliability (linear and squared) and log source magnitude.
  Reference: Jurkat. Cluster bootstrap 2,000, seed `[20261006, 43]`.
* **Adj-3.** Reliability-matched perturbations, |ρ_VIP − ρ_Jurkat| ≤ 0.05; median r_VIP − r_Jurkat; bootstrap
  seed `[20261006, 44]`.
* **Adj-4 (depth).** The positive control K562_GWPS_vipdepth: the same-study K562 source at VIPerturb's cell
  counts.

**Secondary panel B (not decisive).** Perturbations with VIPerturb per-perturbation split-half reliability
≥ 0.10, identical for every source. Report C_S, C_VIP − C_Jurkat and f (seed `[20261006, 42]`).

## 4. Reliability gate (Outcome D), fixed now

N6 is **uninformative (Outcome D)** if any of the following holds:

| id | condition |
|---|---|
| d1 | VIPerturb pooled split-half reliability (mean over repeats of the pooled Pearson of the centred halves) < **0.10** |
| d2 | the bootstrap 2.5th percentile of VIPerturb's summed reliable energy ≤ 0 |
| d3 | the positive control is unstable at VIPerturb depth: \|C_GWPS_vipdepth − C_GWPS\| > **0.10**, or the width of C_GWPS_vipdepth's 95 % interval > **0.15** |
| d4 | the panel has < **300** perturbations |

The thresholds were set from N5's observed stability (depth-matched GWPS moved C by 0.005; interval widths
≈ 0.03) and before any VIPerturb quantity existed.

The source-only gate quantities are written to `outputs/n6/source_gate.json` **before the target is read**.

## 5. Outcomes (mechanical; `analyse_n6.classify`)

* If the gate fails: **D**.
* Else, if the C_GWPS − C_Jurkat interval includes 0: **inconclusive** (the reference is not separated).
* Else:

| outcome | criterion | interpretation |
|---|---|---|
| **A: VIPerturb approaches GWPS** | C_VIP − C_Jurkat interval > 0, **and** the Adj-2 β_VIP − β_Jurkat interval > 0, **and** the f 2.5th percentile ≥ **0.67** | The same-cell transfer advantage survives independent lab/platform variation. It strengthens, without proving, the view that cellular context is a major determinant of transferability. Lab effects are not claimed absent |
| **B: VIPerturb resembles non-K562 sources** | the f 97.5th percentile ≤ **0.33** | The N5 advantage was substantially tied to study/lab compatibility. Pivot toward cross-screen technical variation and its contribution to apparent context specificity |
| **C: intermediate** | C_VIP − C_Jurkat interval > 0, **and** the Adj-2 interval > 0, but neither A nor B | Both cell identity and study/platform compatibility contribute materially. Next: variance/error decomposition across biological and technical axes |
| otherwise | — | **inconclusive**: not robust to reliability adjustment, or the interval spans the B/C boundary |

The thresholds 0.67 and 0.33 split the comparator-to-reference interval into thirds. They are fixed now.

## 6. Leakage controls

* Nothing is fitted to the target.
* The axes, the alias map, the positive-control subsampling counts (VIPerturb cell counts) and the gate use
  only identifiers and source data. The gate is written before the target is loaded.
* VIPerturb and the target are separate experiments from separate labs.
* Tests: `tests/test_n6.py` covers the outcome rules on synthetic bootstraps and the VIPerturb builder against a
  manual computation. `tests/test_n5.py` covers the endpoint's planted-value recovery.

## 7. Outputs

* `outputs/n6/` (git-ignored): the data object and gates, `source_gate.json`, per-perturbation terms,
  bootstrap tables, `n6_summary.json`, `n6_decision.json`, a manifest.
* `reports/n6_results.md`, written after the run.
* Stop after N6.

## 8. Execution

```
bash scripts/research_v3/run_n6_pipeline.sh
```

It verifies the digest, then runs `build_n6.py`, `run_n6.py` and `analyse_n6.py`.
