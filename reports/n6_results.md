# N6 results: independent-lab K562 source (VIPerturb-seq), Outcome D

Date: 2026-10-07. Protocol: `reports/n6_protocol.md` (SHA-256 `86a26384…`). Audits:
`reports/n6_coverage_audit.md`, `reports/viperturb_conversion_audit.md`.

> **Status: Outcome D (uninformative), by the preregistered reliability gate.** Gate item d1 failed:
> VIPerturb's pooled split-half reliability is 0.0866, below the frozen threshold of 0.10.
>
> **No compatibility comparison in this report can be interpreted.** This covers every C_S value, every
> contrast, f, Adj-2, Adj-3 and panel B. The numbers are kept in full for the record. They are **not** evidence
> for or against a same-cell advantage across studies. They do not support any statement about VIPerturb
> versus Jurkat, RPE1 or GWPS, in either direction.

**Integrity, verified before writing:**

* All 23 entries in `data/provenance/research_v3/n6_protocol_digest.txt` re-hash OK. These are the protocol,
  audits, N6 code, axes, conversion files, the N5 freeze file and the converted h5ad (`04ef7c27…`).
* All 16 entries in `outputs/n6/manifest_sha256.txt` re-hash OK. The manifest covers every file in
  `outputs/n6`.
* All 8 entries in `outputs/n6/data/sha256.txt` re-hash OK.
* `code_sha256` in `n6_decision.json` matches the digest for the protocol, `build_n6.py`, `run_n6.py`,
  `analyse_n6.py` and `source_compatibility.py`.
* The N5 freeze (`n5_freeze_sha256.txt`, 33 entries) and the N5 protocol digest re-hash OK.
* N6 was not re-run, tuned or modified. The tables come from `outputs/n6/{source_gate.json, n6_summary.json,
  n6_decision.json, data/build_gates.json}`. §4 is the only new computation. It is exploratory, read-only and
  labelled as such.

**Build gates** (`outputs/n6/data/build_gates.json`):

| gate | VIPerturb_K562 | K562_GWPS_vipdepth |
|---|---|---|
| G3: median on-target canonical delta (n = 470) | −0.281 (pass) | −0.259 (pass) |
| cells | 39,939 selected and 39,939 kept; 7,949 `NO-TARGET` controls | — |
| median cells per perturbation | 47 | 46 |

G4: part sizes asserted (the run completed).

**Attribution.**

* VIPerturb-seq: Bradu et al., bioRxiv 2026.02.12.705613; Zenodo 10.5281/zenodo.18460279; CC BY 4.0.
  Changes: converted from Seurat, pseudobulk statistics derived.
* Replogle et al. 2022, *Cell* (figshare+ 10.25452/figshare.plus.20029387); CC BY 4.0.

---

## 1. Reliability gate (Outcome D), decisive

The source-only quantities were written to `outputs/n6/source_gate.json` before the target was read.

| id | condition (frozen) | observed | fails? |
|---|---|---|---|
| **d1** | VIPerturb pooled split-half reliability < 0.10 | **0.0866** | **YES** |
| d2 | 2.5th percentile of VIPerturb summed reliable energy ≤ 0 | 2,481 (point 2,794) | no |
| d3 | \|C_GWPS_vipdepth − C_GWPS\| > 0.10, or C_GWPS_vipdepth 95 % interval width > 0.15 | \|0.811 − 0.818\| = **0.007**; width **0.049** | no |
| d4 | panel < 300 perturbations | 637 | no |

**Mechanical outcome (`analyse_n6.classify`): D.** One failing item is enough under §4–5 of the protocol. The
other branches (A/B/C/inconclusive) are not evaluated.

### 1.1 The positive control passed: low cell count alone does not explain the failure

The GWPS-at-VIPerturb-depth positive control (d3) **passed**, with a large margin:

* **The endpoint holds up at this depth.** GWPS subsampled to VIPerturb's per-perturbation cell counts
  (median 46 vs 47) gives C = 0.811 [0.785, 0.834], against 0.818 [0.796, 0.838] at full depth. That is a
  shift of 0.007, far inside the 0.10 tolerance.
* **The two have similar per-perturbation reliability.** The medians are 0.052 for GWPS_vipdepth and 0.055
  for VIPerturb.

So the same-study K562 screen, at the same cell counts and with similar per-perturbation reliability, still
gives a stable, high C_S. **Low cell count alone does not explain why VIPerturb failed the reliability gate,
and it would not, on its own, have stopped the endpoint from working.** What separates VIPerturb from the
positive control at equal depth is a property of the VIPerturb data. That property is not identified here.

**Caveat (exploratory, §4).** At matched depth, GWPS's pooled reliability is 0.109. That clears the d1
threshold, but only narrowly. Depth therefore accounts for most of the gap between full-depth GWPS (0.341) and
VIPerturb (0.087). The residual VIPerturb shortfall at equal depth is real but modest. This qualifies the
statement above without contradicting it.

## 2. Full numerical outputs: NON-INTERPRETABLE under the preregistered protocol

> Every quantity in this section was computed by the frozen pipeline after the gate had already failed. Under
> §5 of the protocol, Outcome D makes them **non-interpretable**. They are listed only for completeness and
> auditability. No contrast, ranking or sign below may be cited as a finding.

### 2.1 Primary endpoint C_S (637 perturbations × 6,083 genes; paired bootstrap, 2,000 draws)

| source | point | bootstrap median | 95 % interval |
|---|---|---|---|
| VIPerturb_K562 *(non-interpretable)* | 0.285 | 0.285 | [0.254, 0.317] |
| K562_GWPS | 0.818 | 0.818 | [0.796, 0.838] |
| K562_GWPS_vipdepth (positive control) | 0.811 | 0.811 | [0.785, 0.834] |
| RPE1 | 0.316 | 0.316 | [0.286, 0.346] |
| Jurkat | 0.522 | 0.521 | [0.496, 0.547] |

### 2.2 Contrasts *(all non-interpretable)*

| contrast | median | 95 % interval |
|---|---|---|
| C_VIP − C_Jurkat (preregistered primary) | −0.237 | [−0.270, −0.204] |
| C_GWPS − C_Jurkat (reference separation) | +0.296 | [+0.261, +0.329] |
| f = (C_VIP − C_Jurkat)/(C_GWPS − C_Jurkat) | −0.802 | [−0.997, −0.640] |
| C_VIP − C_GWPS | −0.533 | [−0.565, −0.497] |
| C_VIP − C_RPE1 | −0.031 | [−0.070, +0.007] |

f lies outside [0, 1]. The protocol's scale for f assumes VIPerturb lies between the comparator and the
reference, and a value outside that range is one reason not to read f under gate failure.

### 2.3 Reliability-adjusted analyses *(non-interpretable)*

**Adj-2.** Perturbation fixed effects; split-half reliability (linear and squared) and log magnitude as
covariates; reference Jurkat; cluster bootstrap 2,000.

| source | β − β_Jurkat | 95 % interval |
|---|---|---|
| VIPerturb_K562 | −0.033 | [−0.043, −0.022] |
| K562_GWPS | +0.115 | [+0.101, +0.130] |
| RPE1 | −0.113 | [−0.127, −0.101] |

**Adj-3.** Reliability-matched pairs, \|ρ_VIP − ρ_Jurkat\| ≤ 0.05, n = 149. The median r_VIP − r_Jurkat is
−0.037 [−0.048, −0.026].

### 2.4 Secondary panel B *(not decisive by design; also non-interpretable)*

Panel B is the 137 perturbations with VIPerturb per-perturbation reliability ≥ 0.10.

| source | C_S |
|---|---|
| VIPerturb_K562 | 0.409 |
| K562_GWPS | 0.831 |
| K562_GWPS_vipdepth | 0.827 |
| RPE1 | 0.382 |
| Jurkat | 0.570 |

* C_VIP − C_Jurkat: −0.161 [−0.202, −0.115].
* f: −0.620 [−0.853, −0.393].

Panel B selects on VIPerturb's own noise and cannot rescue the gate (protocol §1, §3).

### 2.5 Raw reliability controls *(descriptive)*

| source | R1 median per-perturbation Pearson | R2 pooled raw cosine | R3 directional accuracy | median split-half reliability | median cells |
|---|---|---|---|---|---|
| VIPerturb_K562 | 0.040 | 0.073 | 0.526 | 0.055 | 47 |
| K562_GWPS | 0.230 | 0.371 | 0.597 | 0.212 | 255 |
| K562_GWPS_vipdepth | 0.118 | 0.231 | 0.560 | 0.052 | 46 |
| RPE1 | 0.118 | 0.166 | 0.544 | 0.453 | 105 |
| Jurkat | 0.124 | 0.225 | 0.553 | 0.197 | 112 |

## 3. What N6 does and does not establish

**Does establish (from gate and control quantities only):**

* The converted VIPerturb data carry an on-target knockdown signal: G3 gives −0.281, and the reliable energy
  is clearly positive (d2).
* On the N6 panel, VIPerturb's pooled split-half reliability falls below the preregistered floor.
* The endpoint itself is stable at VIPerturb's cell counts, as the GWPS positive control shows.

**Does not establish:**

* Whether an independent-lab K562 screen keeps or loses the N5 same-cell transfer advantage. The question N6
  was designed to answer is **unanswered**, not answered negatively.
* Any attribution to lab, library, effector, Flex chemistry or sequencer.

**N5 is unchanged.** It still shows within-study source compatibility, with cell identity not separable from
same-lab compatibility.

## 4. Exploratory post-hoc diagnostic (not preregistered; read-only)

This diagnostic recomputes the d1 statistic for the K562 sources from the saved N6 and N5 part means. It uses
the identical formula in `run_n6.py`: the mean over repeats of the pooled Pearson of the centred split halves.
It was run from the session scratchpad; nothing in `outputs/n6` or the repository was written. The VIPerturb
value reproduces the gate value exactly.

| source | pooled split-half reliability | reliable / total centred energy | median centred response norm |
|---|---|---|---|
| VIPerturb_K562 | 0.0866 (= gate) | 0.163 | 5.00 |
| K562_GWPS_vipdepth | 0.109 | 0.201 | 4.93 |
| K562_GWPS (full depth) | 0.341 | 0.510 | 2.55 |

**Further diagnostics:**

* Fraction of perturbations with per-perturbation reliability ≥ 0.10: VIPerturb 0.215, GWPS_vipdepth 0.267,
  GWPS 0.824, Jurkat 0.779, RPE1 0.975.
* Spearman correlation between VIPerturb and GWPS_vipdepth per-perturbation reliability: 0.35.

**Reading (exploratory):**

* At equal depth, VIPerturb keeps about 20 % less reproducible energy than GWPS.
* By a rough Spearman–Brown argument, VIPerturb would need about 17 % more cells per perturbation (median
  ≈ 55) to reach 0.10.
* d1 therefore failed by a modest margin.
* This must **not** be used to re-open N6. Lowering the threshold, re-panelling or adding depth after the
  outcome would be a post-hoc protocol change.

## 5. Can another independent K562 Perturb-seq dataset answer the question?

**Requirements.** To answer the same-cell / different-study question under the frozen N5/N6 endpoint, a
candidate needs all of the following:

* K562.
* A different lab or study from Replogle/Weissman.
* Single-perturbation (low-MOI) design.
* Transcriptome-wide readout covering most of the 6,083–6,408 axis genes.
* At least 300 panel perturbations with at least 30 cells.
* Expected pooled reliability ≥ 0.10.
* A licence permitting research use.

The search below used the literature and preprints; no data was downloaded.

| candidate | lab / design | overlap with the essential panel | reliability prospects | verdict |
|---|---|---|---|---|
| **VIPerturb-seq K562** (used) | Satija; KRAB-MeCP2; Flex; Ultima | 637 at ≥ 30 cells | d1 failed (0.087) | exhausted under this protocol |
| **Gasperini et al. 2019** (GSE120861) | Shendure; CRISPRi of ~5,900 enhancers; **high MOI (~28 guides/cell)** | Only TSS positive controls (on the order of a few hundred pairs; exact count unverified). Highly expressed genes, overlap with essential panel unknown | Many cells per guide, but every cell carries ~28 perturbations; per-perturbation deltas need regression and are not comparable to the canonical delta | **No**: design incompatible with the endpoint; overlap probably < 300 |
| **Jiang, Dalgarno et al. 2025** (*Nat Cell Biol*) | Satija; KRAB-MeCP2; 6 lines incl. K562; signalling regulators under stimulation | >1,500 perturbations over six lines and five pathway contexts. The K562 share is unverified and enriched for pathway genes, not essential genes | Unknown; stimulation contexts add a condition axis | **Unlikely**: overlap probably far below 300. Same lab and effector as VIPerturb, so it is not independent of VIPerturb |
| **K562 CRISPRi CROP-seq preprint** (bioRxiv 2025.10.23.684112) | Pooled CRISPRi | 377 sgRNAs, of which 56 knock down efficiently, per the preprint | Low | **No** |
| **Adamson 2016, Norman 2019, Jost 2020, Replogle 2020** | Weissman lab (the target's lab); UPR / CRISPRa / titration | Small panels, or CRISPRa | — | **No**: same lab as the target, so not independent |
| **Dixit et al. 2016** | Regev; Cas9 knockout | Tens of transcription factors | — | **No**: too few, and a different modality |
| **X-Atlas/Orion; X-Atlas/Pisces** (Xaira) | FiCS CRISPRi | No K562 context (Orion: HCT116, HEK293T; Pisces: HCT116, HEK293T, HepG2, iPSC, Jurkat ±stim, iPSC multi-differentiation) | — | **No K562** |

**Conclusion: no.** No public K562 Perturb-seq dataset was identified that meets the bar: an independent study,
low MOI, ≥ 300 panel perturbations at ≥ 30 cells, and plausible reliability ≥ 0.10. VIPerturb was the only
candidate that cleared the overlap bar, and it failed reliability.

**Caveats:**

* The search was literature and preprint based, not an exhaustive scan of GEO/SRA.
* Several figures (Gasperini TSS-control count, Jiang K562 share) are unverified and flagged as such.

**What this means.** With current public data, the same-cell / different-study rung of the ladder cannot be
measured. The paper should state the cell-identity vs study-compatibility confound as an explicit limitation,
not an open result.

## 6. Recommendation (not implemented)

1. **Do not** re-run, re-threshold or re-panel N6. Outcome D stands.
2. **Next experiment, if any: "N7-desk".** This is a short, preregistered, metadata-only go/no-go audit. Its
   aim is to close the same-cell/different-study rung formally, with no new outcome data.
   * Check GEO/SRA/Zenodo for K562 single-cell CRISPR screens from 2023–2026 by any lab other than Weissman's.
     Pre-commit the N6 eligibility rules (§5 above) and a Spearman–Brown reliability projection, anchored on
     the GWPS_vipdepth reference.
   * Include GSE120861 and the Jiang 2025 K562 arm explicitly.
   * Go only if at least one candidate passes **every** rule on metadata alone.
   * The expected outcome is **no-go**. In that case, record the rung as unidentifiable with public data and
     move to paper writing on N3/N5. The claims would be within-study compatibility, the RPE1 ordering, and the
     calibration-budget limits.
3. **Alternative.** A new K562 screen generated for the purpose is the only route to a clean same-cell /
   different-study test. Matching the effector and panel to Replogle while changing lab and platform would be
   needed. This is out of scope for the current project.
