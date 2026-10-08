# N5 coverage and provenance audit (read-only)

Date: 2026-10-07. Scope:

* Metadata, identifiers and storage layout only.
* **No N5 transfer quantity was computed or inspected.**
* The only expression values read: 3 cells × 2,000 genes of the GWPS matrix, to confirm integer raw counts
  (a format check).

Manifests:

* `data/splits/n5_k562/` (axes, exclusions, SHA-256 in `sha256.txt`);
* `data/provenance/research_v3/n5_source_metadata.json` (SHA-256 `823ed569…`).

## 1. Candidate datasets

"Per publication" means the fact comes from the cited paper and is **not encoded in any local metadata**.
"Not recorded" means it is neither local nor verified.

### 1.1 Datasets with the core comparisons

| | K562 essential (**target**) | K562 GWPS | RPE1 essential | HepG2 | Jurkat | HCT116 | HEK293T |
|---|---|---|---|---|---|---|---|
| local path | `data/raw/scperteval/replogle22k562_processed_complete.h5ad` | `data/raw/competition_v2/K562_gwps_raw_singlecell_01.h5ad` | `…/replogle22rpe1_processed_complete.h5ad` | `…/nadig25hepg2_…` | `…/nadig25jurkat_…` | `…/xatlas_orion/data/HCT116` | `…/HEK293T` |
| source | scPertEval of Replogle 2022 (figshare+ 20029387) | Replogle 2022, figshare+ 20029387, file 35775507 | scPertEval of Replogle 2022 | scPertEval of Nadig 2025 (GSE264667) | same | X-Atlas/Orion (Xaira), SLAF @598aa544 | same |
| study / lab | Replogle 2022 / Weissman | **same study** | **same study** | Nadig 2025 / Weissman (companion) | same | Xaira | Xaira |
| licence | CC BY 4.0 | CC BY 4.0 | CC BY 4.0 | public GEO deposit (as registered) | same | CC BY-NC-SA 4.0, research-only | same |
| assay | CRISPRi Perturb-seq, 10x 3′ | same | same | same | same | CRISPRi FiCS (fixed cells) | same |
| library | essential dual-sgRNA (per publication; no guide IDs locally) | genome-wide dual-sgRNA (`sgID_AB` verified) | essential (per publication) | essential protocol (per publication) | same | genome-wide dual-guide (`P1P2-1\|P1P2-2` verified) | same |
| timepoint | ~day 6 (per publication) | ~day 8 (per publication) | ~day 7 (per publication) | not recorded | not recorded | not recorded | not recorded |
| sequencing platform | not recorded | not recorded | not recorded | not recorded | not recorded | not recorded | not recorded |
| cells | 308,646 | 1,989,578 | 240,774 | 133,757 | 258,202 | 3,409,169 | 4,534,299 |
| perturbations (≥ 30 cells) | 1,971 | 9,675 (of 9,866) | 2,016 | 1,818 | 2,137 | 18,293 labels | 18,311 labels |
| controls | 10,691 | 75,328 (all 267 gem groups) | 11,485 | 4,976 | 12,013 | 165,777 | 218,838 |
| genes | 8,563 | 8,248 (Ensembl-indexed) | 8,749 | 9,623 | 8,881 | 38,606 | 38,606 |
| naming | symbol; suffixes stripped; `control` | `obs.gene` symbol + `transcript` (P1P2/P1/P2/ENST) | as target | as target | as target | `gene_target` symbol | same |
| raw counts | **no** (log1p CP10K) | **yes** (integer, float32 dense) | no | no | no | yes | yes |
| usable here | yes | yes | yes | yes | yes | yes (research-only memo) | yes |

### 1.2 Requested but not available

| requested | status |
|---|---|
| K562 KD8 Illumina / K562 KD8 Ultima (technical re-measurement pair) | **No local file and no provenance record.** Cannot be audited |
| A separate K562 essential "day 6" raw file | Only the scPertEval-processed version is local; the raw figshare file was not downloaded |
| Same cell line, different lab: VIPerturb-seq K562 | GREEN (CC BY 4.0) but **not downloaded**. Seurat `.rds` (R not installed). Manifest: 637 of the N3 panel have ≥ 30 cells, median 35 cells |
| Same cell line (RPE1), different lab: Kaden 2025 | Local, but only **103** panel perturbations, with per-perturbation split-half reliability 0.12–0.17 (`external_unseen_perturbation_validation_v1.md`) |

## 2. Overlap (identifiers only)

### 2.1 Perturbations

| set | count |
|---|---|
| N3 six-context panel | 1,062 |
| N3 panel ∩ GWPS (any cells) | 1,062 (100 %) |
| **N3 panel ∩ GWPS ≥ 30 cells: the N5 panel** | **1,054 (99.2 %)** |
| N3 perturbations below 30 GWPS cells (excluded) | 8: ARL2, DNAJC17, PDCD11, PLK1, RPS27A, SEM1, THOC5, WDR70 |
| four-context panel (1,264) ∩ GWPS ≥ 30 | 1,253 |

* All N5 perturbations are, by construction of the N3 panel, present with ≥ 30 cells in K562 essential, RPE1,
  HepG2, Jurkat, HCT116 and HEK293T.
* GWPS cells per N5 perturbation: median 223, minimum 30.
* **Naming.**
  * 0 symbol mismatches; no alias resolution needed. Every N3 symbol is an exact `obs.gene` label in GWPS.
  * 0 case-only matches.
* **Duplicate mappings.** 110 of the 1,054 have two or more TSS-specific guide pairs in GWPS (`transcript` ≠
  P1P2).
  * scPertEval's essential labels strip guide suffixes, so pooling all of a gene's cells is the matching
    definition.
  * Whether the essential screen also had multiple TSS guides for these genes is **not verifiable locally**.
* **Ambiguous symbols.** The GWPS gene axis has 2 duplicated symbols (HSPA14, TBCE). Both are excluded from the
  gene axis, and neither is an N5 perturbation.

### 2.2 Genes

* N3 genes (6,499) ∩ GWPS unique symbols: **6,408**. This is the N5 gene axis.
* The response vector definition can be reproduced **without changing the scientific object**:
  * scPertEval applied `normalize_total(1e4)` over each dataset's own filtered gene axis, then `log1p`,
    `filter_cells(min_genes=200)` and `filter_genes(min_cells=3)`.
  * GWPS raw counts on its 8,248-gene axis allow exactly this recipe: CP10K over the 8,248 genes, log1p, the
    same cell filter. Then restrict to the 6,408 genes.
  * The residual asymmetry: upstream guide assignment and cell QC are the depositors'. The target's QC was
    applied by scPertEval to the essential screen; GWPS QC is as deposited.

### 2.3 Reliability feasibility

| dataset | how reliability can be estimated |
|---|---|
| target and the five non-GWPS sources | Already split into disjoint F/E1/E2 parts, perturbed **and** control cells, 5 repeats (the N3 object). Subsetting to the N5 axes is exact |
| GWPS | Raw cells, ≥ 30 per perturbation (median 223), 75,328 controls spread over all 267 gem groups. The identical three-way disjoint split is feasible. Perturbations span a median of 149 gem groups, so a pooled control mean (as scPertEval) is appropriate |

Technical-replicate reliability (the same cells re-sequenced) is **not available** for any dataset.

## 3. Compatibility ladder: what is supported

Target: K562 essential.

| level | definition | available source(s) | supported? |
|---|---|---|---|
| A | same experiment, different technical measurement | none | **No** |
| B | same cell line, different screen | K562 GWPS. Same study and lab; differs in library scale, timepoint (~d8 vs ~d6) and batches | **Yes**, but "different screen" here is *within one study* |
| C | different cell line, same study | RPE1 essential. Same study, lab and essential library; ~d7 | **Yes** |
| C′ | different cell line, same lab, companion study | HepG2, Jurkat (Nadig 2025) | Yes (descriptive rung) |
| D | different cell line, different study/lab | HCT116, HEK293T (Xaira; FiCS; genome-wide library) | **Yes** |
| — | same cell line, different study/lab | none usable (VIPerturb not local or too shallow; Kaden too few and unreliable) | **No** |

**Which causal comparisons are supported:**

* **B vs C (primary).** This contrasts same cell vs different cell with study, lab and chemistry held
  fixed.
  * The remaining differences are timepoint and library scale.
  * **They penalise the same-cell source.** GWPS differs from the target in library and timepoint; RPE1 shares
    the target's essential library.
  * An advantage of GWPS over RPE1 therefore cannot be explained by library or timepoint compatibility.
  * An advantage of RPE1 over GWPS could be either cell-independent screen compatibility or chance. It would
    be uninterpretable as "biology does not matter".
* **B vs D, and C vs D (secondary).** These confound lab, chemistry, library and cell. They are ordinal
  descriptions, not causal contrasts.

**Not supported:**

* N5-A (technical ceiling). The within-target split-half reliability is the only ceiling available; it is
  measurement noise, not platform disagreement.
* The decisive N5-C contrast (same cell / different study vs different cell / same study). It needs a K562
  source from another study. None is usable locally.

## 4. Feasibility verdict

**PARTIAL (strong).**

* Levels B, C, C′ and D are supported on **1,054** shared perturbations × **6,408** genes, with
  disjoint-split reliability estimable for every source.
* Level A and "same cell / different study" are not supported.

**Most serious confound.** The only same-cell source (GWPS) comes from the **same study and lab** as the target.

* Cell-identity effects can be identified *within a study*.
* "Biological match" cannot be separated from "study/lab compatibility" in the sense of N5-C.
* The library/timepoint differences between GWPS and the target run against the same-cell hypothesis. This
  makes a positive result conservative and a negative result ambiguous.

**Next-best datasets for the missing levels:**

* **Same cell, different lab:** VIPerturb-seq K562 (GREEN; needs an R/Seurat conversion; shallow, so a
  depth-matched design).
* **Technical level:** a public K562 dataset re-sequenced on two platforms, if one exists. None was identified
  locally, so this needs a search.
