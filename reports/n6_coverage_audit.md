# N6 coverage audit: VIPerturb-seq K562 as an independent-lab K562 source

Date: 2026-10-07. Scope: identifiers, QC counts and provenance only. **No source→target compatibility, and no
VIPerturb response vector, was computed.**

* Conversion: `reports/viperturb_conversion_audit.md`.
* Frozen axes: `data/splits/n6_k562/` (`sha256.txt`).

## 0. N5 frozen first

* The N5 protocol digest and the output manifest re-verified with 0 failures.
* `reports/n5_results.md` carries the mechanical Outcome-A verdict and the within-study interpretation.
* All N5 artefacts are hashed in `data/provenance/research_v3/n5_freeze_sha256.txt`. `outputs/n5` is now
  read-only.
* N6 writes only to `outputs/n6`, `data/splits/n6_k562` and new scripts.

## 1. Dataset (details in the conversion audit)

| item | VIPerturb-seq K562 | Replogle K562 essential (target) |
|---|---|---|
| lab | Satija (NYGC/NYU) | Weissman (Whitehead) |
| effector | dCas9-KRAB-**MeCP2** | dCas9-KRAB (per publication) |
| library | GuEST-List, **3 single sgRNAs per gene** | essential **dual-sgRNA** pairs |
| timepoint | day 6 | ~day 6 (per publication) |
| chemistry | **10x Flex v2: probe-based, fixed cells** | 10x 3′ (poly-A capture, live cells) |
| sequencing | **Ultima UG100** | not recorded locally |
| cells / controls | 906,837 / 7,949 `NO-TARGET` | 308,646 / 10,691 |
| genes measured | 19,068 probe targets | 8,563 |
| raw integer counts | yes (converted, lossless) | no (log1p CP10K) |
| licence | CC BY 4.0 | CC BY 4.0 |

## 2. Perturbation overlap

Counts use cells passing the scPertEval cell filter (all VIPerturb cells pass) and the frozen alias map:
GET1 → WRB, GET3 → ASNA1, MICOS10 → MINOS1.

| reference set | size | in VIPerturb (any cells) | ≥ 30 VIPerturb cells | median cells among ≥ 30 |
|---|---|---|---|---|
| **N5 panel (frozen; target ∩ GWPS ∩ six contexts)** | 1,054 | 1,040 (98.7 %) | **637 (60.4 %)** | **47** (range 30–119) |
| N3 panel | 1,062 | 1,045 | 637 | 47 |
| Replogle K562 essential (≥ 30 cells) | 1,971 | 1,937 | 1,039 | — |
| GWPS K562 (≥ 30 cells) | 9,675 | 9,447 | 7,293 | — |
| RPE1 essential | 2,016 | 1,964 | 1,162 | — |
| HepG2 | 1,818 | 1,778 | 1,017 | — |
| Jurkat | 2,137 | 2,084 | 1,233 | — |

The non-panel rows come from the manifest, which matches the converted object exactly, and do not apply the
alias map. They are indicative only.

**Naming and mapping:**

* VIPerturb uses some **older HGNC symbols**. Three panel genes are pure renames (GET1/WRB, GET3/ASNA1,
  MICOS10/MINOS1). The map is frozen with exactly these three entries and applied to perturbation labels. The
  measured-gene (probe) axis uses current symbols; the alias is a fallback only.
  * GET1 (47 cells) and GET3 (31) are included.
  * MICOS10 (25 cells) falls below the floor.
* 14 panel genes are absent from VIPerturb: ADAT3, CHEK1, HAUS7, PRELID3B, PRPF19, PTCD1, RPA1, RPL11,
  SIGLEC14, SPC25, TBC1D3, TMEM199, ZNF284, ZNF763. These are mostly strongly essential genes and are not
  resolvable by renaming.
* **Guide structure.** Of the 637 included perturbations, 583 have 3 guides, 53 have 2 and 1 has 1. All of a
  gene's cells are pooled, matching the suffix-stripped labels of the target.
* **Duplicate representations.** None within a gene label. The only duplicated barcodes are the controls
  (resolved in conversion).
* **Controls.** 7,949 `NO-TARGET` cells, pooled across wells (as in scPertEval and N5).

**Selection caveat (source-only, but not neutral).**

* 403 panel perturbations have 1–29 VIPerturb cells. Low counts at day 6 plausibly reflect fitness depletion
  of strongly essential knockdowns.
* The N6 panel is therefore **enriched for milder perturbations** relative to N5. Every source is evaluated on
  the same 637, so contrasts are fair, but N6 values are not directly comparable to the N5 values (the panel
  and centring differ).

## 3. Gene overlap

* N5 genes (6,408) present on the VIPerturb probe axis: **6,083**. This is the N6 gene axis. The probe axis uses current
  symbols, so no alias is needed for genes; only perturbation labels use older symbols.
* The N5 response definition is reproducible on VIPerturb raw counts: cell filter ≥ 200 genes, CP10K over the
  object's full 19,068-gene axis, log1p, restriction to the axis.
* The **scientific object differs in one unavoidable way.** Flex measures transcripts through **probes**, not
  poly-A capture, so per-gene detection efficiency and the normalisation denominator differ from 3′ data.
  * Template removal, which centres each gene over perturbations, absorbs gene-wise detection offsets.
  * It does not absorb gene-wise scaling differences. Those can lower measured compatibility. This is part of
    the "independent platform" bundle under test, not something to correct away.

## 4. Reliability feasibility

* At 30–119 cells per perturbation (median 47), the N1 disjoint F/E1/E2 split gives about 23 / 11 / 11 cells at
  the median, and at least 15 / 7 / 7. That matches the minimum used throughout N1–N5.
* 7,949 controls give a three-way control split of about 3,974 / 1,987 / 1,987.
* Per-perturbation reliability will be lower than GWPS (median 223 cells), similar to the depth-matched GWPS
  of N5 (median 89 cells; endpoint stable).
* Whether pooled reliability is adequate is decided by the **preregistered gate** (protocol §4), not here.

## 5. Feasibility verdict

**PARTIAL (strong): usable, but materially restricted.**

**Supports.** A same-cell-line source from an **independent lab, library, effector, chemistry and sequencer**,
tested against the frozen comparators on **637** shared perturbations × **6,083** genes. That is enough for
paired-bootstrap uncertainty of the kind obtained in N5.

**Restrictions:**

* the 60 % panel, enriched for milder perturbations;
* shallow depth (median 47 cells);
* reliability that is unknown until the gate runs;
* all technical axes changing **together**. N6 can show whether same-cell compatibility survives an entirely
  independent screen, but it cannot attribute a loss to any single technical factor.

**Largest technical difference from Replogle.** The **probe-based, fixed-cell 10x Flex chemistry, read on
Ultima**, versus live-cell 3′ capture. It changes which molecules are counted and their per-gene efficiency,
on top of a different effector (KRAB-MeCP2) and single- vs dual-guide design.

Note the relevance: the Arc 2026 data are 10x Flex on Ultima.

**Identifiability.** Sufficient to run N6 as preregistered. A strong VIPerturb result identifies that same-cell
compatibility survives an independent lab/platform. A weak result is attributable to "screen/platform
compatibility" only if the depth-matched positive control passes. Otherwise it is Outcome D.

**Alternatives** if N6 proved uninformative: another independent K562 screen with raw counts, for example from
the Arc 2026 or other Flex-based releases. None is identified locally.
