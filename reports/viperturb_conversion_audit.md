# VIPerturb-seq conversion audit (Seurat → h5ad)

Date: 2026-10-07. No compatibility or transfer quantity was computed. Only file integrity, structure,
identifiers and QC counts were examined.

## 1. Source

| item | value |
|---|---|
| publication | Bradu A., Blair J.D., Grabski I.N., Mascio I., Lee J., McCormick C., Satija R. "Genome-wide single-cell perturbation screens with VIPerturb-seq." bioRxiv 2026, doi:10.64898/2026.02.12.705613 (Satija lab, NYGC/NYU) |
| dataset | Zenodo record **10.5281/zenodo.18460279** (concept 10.5281/zenodo.18460278), published 2026-02-02; the latest version (checked 2026-10-07) |
| licence | **CC BY 4.0** (Zenodo `metadata.license.id`). Attribution required |
| files used | `genome_wide_binA.RDS` (3,558,351,019 B, MD5 `31d69084…`), `genome_wide_binB.RDS` (3,781,446,726 B, MD5 `0b244486…`), `genome_wide_binC.RDS` (2,885,854,190 B, MD5 `af9e5908…`), `genome_wide_manifest.txt` (398,905 B, MD5 `c095c8f2…`) |
| files not used | `genome_wide_filtered.rds` (pre-filtered on the authors' perturbation "fingerprints", a phenotype-based selection); `vimentin_screen.rds`; `multimodal_cell_line_mixing_pilot.rds` |
| download | `data/raw/viperturb/download.py`. Every MD5 equals Zenodo's (`data/raw/viperturb/download_log.jsonl`) |

**Biology and technology** (from the preprint; not encoded in the objects):

* K562 expressing **dCas9-KRAB-MeCP2** (CRISPRi).
* **GuEST-List** library: 57,050 sgRNAs, **3 single sgRNAs per gene**, adapted from Dolcetto set A.
* MOI 0.3–0.5; harvested **day 6** after puromycin selection.
* **10x Flex v2 (Apex)**: probe-based, fixed cells, with combinatorial indexing over 24 wells.
* Split-probe guide barcodes, demultiplexed with MULTISeqDemux.
* Sequenced on **Ultima UG100**.
* About 880,000 cell barcodes; a median of 42 cells per targeted gene reported.

## 2. Environment (isolated; nothing installed into base)

* Dedicated mamba environment `vcc-viperturb-convert` (miniforge): **R 4.4.3, SeuratObject 5.4.0, Matrix
  1.7-6, jsonlite 2.0.0, digest 0.6.39**.
* Full spec: `scripts/research_v3/viperturb/environment_vcc-viperturb-convert.yml`.
* Seurat itself is not needed: objects are read with SeuratObject.

## 3. Procedure

1. **`scripts/research_v3/viperturb/export_seurat.R`.** Per bin:
   * `readRDS`;
   * record class, assays, layers and metadata columns;
   * take `LayerData(obj, "RNA", "counts")` (dgCMatrix); assert integer and non-negative;
   * write the exact CSC slots `i`, `p` and `x` as little-endian int32, plus gene names, barcodes and the full
     `meta.data` CSV;
   * write validation statistics computed **in R on the source object**: dimensions, nnz, total counts, the
     column sums of 2,000 random cells, the row sums of 500 random genes, and 5,000 random entries
     (seed 20261006).
2. **`scripts/research_v3/viperturb/assemble_h5ad.py`.** Rebuilds each bin in Python and re-checks every R
   statistic. Then it:
   * requires an identical gene axis across bins;
   * requires every barcode present in more than one bin to have **identical** counts, and keeps it once;
   * writes `data/raw/viperturb/viperturb_k562_genome_wide.h5ad` (CSR int32 counts, gzip).
3. Chain: `data/raw/viperturb/convert_rest.sh`. Log: `convert_rest.log`.

## 4. Structure of the source objects (identical across bins)

* Class: `Seurat`. Assays: `RNA` (layer `counts`) and `GDO` (guide counts, 57,050 features; layers `counts`,
  `data`). Default assay RNA.
* `meta.data`: `orig.ident, nCount_RNA, nFeature_RNA, nCount_GDO, guide, gene, sample, nFeature_GDO`.
  * `gene` is the target symbol (`NO-TARGET` for controls). `guide` is e.g. `PRMT7g0`. `sample` is the
    well/lane.
* Peak R memory: 54 / 57 / 44 GB (bins A / B / C).

## 5. Validation results (`data/provenance/research_v3/viperturb_conversion_validation.json`)

| check | bin A | bin B | bin C |
|---|---|---|---|
| dimensions (genes × cells) | 19,068 × 321,022 | 19,068 × 341,182 | 19,068 × 260,531 |
| nnz | 1,676,324,776 | 1,780,612,863 | 1,357,917,183 |
| total counts | 5,118,746,603 | 5,436,366,789 | 4,137,263,374 |
| dims, nnz, totals, 2,000 column sums, 500 row sums, 5,000 entries, barcodes = metadata rows, integer ≥ 0 | all pass | all pass | all pass |
| gene axis identical to bin A | — | pass | pass |

**Merged object:**

* 922,735 bin cells − 15,898 duplicate barcodes = **906,837 cells** × **19,068 genes**.
* nnz 4,731,586,400, which equals the bin sum minus the two extra copies of the control nnz (exact).
* Total counts 14,436,147,174.
* **Duplicated barcodes:** 15,898, all **byte-identical** across bins. These are the 7,949 `NO-TARGET`
  controls, present in all three bins by design.
* **Per-label cell counts:** they match the Zenodo manifest **exactly for all 18,886 labels**.
* **Guides per targeted gene:** 3 for 16,598, 2 for 2,115, 1 for 170 (and 2 labels with 6).
* Every cell has ≥ 200 detected genes.

**Artefact.**

* `data/raw/viperturb/viperturb_k562_genome_wide.h5ad`, SHA-256
  **`04ef7c27f3be73cde138d4b91871c7527671a7f2632d48a187d7f0b66c7d4372`**.
* Recorded in the validation JSON and the N6 digest.

**Verdict: conversion lossless for raw counts, barcodes, perturbation and control identity, gene identifiers
and metadata.**

## 6. Caveats

* The objects are the authors' processed release: cell calling, guide assignment and their QC. Upstream
  FASTQ-level choices are theirs.
* The `GDO` guide-count assay was not exported. Perturbation identity comes from the authors' `guide`/`gene`
  assignment.
* Gene symbols follow the authors' annotation. Some are older HGNC symbols (e.g. `WRB` for GET1), handled by
  the frozen alias map in N6.
