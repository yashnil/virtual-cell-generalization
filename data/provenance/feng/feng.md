# Feng et al. targeted CRISPRi screen — per-line log fold changes

Acquired 2026-09-21, approved for this file only. **The count matrix was not
downloaded.**

| | |
|---|---|
| publication | Feng et al., *Cell Genomics* 2026 — "A genome-scale single-cell CRISPRi map of trans gene regulation across human pluripotent stem cell lines" |
| identifiers | PMC12903452 · PubMed 41330380 · bioRxiv 2024.11.28.625833 |
| raw sequencing | ENA **ERP165335** |
| processed data | Figshare DOI `10.6084/m9.figshare.26819743` |
| count data | Figshare DOI `10.6084/m9.figshare.27989294` (**not downloaded**) |
| licence | MIT |
| retrieved | 2026-09-21 |

## File

| | |
|---|---|
| name | `TargetedScreen_LFC_byGene-perLine.tsv.gz` |
| figshare file id | 49291675 |
| bytes | **1,919,661,029** (matches the API's declared size) |
| upstream md5 | `27ae8d0109b9c42f1f972c45987db77e` (**verified after download**) |
| local sha256 | see `feng_sha256.txt` |
| transfer | resumable HTTPS (`curl -C -`) |

## Schema

Long format, one row per (cell line, target, expressed gene):
`Target`, `Cell_Line`, `Expressed_Gene_Symbol`, `Expressed_Gene_Ens_ID`,
`lfc`, `pval_adj` — confirmed locally, see
`outputs/feng_multicontext_v1/schema.json`.

**LFC definition.** A log fold change of the targeted condition against the
screen's non-targeting control guides, computed per cell line by the authors'
own model. This is *not* the same estimator as our source deltas, which are
differences of log-normalised pseudobulk means. Both are log-space contrasts
against a non-targeting control, so they are comparable in kind, but a
correlation computed across the two inherits the estimator difference. This
limitation is stated in the report and is not removable without the count
matrix.

## Design, from the companion cell metadata

`TargetedScreen_Cell-Metadata.tsv.gz` (13.5 MB, fetched as an incidental audit
file) gives 1,161,864 cells with `Cell_ID`, `Batch`, `Guide_Call`, `Cell_Line`.

- **19 cell lines** (the publication describes 20 for the targeted screen; the
  released table carries 19), 47,999–76,059 cells each.
- 446 distinct guide-call targets including `NonTarget`; 444 gene targets.
- **Median 74 cells per (line, target)**; q10 = 22, q90 = 123. Only 25.1% of
  (line, target) pairs reach 100 cells; 73.6% reach 50.
- **Non-targeting control cells are scarce: 8,241 in total (0.7%)**, median 488
  per line, **minimum 86**. Every LFC in a line shares that line's control
  estimate, so a line with few control cells carries a correlated noise floor
  across all of its targets. This is recorded here because it bears directly on
  the reliability qualification.

## Leakage assessment

**No leakage risk for the Arc task.** These are perturbation outcomes, but in
iPSC lines unrelated to Arc's A/B/C, and they are used only to *evaluate* a
frozen predictor. Nothing in this phase is fitted to them.
