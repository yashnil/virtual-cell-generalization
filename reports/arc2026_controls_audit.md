# Audit of the official Arc Virtual Cell Challenge 2026 validation controls

Date: 2026-09-18
Status: complete, all invariants passed
Scope: **read-only**. No raw file was modified, no model was trained, and the
A/B/C contexts were not used as perturbation-response training data.

Reproduce with:

```bash
uv run python scripts/audit_arc2026_controls.py
```

Code: `src/virtual_cell/data/arc2026.py`.
Tests: `tests/test_arc2026_controls.py` (29 tests).
Generated tables and figures: `outputs/arc2026_controls_audit/` (git-ignored).

---

## 1. Source and provenance

| item | value |
|---|---|
| Source | https://virtualcellchallenge.org/app/datasets |
| Downloaded | 2026-09-18 |
| Working copy | `data/raw/arc2026/controls/` |
| Provenance note | `data/provenance/arc2026_controls.md` |
| Recorded checksums | `data/provenance/arc2026_controls_sha256.txt` |
| Manifest copy | `data/provenance/arc2026_manifest.json` |

All six files match their recorded SHA-256 digests before and after the audit
(`shasum -a 256 -c data/provenance/arc2026_controls_sha256.txt`, verified in
`test_official_files_match_their_recorded_checksums`). The working copy of
`manifest.json` is byte-identical in content to the provenance copy.

The raw bundle is excluded from git by `.gitignore` (`data/raw/*`, `*.h5ad`),
and so is `outputs/`. The figures referenced below are therefore reproducible
artefacts, not committed files.

---

## 2. Manifest

`data/raw/arc2026/controls/manifest.json`:

```json
{
  "season": "2026",
  "partition": "val",
  "panel_id": "vcc2026-val-1",
  "contexts": ["A", "B", "C"],
  "pert_col": "target_gene",
  "context_col": "context",
  "control_label": "non-targeting",
  "n_genes": 18533,
  "n_constructs": 300,
  "per_context": {
    "A": {"n_perturbations": 300, "control_cells": 18400,
          "ground_truth_cells": 138400, "n_ntc_ids": 46},
    "B": {"n_perturbations": 300, "control_cells": 18400,
          "ground_truth_cells": 138400, "n_ntc_ids": 46},
    "C": {"n_perturbations": 300, "control_cells": 18400,
          "ground_truth_cells": 138400, "n_ntc_ids": 46}
  },
  "cells_per_pert": 400
}
```

Field meanings, each confirmed against the data rather than assumed:

| field | meaning | verified against data |
|---|---|---|
| `season` / `partition` | 2026 season, validation split (not the final D/E/F split) | n/a |
| `panel_id` | identifier of this gene/perturbation panel | n/a |
| `contexts` | the three opaque validation contexts | matches the three `context_*.h5ad` labels exactly |
| `pert_col` | `obs` column naming the perturbation | present in all three files |
| `context_col` | `obs` column naming the context | present in all three files |
| `control_label` | value of `pert_col` marking unperturbed cells | the only value present in all three files |
| `n_genes` | genes in the expression panel | equals `n_vars` in all three files and the row count of `gene_names.csv` |
| `n_constructs` | targeted perturbations to predict | equals the row count of `pert_counts.csv` |
| `cells_per_pert` | cells per perturbation in a submission | equals the cells behind every `ntc_id` (400) |
| `control_cells` | control cells supplied per context | equals `n_obs` (18,400) in all three files |
| `n_ntc_ids` | distinct non-targeting guides | equals the 46 `ntc_id` categories in all three files |
| `ground_truth_cells` | cells in the hidden ground truth | `= control_cells + n_constructs x cells_per_pert = 18,400 + 300 x 400 = 138,400`; exact |

The manifest is internally self-consistent: `n_ntc_ids x cells_per_pert = 46 x 400
= 18,400 = control_cells`, and the `ground_truth_cells` identity holds exactly.
The latter shows that the hidden ground truth **includes** the control cells,
which the six-metric scoring will exclude from the submission itself.

## 3. Side-car tables

### `gene_names.csv`

| property | value |
|---|---|
| rows | 18,533 |
| columns | `gene_name` (single column) |
| dtype | string / `object` |
| duplicates | 0 |
| missing values | 0 |
| first entries | `TSPAN6, TNMD, DPM1, SCYL3, C1orf112` |
| last entries | `AC007244.1, GABARAPL2, SLC39A4, TKT, PRH1` |

**Ordering matters.** This file defines the column order of a submission, and
it matches the `var_names` order of all three `.h5ad` files exactly,
element-for-element. The file is *not* alphabetically sorted; it must be used
as given.

### `pert_counts.csv`

| property | value |
|---|---|
| rows | 300 |
| columns | `target_gene` (single column) |
| dtype | string / `object` |
| duplicate rows | 0 |
| distinct target genes | 300 |
| contains `non-targeting` | no |

All 300 target genes are present in the 18,533-gene expression panel, so every
perturbation is measurable in the space it is scored on. Despite the file name,
this release carries **no count column** — it is a list of target genes only, so
the per-perturbation cell budget comes from `cells_per_pert` (400) in the
manifest, uniformly across all 300 targets.

Ordering of `pert_counts.csv` is not itself a submission constraint, but the
row set is: 300 targets x 400 cells x 3 contexts = 360,000 predicted cells per
submission phase.

---

## 4. Per-context statistics

All three files: `AnnData`, `X` stored as **CSR sparse**, dtype **float32**,
`obs` columns `['target_gene', 'context', 'ntc_id']`, `var` columns `[]` (index
only), no `layers`, `obsm`, `varm`, `obsp`, `varp` or `uns` entries.

| | A | B | C |
|---|---|---|---|
| shape (cells x genes) | 18,400 x 18,533 | 18,400 x 18,533 | 18,400 x 18,533 |
| X sparse / type | yes / `csr_matrix` | yes / `csr_matrix` | yes / `csr_matrix` |
| X dtype | float32 | float32 | float32 |
| stored entries | 109,910,478 | 101,546,777 | 108,213,983 |
| explicit stored zeros | 0 | 0 | 0 |
| sparsity | 0.6777 | 0.7022 | 0.6827 |
| min count | 0 | 0 | 0 |
| max count | 977 | 1,753 | 2,462 |
| all finite | yes | yes | yes |
| all non-negative | yes | yes | yes |
| all integer-valued | yes | yes | yes |
| obs_names unique | yes (0 dupes) | yes (0 dupes) | yes (0 dupes) |
| var_names unique | yes (0 dupes) | yes (0 dupes) | yes (0 dupes) |
| duplicate cells | none | none | none |
| duplicate genes | none | none | none |

Counts are integer-valued but stored as float32; anything that assumes an
integer dtype must cast explicitly.

### Library size (total UMI counts per cell)

| | A | B | C |
|---|---|---|---|
| min | 3,275 | **710** | 3,447 |
| median | 20,109 | 19,946 | 20,034 |
| mean | 21,133.5 | 19,989.9 | 21,156.5 |
| max | 52,420 | 42,666 | 50,965 |
| CV | 0.459 | 0.385 | 0.438 |
| 0.5th percentile | 3,930 | **858** | 4,119 |
| cells < 2,000 counts | 0 | **458 (2.49%)** | 0 |
| cells < 5,000 counts | 323 (1.76%) | 874 (4.75%) | 247 (1.34%) |

### Genes detected per cell

| | A | B | C |
|---|---|---|---|
| min | 2,205 | **580** | 1,937 |
| median | 6,147 | 5,756 | 6,006 |
| mean | 5,973.4 | 5,518.8 | 5,881.2 |
| max | 8,453 | 7,649 | 8,482 |

### Metadata

Identical in structure across the three contexts:

- `target_gene`: categorical, one category, `'non-targeting'`, 18,400 cells.
  **These are control cells only; no perturbed cells are present.**
- `context`: categorical, one category, exactly `'A'` / `'B'` / `'C'`.
- `ntc_id`: categorical, 46 categories, all of the form `non-targeting-<n>`,
  each with exactly **400 cells**. The same 46 guide identifiers appear in all
  three contexts, so the non-targeting guide panel is shared.
- `obs_names`: `<context>_ctrl_<6-digit index>`, e.g. `A_ctrl_000000`.
- `var_names`: HGNC-style gene symbols, no `var` columns.

---

## 5. Cross-context invariant checks

44 invariants, **all passed** (`outputs/arc2026_controls_audit/invariant_checks.csv`):

| invariant | result |
|---|---|
| context labels match the manifest (`['A','B','C']`) | PASS |
| each file stores its label exactly (no case/whitespace drift) | PASS |
| no accidental context-label swaps (one label per file, matching the filename) | PASS |
| gene sets identical across A/B/C | PASS (0 genes differ) |
| gene **order** identical across A/B/C | PASS |
| `var_names` match `gene_names.csv` exactly in order | PASS (all 18,533) |
| cell counts match the manifest | PASS (18,400 each) |
| gene counts match the manifest | PASS (18,533 each) |
| only control cells present (`target_gene == 'non-targeting'`) | PASS |
| guide count matches the manifest (46) | PASS |
| every guide has exactly `cells_per_pert` (400) cells | PASS |
| manifest cell arithmetic self-consistent | PASS |
| non-targeting guide panel identical across contexts | PASS |
| `obs_names` unique within each context | PASS |
| cell ids disjoint across contexts | PASS (0 shared, all pairs) |
| `var_names` unique within each context | PASS |
| counts finite, non-negative, integer-valued | PASS (all three) |
| raw count structure consistent across contexts | PASS (CSR/float32/no explicit zeros everywhere) |

Cell identifiers are namespaced by context (`A_ctrl_*`, `B_ctrl_*`,
`C_ctrl_*`), so the zero overlap is by construction rather than coincidence.

### One failure raised and resolved during the audit

The first run reported `var_names match gene_names.csv exactly in order` as
FAILED for all three contexts. This was **not** a data discrepancy: the `.h5ad`
`var` index carries the pandas `string` extension dtype while the CSV column
loads as `object`, and `pandas.Index.equals` returns `False` across those dtypes
even when every label matches. All 18,533 names were verified identical
element-for-element in identical order. The comparison now goes through
`arc2026.same_labels`, which casts both sides to `object` and compares values
and order only — no whitespace, case, or unicode normalisation, so a genuine
mismatch still fails (covered by `test_same_labels_is_order_and_value_sensitive`).

---

## 6. Basal-state exploratory analysis

Pseudobulk = per-gene mean of `log1p(1e4 * count / library_size)` over all
18,400 control cells, matching
`virtual_cell.preprocessing.pseudobulk.mean_expression`. Raw-count means are
also written, to `pseudobulk_raw_mean.csv`.

### Pairwise correlations (18,533 genes)

| pair | Pearson | Spearman |
|---|---|---|
| A vs B | **0.6882** | 0.7879 |
| A vs C | **0.6022** | 0.7579 |
| B vs C | **0.7320** | 0.8513 |

Restricted to the 14,153 genes non-zero in all three contexts, Pearson is
0.6733 / 0.5686 / 0.7146 for the same pairs — the correlations are not an
artefact of shared zeros. On `log1p` of the pseudobulk they are 0.7268 / 0.6566
/ 0.7732.

So the ordering is consistent across every variant: **B and C are the closest
pair at baseline, A and C the most distant, and A is the most distinctive
context of the three.**

### Gene coverage

| | A | B | C |
|---|---|---|---|
| genes with zero total counts | 2,398 | 2,853 | 2,317 |

861 genes are zero in all three contexts; 3,519 are zero in at least one
context but not all; 14,153 are expressed in all three.

### Top genes contributing to basal differences

Descriptive statistic only: `spread = max - min` of the per-context pseudobulk
value. No test, ranking model, or identity inference. Distribution of spread
across all 18,533 genes: median 0.113, mean 0.241, max 4.139; 2,861 genes
exceed 0.5 and 725 exceed 1.0.

| gene | A | B | C | spread | max in | min in |
|---|---|---|---|---|---|---|
| CLU | 0.032 | 4.171 | 1.269 | 4.139 | B | A |
| S100A6 | 0.000 | 4.057 | 2.788 | 4.056 | B | A |
| VIM | 2.754 | 4.068 | 0.013 | 4.056 | B | C |
| ACTB | 0.269 | 0.414 | 4.016 | 3.746 | C | A |
| BASP1 | 0.000 | 3.509 | 1.372 | 3.509 | B | A |
| ARHGDIB | 3.520 | 0.022 | 0.376 | 3.499 | A | B |
| TRBC1 | 3.467 | 0.000 | 0.000 | 3.467 | A | C |
| HIST1H1B | 3.337 | 1.145 | 0.000 | 3.336 | A | C |
| ANXA2 | 0.300 | 3.564 | 2.948 | 3.264 | B | A |
| TACSTD2 | 0.000 | 0.522 | 3.212 | 3.212 | C | A |
| CDKN2A | 0.000 | 3.173 | 0.369 | 3.173 | B | A |
| S100A11 | 0.001 | 1.725 | 3.153 | 3.152 | C | A |
| KRT8 | 0.027 | 3.157 | 2.642 | 3.129 | B | A |
| KRT7 | 0.000 | 3.064 | 0.053 | 3.063 | B | A |
| ADA | 3.091 | 0.069 | 1.329 | 3.022 | A | B |

Full 50-gene table: `outputs/arc2026_controls_audit/top_basal_difference_genes.csv`.

### Figures

All under `outputs/arc2026_controls_audit/` (git-ignored, regenerated by the script):

| figure | content |
|---|---|
| `library_size_distributions.png` | library size per context, linear and log10 |
| `genes_detected_distributions.png` | genes detected per cell per context |
| `pseudobulk_scatter.png` | A-B, A-C, B-C pseudobulk scatter with identity line and r |
| `control_cell_pca.png` | PCA of 2,000 control cells per context (2,000 most expressed genes, seed 0) |

The PCA separates the three contexts into three tight, non-overlapping clusters
on PC1 and PC2 (26.8% and 15.3% of variance), with within-context spread far
smaller than between-context distance. PC3 onwards drop to 2.5% and below.

---

## 7. Findings that are not invariant failures

None of the following blocks the reproduction gate, but each is recorded
because it constrains later work.

1. **Context B has a genuine low-depth tail that A and C do not.** 458 B cells
   (2.49%) fall below 2,000 UMIs and 178 below 1,000, against zero such cells in
   A and C, whose minima are 3,275 and 3,447. B also has the lowest median genes
   detected. Any per-cell filtering threshold will therefore remove cells from B
   and almost none from A or C, which silently changes the effective control
   sample size per context. Decide the threshold once, apply it identically, and
   report the per-context cell loss.

2. **The gene panel excludes the highest-abundance transcript families.** The
   18,533-gene panel contains **zero** cytoplasmic ribosomal protein genes
   (`RPL*`/`RPS*`), **zero** mitochondrial ribosomal genes (`MRPL*`/`MRPS*`) and
   **zero** mitochondrial rRNA (`MT-RNR*`); it keeps the 12 protein-coding `MT-`
   genes. GAPDH, EEF1A1, PTMA, MALAT1, NEAT1, LDHA and XIST are also absent.
   Consequences: the expression profile is much flatter than raw 10x data (the
   top 10 genes hold only 3.0-4.7% of the library, versus 15-30% typical for
   unfiltered 10x), the mitochondrial fraction is 0.30-0.60%, and **absolute
   expression levels are not directly comparable to unfiltered public data**.
   Any Replogle/Nadig comparison must intersect gene sets first and expect to
   lose the ribosomal block entirely. Pinned by
   `test_panel_excludes_the_highest_abundance_transcript_families` so a changed
   panel in a future download fails loudly.

3. **The counts behave like real scRNA-seq, not like a Poisson simulation.**
   This was checked because the flat profile in (2) is unusual. Per-gene
   variance/mean over cells has median 1.45-1.63 with a tail to 208-377, i.e.
   genuinely overdispersed. Gene labels also carry real covariance structure:
   the 12 `MT-` genes have mean pairwise r = 0.42 across cells against 0.01 for
   random expressed genes, and the HIST1 cluster is likewise elevated
   (mean r = 0.11, max 0.93). The 300 CRISPRi targets are expressed in the
   controls (median mean-count 0.86-1.23 versus 0.13-0.20 for a typical panel
   gene, and none are zero) — consistent with a panel designed to be knocked
   down. Gene symbols are therefore biologically meaningful and gene-keyed
   priors (DepMap, GO, network features) can legitimately be joined onto them.

4. **ACTB is near-absent in A and B (mean 0.88 and 1.28 counts per cell) but is
   the top gene in C (120 counts).** Unlike (2), this is not explained by panel
   filtering, since ACTB is in the panel. It is recorded as an unresolved
   observation. No identity inference was attempted, and none should be read
   into it.

5. **Basal correlation between contexts is low.** Pearson 0.60-0.73 is below
   what is typically seen between human cell lines on log-normalised pseudobulk
   and below the 0.82-0.85 of this repository's own synthetic contexts. Combined
   with the clean PCA separation, the three validation contexts are far apart at
   baseline. This raises, not lowers, the difficulty of the zero-shot transfer
   the project is built around, and it means "nearest context" style baselines
   have little to lean on.

6. **`pert_counts.csv` carries no count column** in this release, only
   `target_gene`. The per-perturbation cell budget must be taken from
   `cells_per_pert = 400` in the manifest. This resolves, for the validation
   phase, the cells-per-perturbation ambiguity recorded in the research log:
   the official bundle is uniform at 400.

---

## 8. Interpretation

Strictly descriptive, limited to baseline state:

- The three validation contexts are **internally clean and mutually
  consistent**: same genes in the same order, same guide panel, same cell
  budget, disjoint cell ids, raw integer counts throughout.
- They are **well separated at baseline**. Pseudobulk correlations of 0.60-0.73
  and complete PCA separation say the contexts differ substantially in basal
  transcriptional state, with B and C closer to each other than either is to A.
- Sequencing depth is **comparable in the median** (about 20,000 UMIs and
  5,500-6,100 genes detected per cell in all three) but **not in the lower
  tail**, where B is distinctly shallower.

No attempt was made, and none should be made, to identify the biological
identity of A, B or C. The purpose of this section is only to characterise how
different the unseen contexts are at baseline.

---

## 9. Standing constraint

**A/B/C are evaluation inputs only.** They contain unperturbed control cells
exclusively; their perturbation responses are hidden. They may be used to
derive context representations (basal state, depth, gene coverage) and to
format submissions. They must never be used as perturbation-response training
data, and no held-out context's responses may enter training or model
selection.

**No novel model work follows from this audit.** The Molina & Zhang
reproduction gate remains mandatory before any interaction model or other novel
architecture is built.
