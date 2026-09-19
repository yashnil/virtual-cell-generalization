# Research log

Running record of what was done, what was decided, and what is still unknown.
Newest entries at the bottom. Dates are absolute (YYYY-MM-DD).

Conventions:

- **Know** = supported by a cited source or by a test/experiment in this repo.
- **Suspect** = hypothesis, not yet tested.
- **Don't understand** = open question that needs reading or an experiment.
- **Disproved** = something an experiment or source showed to be wrong.

---

## 2026-09-06 — Phase 1: project understanding and data infrastructure

### Context

The repository contained `plans.MD` (the full project specification), a `uv`
project with Python 3.11, empty package stubs under `src/virtual_cell/`, a
one-file synthetic-data script, three generated synthetic `.h5ad` files, and
empty report files. No commits existed yet.

### What was implemented

Data infrastructure only. No model training, no real data download.

- `src/virtual_cell/data/io.py`
  - `load_context(path, context_key="context", expected_context=None,
    require_integer=True)`: loads one `.h5ad` and runs integrity checks.
  - `check_integrity` / `check_counts` / `check_context_labels`: reusable
    checks (non-empty; unique gene names; unique cell ids; finite,
    non-negative, integer-valued `X`; exactly one context label per file;
    exact-string label match; no silent case or whitespace normalisation).
  - `load_contexts(paths)`: loads several files keyed by their context label,
    verifying labels when a mapping is given and rejecting duplicate labels.
  - `DataIntegrityError` for all assumption violations.
- `src/virtual_cell/data/summary.py`
  - `library_sizes`, `sparsity`, `genes_detected_per_cell`.
  - `ContextSummary` dataclass and `summarize_context` / `summarize_contexts`
    (cells, genes, mean/median/min/max library size, sparsity, mean genes
    detected). Sparse and dense inputs give identical results; explicitly
    stored zeros in sparse matrices are not counted as detected.
- `src/virtual_cell/data/synthetic.py`
  - `make_synthetic_context`, `make_synthetic_contexts`,
    `write_synthetic_contexts`. Deterministic given a seed. The default now
    shares one per-gene rate vector across contexts with log-normal
    per-context modulation, so basal profiles are correlated (r about 0.82 to
    0.85 between synthetic contexts) as real cell lines are. The original
    independent-rates behaviour is available with `shared_base_rate=False`.
- `src/virtual_cell/preprocessing/pseudobulk.py`
  - `normalized_matrix` (scale to target sum, optional log1p, zero-count cells
    stay zero), `mean_expression`, `shared_genes`, `basal_mean_table`,
    `basal_mean_correlation`, `basal_log_fold_change`.
- `scripts/make_synthetic_controls.py`: thin CLI over the synthetic module;
  path is resolved relative to the repo root so it works from any cwd.
- `scripts/explore_synthetic_contexts.py`: prints AnnData structure, count
  matrix corner, per-context summary table, basal-mean table and correlation
  matrix, log2 fold changes, and saves `outputs/exploration/*.csv` and
  `synthetic_context_overview.png`.
- `tests/` (31 tests): loading synthetic contexts from disk and memory,
  non-negative integer counts, exact context-label preservation (including
  case), duplicate genes, duplicate cells, negative / NaN / non-integer
  values, missing or multiple labels, empty matrices, missing files, duplicate
  labels across files, summary statistics against hand-computed values,
  sparse vs dense agreement, explicit-zero handling, normalisation row sums,
  gene alignment across contexts, synthetic determinism and correlation.
- `pyproject.toml`: build backend now packages both `virtual_cell` and
  `virtual_cell_generalization` (previously only the hello-world module would
  have been packaged; the editable install hid this). Added pytest and ruff
  configuration.
- `.gitignore`: added `outputs/` and `reports/figures/`; fixed a missing
  trailing newline that had glued `.DS_Store` to the next line.
- `README.md` rewritten; `reports/literature_notes.md` written from web
  research (see that file for sources).

### Commands run

```
uv sync
uv run pytest            # 31 passed
uv run ruff check .      # All checks passed
uv run ruff format .
uv run python scripts/make_synthetic_controls.py
uv run python scripts/explore_synthetic_contexts.py
uv build                 # wheel contains both packages
```

### Synthetic exploration result (no biological meaning)

| context | cells | genes | mean lib | median lib | sparsity |
|---|---|---|---|---|---|
| A | 200 | 100 | 158.7 | 158 | 0.306 |
| B | 200 | 100 | 186.0 | 185.5 | 0.267 |
| C | 200 | 100 | 221.0 | 221 | 0.216 |

Pairwise Pearson correlation of normalised basal means: 0.82 to 0.85.
Before the shared-baseline change it was 0.12 or lower, which would have made
"nearest context" and basal-state conditioning meaningless on synthetic data.

### Discrepancies and risks found

- **Synthetic label collision with Arc.** Arc's validation contexts are also
  named `A`, `B`, `C` (and final contexts `D`, `E`, `F`). Synthetic files must
  never be confused with the Arc control bundle. Keep them in
  `data/raw/synthetic/` and consider renaming synthetic labels (for example
  `SYN_A`) before the Arc bundle is downloaded.
- **Cells-per-perturbation ambiguity in Arc docs.** The cell-eval2 metrics
  brief (2026-08-19) says the number of cells per perturbation is
  unconstrained, but the `vcc` CLI wiki (v0.2.0, 2026-09-01) enforces exactly
  400. Treat the CLI as authoritative and re-check before submitting.
- **Arc rules page could not be fetched** (client-side rendered). Team size,
  publication obligations, and the daily submission limit remain unverified.
- **Synthetic data are Poisson**, so they have no overdispersion, no
  dropout beyond Poisson zeros, and no gene-gene correlation. They validate
  code paths only. Any summary statistic on them (sparsity about 0.2 to 0.3,
  library size about 200) is far from real 10x data (sparsity above 0.9,
  library sizes in the thousands).
- **The package `virtual_cell_generalization`** (hello-world entry point) is
  redundant with `virtual_cell`. Left in place; can be removed later.
- **Correction (same day): the 2026 decomposition paper exists.** Molina and
  Zhang, bioRxiv 10.64898/2026.07.24.740459 (posted 2026-07-27). The first
  literature pass missed it and this log wrongly said it did not exist. Its
  decomposition is exactly plans.MD section 9. COMPASS (2026-08-06) is a
  separate, related paper. See the revision entry below.
- **State reports two cross-context tasks**: an underrepresented-context
  task (30% of the held-out context's perturbations in training, where the
  headline gains are) and a separate zero-shot context task (ST+SE). Only the
  zero-shot task is comparable to a strict leave-one-context-out benchmark.
- **PrePR-CT** uses chemical/cytokine perturbations and R^2 metrics that the
  evaluation literature criticises; it still establishes that
  control-derived context graph priors are not novel.
- The "27 methods on 29 datasets" benchmark is scPerturBench (Nature
  Methods, Feb 2026), confirmed. "Virtual Cells Need Context, Not Just
  Scale" (Dibaeinia et al., bioRxiv 2026-02-09) is confirmed but is a
  position paper with single-dataset evidence.
- Arc's scorer uses cell-level Wilcoxon DE, which the pseudobulk-DE
  literature considers pseudoreplicated. We need Arc's definition for the
  challenge and a pseudobulk DE for the research track.

### Know / Suspect / Don't understand / Disproved

**Know**

- Arc 2026 task format, dates, metrics, and normalisation (primary sources in
  the literature notes).
- The loader and summary code behave as specified on synthetic and
  hand-built edge cases (tests).

**Suspect**

- Basal mean expression (pseudobulk) will be a weak but non-trivial context
  representation; PCA or pathway scores will be needed for the interaction
  model.
- The conserved-effect baseline (per-perturbation mean response across
  source contexts) will be hard to beat on Arc's metrics because most
  transferable signal is in beta. Note it is not Arc's zero anchor: Arc's
  baseline is a single context-level mean response, so a beta-only model can
  score above 0.

**Don't understand**

- How large the interaction component is in the Replogle/Nadig four-line data
  once measurement noise is removed, and how the 2026 decomposition paper
  estimated it.
- Whether Arc's 18,533-gene panel intersects cleanly with Replogle/Nadig gene
  sets and how much is lost at the intersection.
- What the Arc contexts are (undisclosed cell lines from different
  tissues); whether they resemble any public cell line enough for
  nearest-context transfer.
- How the split-half reproducibility ceiling behaves for weak perturbations.

**Disproved**

- Nothing yet. No experiments have been run on real data.

### Next steps (ordered)

1. Install `vcc-cli`, log in, download the 2026 controls bundle, and write a
   loader for the official `context_{A,B,C}.h5ad` + `gene_names.csv` +
   `pert_counts.csv` layout with the same integrity checks.
2. Download Replogle 2022 (K562, RPE1) and Nadig 2025 (HepG2, Jurkat) CRISPRi
   data; build gene intersection, per-context pseudobulk per perturbation,
   and the leave-one-context-out split files under `data/splits/`.
3. Implement differential expression (Wilcoxon on CPM, BH, matching the Arc
   definition) and the six Arc metrics locally via `cell-eval2`.
4. Implement Baselines 0 to 5 from `plans.MD` and produce the first
   leave-one-context-out benchmark table.
5. Reproduce the response decomposition on the four public cell lines and
   quantify variance explained by each component (Phase 3 gate).

---

## 2026-09-06 (revision) — Framing corrections after reading Molina and Zhang

### What changed

- Read Molina and Zhang (2026) in full (PDF and full text). Recorded in
  `reports/literature_notes.md` section 5. Corrected the false "not found"
  statement in the notes, this log, and memory.
- Corrected the Arc normalisation interpretation: 0 is the official
  mean-response baseline (a single context-level response), 1 an experimental
  replicate anchor. A perturbation-specific beta-only predictor is a different
  model and can score above 0 with gamma = 0.
- Distinguished State's underrepresented-context (30%) task from its
  zero-shot context task.
- Arc contexts are described as six cell lines from different tissues of
  origin until official metadata says otherwise.
- Appended a framing addendum to `plans.MD` and updated `README.md`.

### The negative result we now build around

Molina and Zhang: across K562/RPE1/HepG2/Jurkat, reproducible variance splits
into template 27.8%, conserved beta 29.4%, interaction gamma 23.5%, noise
19.3%. Response-aligned Ridge/MLP on DepMap recover beta in held-out lines
(r = 0.25 to 0.39) but no model, including STATE and MORPH, recovered gamma
without measured perturbations from the target context; adding target
control expression did not help. Gamma became learnable with 30% of target
perturbations measured.

### Revised research questions

- **A.** How much context-specific response is identifiable under truly
  zero-shot context shift?
- **B.** Can richer context-conditioned biological priors recover any
  predictable portion of gamma beyond the sources already tested?
  (High-risk hypothesis; prior evidence negative.)
- **C.** When gamma cannot be predicted, can we identify perturbation/context
  pairs for which beta is sufficiently transferable?

Added direction: predict whether a perturbation is likely conserved versus
strongly context-dependent in an unseen context (transferability
prediction), which is useful even if gamma itself stays unpredictable.

### Know / Suspect / Don't understand / Disproved (update)

**Know**: the decomposition and the zero-shot negative result above [P];
beta transfers; the template is largely in control-expression PCs but the
residual is not.

**Suspect**: the identifiable part of gamma zero-shot is small but not zero
for a minority of perturbations with lineage-specific biology (e.g. GATA1 in
K562); those may be flagged by gene priors even if their gamma vectors
cannot be predicted.

**Don't understand**: whether Molina and Zhang's control-expression test
(PCs, MORPH fine-tuning on target controls) exhausts what basal state can
say; whether the noise floor (19.3% overall, 37.4% after template removal)
hides recoverable gamma in low-cell-count perturbations.

**Disproved (by the literature, not by us)**: the original assumption that
basal state x gene priors would recover gamma zero-shot with existing
information sources.

### Next experimental milestone (Phase 3 gate, revised)

Reproduce the Molina/Zhang decomposition on our own pseudobulk pipeline for
the four lines, with the noise correction, and report per-component variance
and the zero-shot ceiling for gamma. Add the beta-only and response-aligned
Ridge/MLP baselines and score them with cell-eval2 on a leave-one-line-out
split, so we know what a beta-only Arc submission is worth before any
interaction model is built. Do not implement the interaction model until this
milestone is complete.

---

## 2026-09-18 — Read-only audit of the official Arc 2026 validation controls

### Context

The official Arc Virtual Cell Challenge 2026 validation control bundle was
downloaded to `data/raw/arc2026/controls/` (six files, provenance and SHA-256
digests under `data/provenance/`). This entry covers a strictly read-only audit
of that bundle. No model was trained, no raw file was modified, and A/B/C were
not used as perturbation-response training data. Full report:
`reports/arc2026_controls_audit.md`.

### What was implemented

- `src/virtual_cell/data/arc2026.py` — the official bundle as a module: layout
  constants and `bundle_paths` / `missing_bundle_files`; `Arc2026Manifest` with
  `load_manifest`, `load_gene_names`, `load_pert_counts`; memory-safe HDF5
  access (`stream_row_chunks`, `read_obs`, `read_var_names`, `read_shape`,
  `x_encoding`, `x_dtype`) so a 1e8-entry CSR matrix is never densified;
  `audit_context_file` producing a `ContextAudit` in one streaming pass;
  `audit_summary_table`, `depth_quantile_table`, `basal_mean_table`,
  `top_basal_difference_genes`, `cell_id_overlaps`, `panel_composition`,
  `subsample_context`; `check_cross_context_invariants` returning one
  `InvariantCheck` per invariant so all failures are reported at once;
  `same_labels`, `sha256sum`, `parse_checksum_file`.
- `scripts/audit_arc2026_controls.py` — reproducible audit entry point. Prints
  the full report, writes tables and four figures to
  `outputs/arc2026_controls_audit/` (git-ignored), and exits non-zero naming
  the failures if any invariant fails. Runs in about 25 s.
- `tests/test_arc2026_controls.py` — 29 tests, skipped when the git-ignored
  bundle is absent. Suite is now 60 tests (was 31).

### Result: all 44 invariants passed

A/B/C are each 18,400 cells x 18,533 genes, CSR float32, no explicit stored
zeros, counts finite, non-negative and integer-valued (max 977 / 1,753 /
2,462). Identical gene sets in identical order, matching `gene_names.csv`
element-for-element. Exactly one context label per file, control cells only
(`target_gene == 'non-targeting'`), 46 shared non-targeting guides x 400 cells,
unique and cross-context-disjoint cell ids. The manifest is self-consistent:
`46 x 400 = 18,400 = control_cells` and
`control_cells + 300 x 400 = 138,400 = ground_truth_cells`, so the hidden
ground truth includes the control cells.

Median depth is comparable across contexts (about 20,000 UMIs, 5,500-6,100
genes detected per cell).

| | A | B | C |
|---|---|---|---|
| sparsity | 0.6777 | 0.7022 | 0.6827 |
| library size median / mean | 20,109 / 21,134 | 19,946 / 19,990 | 20,034 / 21,157 |
| library size min | 3,275 | 710 | 3,447 |
| genes detected median | 6,147 | 5,756 | 6,006 |

Basal pseudobulk Pearson: A-B 0.688, A-C 0.602, B-C 0.732 (Spearman 0.788 /
0.758 / 0.851). PCA of control cells separates the three contexts completely
(PC1 26.8%, PC2 15.3%).

### Discrepancies and risks found

- **One invariant failed on the first run and was a false alarm.**
  `var_names == gene_names.csv` reported FAIL because the `.h5ad` `var` index
  is pandas `string` dtype while the CSV loads as `object`, and
  `Index.equals` is False across those dtypes. All 18,533 labels were verified
  identical in identical order. Fixed with `arc2026.same_labels`, which
  compares values and order only and still rejects reordering, case changes and
  whitespace (tested).
- **[Know] The 18,533-gene panel excludes the highest-abundance transcript
  families**: zero `RPL*`/`RPS*`, zero `MRPL*`/`MRPS*`, zero `MT-RNR*`; GAPDH,
  EEF1A1, PTMA, MALAT1, NEAT1, LDHA, XIST also absent; the 12 protein-coding
  `MT-` genes are kept. The profile is therefore much flatter than raw 10x
  (top 10 genes hold 3.0-4.7% of the library) and the mitochondrial fraction is
  0.30-0.60%. **Absolute expression is not comparable to unfiltered public
  data**, and the Replogle/Nadig intersection will lose the ribosomal block
  entirely. Pinned by a test so a changed panel fails loudly.
- **[Know] The counts are real scRNA-seq, not a simulation.** Checked because
  the flat profile above is unusual. Per-gene variance/mean has median
  1.45-1.63 with a tail to 208-377 (overdispersed, not Poisson). The 12 `MT-`
  genes have mean pairwise r = 0.42 across cells versus 0.01 for random
  expressed genes; the HIST1 cluster is likewise elevated. The 300 CRISPRi
  targets are expressed in the controls (median mean-count 0.86-1.23 versus
  0.13-0.20 for a typical panel gene, none zero). Gene symbols carry real
  biology, so gene-keyed priors (DepMap, GO, networks) can be joined on them.
- **[Know] Context B has a low-depth tail that A and C do not**: 458 cells
  (2.49%) below 2,000 UMIs and 178 below 1,000, against zero in A and C. Any
  per-cell QC threshold will remove cells from B and almost none from A or C.
  Choose one threshold, apply it identically, report per-context cell loss.
- **[Don't understand] ACTB is near-absent in A and B** (0.88 and 1.28 mean
  counts per cell) but the top gene in C (120). Not explained by panel
  filtering, since ACTB is in the panel. Recorded as an open observation; no
  identity inference was attempted.
- **[Know] Cells-per-perturbation ambiguity is resolved for the validation
  phase.** `pert_counts.csv` in this release has only a `target_gene` column
  and no count column; the budget is `cells_per_pert = 400` from the manifest,
  uniform across all 300 targets. The earlier conflict between the cell-eval2
  brief ("unconstrained") and the vcc CLI ("exactly 400") resolves in favour of
  400 for this bundle. Re-check for the final D/E/F phase.
- **[Know] The synthetic-vs-Arc label collision is now live.** Synthetic
  contexts are also called A/B/C. They are in `data/raw/synthetic/` and the Arc
  bundle in `data/raw/arc2026/controls/`, and the Arc loaders are a separate
  module, but renaming the synthetic labels (`SYN_A`) is still worth doing.
- **[Suspect, weakened] Nearest-context transfer has little to lean on.** Basal
  Pearson of 0.60-0.73 between validation contexts is below the 0.82-0.85 of
  our own synthetic data and below what is typical between human cell lines,
  and the PCA separation is complete. The contexts are far apart at baseline,
  which raises rather than lowers the difficulty of zero-shot transfer.

### Commands run

```
shasum -a 256 -c data/provenance/arc2026_controls_sha256.txt   # all OK, before and after
uv run python scripts/audit_arc2026_controls.py                # 44/44 invariants passed
uv run pytest -v                                               # 60 passed
uv run ruff check .                                            # All checks passed
uv run ruff format --check .                                   # 19 files already formatted
uv build                                                       # wheel contains both packages
```

### Next steps (ordered, unchanged in intent)

The reproduction gate stands. **No interaction model or other novel
architecture until it is complete.**

1. **Next action:** acquire the four public CRISPRi cell lines used by Molina
   and Zhang — Replogle 2022 (K562, RPE1) and Nadig 2025 (HepG2, Jurkat) —
   pinning dataset versions and recording provenance and SHA-256 digests the
   same way the Arc bundle was, then build the gene intersection against the
   Arc 18,533-gene panel and record how much is lost (expect the ribosomal
   block to drop out entirely).
2. Per-context, per-perturbation pseudobulk on those four lines; frozen
   leave-one-context-out split files under `data/splits/`.
3. Differential expression matching Arc's definition, plus a pseudobulk DE for
   the research track; the six Arc metrics locally via `cell-eval2`.
4. Baselines 0 to 5 from `plans.MD`; first leave-one-context-out benchmark table.
5. Reproduce the Molina/Zhang decomposition with the noise correction and
   report per-component variance (template / beta / gamma / noise) and the
   zero-shot ceiling for gamma. This is the Phase 3 gate.
