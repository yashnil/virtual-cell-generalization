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

---

## 2026-09-18 (later) — Molina & Zhang reproduction gate: Task 1 done, Tasks 2–3 blocked

### Context

Started the mandatory reproduction gate. Audited the authors' official
repository, implemented our own decomposition, and prepared the raw-data
acquisition plan. Full frozen specification:
`reports/molina_zhang_reproduction_spec.md`. Arc A/B/C were not touched.

### Headline: the authors' processed data does not exist publicly

`xinyizhanglab/perturbation-decomposition` @
`a15214780619736d393f40240e56ba992fd416a3` (retrieved 2026-09-18) is 189 KB of
source only. Its README claims "Data is included (processed pseudobulk + DepMap
embeddings, ~90 MB)" and "All predictions are included in `results/` (~150 MB)",
but `.gitignore` excludes `data/processed/`, `data/raw/`, `*.pkl` and
`results/`. There are 0 releases, 0 tags, 0 forks, an empty wiki, and the
organisation has exactly one repository. No Zenodo or figshare deposit was
found. **There is no reference fixture**, so Task 2's numeric verification and
all of Task 3 cannot be run as specified.

Compounding this: 34 of ~45 Python files hardcode
`/mnt/storage01/home/amolina/benchmark_2026`, and `data/prepare_release.py` is
documented "Run ON THE CLUSTER" and additionally reads a private sibling tree
`agop-perturbation/morph_pipeline/` that is not public. The DepMap embeddings
ship from that private pipeline.

### [Know] The authoritative decomposition is in a figure script, not the module

`decomposition/anova.py` — the entry point the README documents — is **not** the
paper's decomposition. The paper's numbers come from `figures/fig1_f.py`. They
differ on the denominator (`np.var(D)` vs uncentred `mean ||D[c,p]||²`), on
whether `mu` is reported at all, on cell-line ordering (sorted vs declaration
order), and critically `anova.py` has **no noise correction**. Its three
fractions do not sum to 1. Running the documented command cannot produce
27.8 / 29.4 / 23.5 / 19.3 even with the data in hand.

Decomposition as actually implemented: plain balanced two-way ANOVA by cell
means; `template` is `mu + alpha` combined; noise is a **residual** left after
subtracting four reproducible components estimated by 50 split-half resamples
(`RandomState(42)`) using **cross-half dot products**, with controls *not* split
and an odd cell dropped. Template removal for the second bar is **projective**,
not subtractive.

### [Know] Defects in the released code

- `evaluation/metrics.py::evaluate_predictions` calls `centroid_accuracy` and
  `median_rank`, neither defined anywhere → `NameError`. This is the README's
  own Quick Start snippet.
- `configs/datasets.yaml`: the `nadig_hepg2` block is mis-indented and parses to
  `None`, leaking 8 keys as siblings. HepG2 is unconfigurable as shipped.
- `configs/datasets.yaml`: the Replogle `raw_url` (figshare 21999546) points to
  an unrelated PLoS ONE figure, *"Characteristics of the Greifer and the
  Axon-Hook"*, and the same wrong ID is given for both K562 and RPE1. Correct
  record is figshare+ 20029387.
- `process.py` writes `{cl}_pseudobulk.pkl`; `anova.py` reads `{cl}.pkl`.
- Min cells per perturbation is 5 in the data scripts but **10** in the figure
  that produces the paper's numbers.

### [Don't understand] The upstream preprocessing is absent

Both data scripts **consume** `obsm["X_hvg"]`; nothing in the repository creates
it. Normalisation, log transform, HVG flavour and batch handling are all
unspecified (`configs` says only `n_hvgs: 2000`). The paper's per-line cell and
perturbation counts (K562 188,590 cells / 1,383 perts; RPE1 173,737 / 1,499;
HepG2 96,616 / 1,340; Jurkat 184,470 / 1,537) are also much smaller than the raw
deposits, so an undocumented filtering step sits in between. A byte-faithful
raw-data reproduction is therefore not currently possible.

### What was implemented

- `src/virtual_cell/decomposition/anova.py`: our own implementation, faithful to
  `fig1_f.py`. `build_response_tensor` (shared-perturbation intersection, frozen
  orders, balance checks), `decompose`, `total_sum_of_squares`,
  `sums_of_squares`, `uncorrected_fractions`, `project_out_template`,
  `cross_half_signal`, `split_half_delta_tensors`, `split_half_signal`,
  `noise_corrected_fractions`, `beta_fraction_per_perturbation`.
- `tests/test_decomposition.py`: 41 tests. Exact reconstruction; all zero-sum
  conditions; SS partition and pairwise orthogonality; **planted-component
  recovery** (components and SS budget planted in synthetic tensors and
  recovered exactly, establishing correctness without the authors' data);
  permutation invariance; split-half unbiasedness, monotone noise growth as
  cells per perturbation fall, and seed reproducibility; frozen perturbation /
  cell-line / gene sets; rejection of unbalanced, malformed and non-finite
  input; and **verbatim equivalence with the reference algebra**, by
  transcribing `fig1_f.py` L78–90, L94–104, L131–144 into the tests and
  asserting agreement to 1e-12.

Nothing was tuned toward the paper's percentages; those remain unverified.

### Raw data acquisition plan (verified, not started)

| cell line | source | exact file | size |
|---|---|---|---|
| K562 | figshare+ 20029387 | `K562_essential_raw_singlecell_01.h5ad` | 10.66 GB |
| RPE1 | figshare+ 20029387 | `rpe1_raw_singlecell_01.h5ad` | 8.70 GB |
| HepG2 | GEO GSE264667 | `GSE264667_hepg2_raw_singlecell_01.h5ad` | 5.2 GB |
| Jurkat | GEO GSE264667 | `GSE264667_jurkat_raw_singlecell_01.h5ad` | 8.7 GB |
| DepMap 2024Q2 | DepMap portal | `CRISPRGeneEffect.csv` | ~0.1 GB |

Total ≈ 33.4 GB. GEO filenames confirmed against the FTP listing; figshare
filenames, sizes and stable file endpoints confirmed against the figshare API.
**Awaiting approval before downloading.**

### Commands run

```
uv run pytest -v          # 101 passed
uv run ruff check .       # All checks passed
uv run ruff format --check .
uv build                  # wheel contains both packages
```

### Gate status

**NOT passed.** Task 1 complete, Task 2 implemented and mathematically verified
but numerically unverified, Task 3 blocked, Task 4 complete and awaiting
approval. The hard rule stands: no interaction model or other novel architecture
until the gate passes.

### Next steps (ordered)

1. **Decision required from the user** on how to unblock the reference data.
   Options: (a) email the corresponding authors for `data/processed/` and
   `results/`; (b) open a GitHub issue on the repo; (c) proceed to raw-data
   acquisition (~33.4 GB) and accept that our preprocessing will be *our own*
   reconstruction, making this an independent re-derivation rather than a
   byte-faithful reproduction; (d) wait for the peer-reviewed version.
2. If (c): download the five files above, reconstruct normalisation + 2,000-HVG
   selection ourselves, document every choice as a deviation, and report the
   decomposition with sensitivity analysis over the unspecified choices.
3. Only then compare against 27.8 / 29.4 / 23.5 / 19.3 and judge the gate.
4. Implement beta-only transfer and beta/gamma-recovery metrics — noting these
   are **our** constructions, absent from the authors' repository.
5. Arc alignment stage, strictly after the gate.

---

## 2026-09-18 (later still) — reproduction track closed; independent four-context track opened

### Decision

Per the user: do **not** contact Arc or the paper's authors, and do **not**
download the ~33.4 GB of raw deposits. The missing artifacts are specific to the
Molina & Zhang reproducibility repository, not to the Virtual Cell Challenge.
Use the public **scPertEval** standardized datasets instead.

Gate amended in two tracks:

- **Exact Molina & Zhang reproduction — BLOCKED, now closed as unachievable
  from public materials.**
- **Independent four-context re-derivation — READY.**

**Naming rule, binding from here on:** this work is the *independent
four-context decomposition*. It is never to be called a Molina & Zhang
reproduction, and its numbers are never to be presented as reproducing theirs.
Recorded in `reports/scperteval_four_context_data_spec.md` §3 and as an
amendment (§7) to `reports/molina_zhang_reproduction_spec.md`.

### Source audit

`Virtual-Cell-Research-Community/scPertEval` @
**`4685f11927e887745737600170da7a655b727553`** (`main`, release v0.2.0,
2026-09-09), retrieved 2026-09-18, MIT. Data in the public read-only bucket
`gs://scperteval/processed/`, also plain HTTPS (all four return HTTP 200).

| dataset | cell line | bytes | GB |
|---|---|---:|---:|
| `replogle22k562` | K562 | 2,430,512,332 | 2.431 |
| `replogle22rpe1` | RPE1 | 1,877,364,555 | 1.877 |
| `nadig25hepg2` | HepG2 | 1,236,448,196 | 1.236 |
| `nadig25jurkat` | Jurkat | 2,004,474,709 | 2.004 |
| **total** | | **7,548,799,792** | **7.549** |

### [Know] Contents verified directly, without downloading

Read the remote HDF5 metadata over HTTP range requests — **24.2 MB fetched of
7,549 MB**, no expression data transferred.

| | K562 | RPE1 | HepG2 | Jurkat |
|---|---|---|---|---|
| cells | 308,646 | 240,774 | 133,757 | 258,202 |
| genes | 8,563 | 8,749 | 9,623 | 8,881 |
| controls | 10,691 | 11,485 | 4,976 | 12,013 |
| perturbations | 1,971 | 2,016 | 1,818 | 2,137 |
| min cells/pert | 30 | 30 | 30 | 30 |

All four: `X` = CSR float32 `log1p(CP10K)`; `obs` has only `perturbation`; `var`
has only the index; `layers`/`obsm`/`uns`/`varm`/`varp` all empty; gene ids are
HGNC symbols; control label is the literal `control`; no `+` combinations.

- **[Know] Raw counts are not retained anywhere** — the raw layer was dropped in
  trimming. Anything needing counts must go back to the original deposits.
- **[Know] The four datasets do not share a gene space** (8,563–9,623). The
  four-way intersection is **6,640 genes**; union 11,909. Each dataset loses
  22–31 % of its genes to the intersection.
- **[Know] Shared perturbations across all four: 1,264** (union 2,374, pairwise
  1,511–1,872). 1,072 of the 1,264 target a gene that is itself in the shared
  response space; 192 do not.
- **[Know] The ≥30 cells-per-perturbation floor is inherited from upstream**, so
  every min-cells threshold from 1 to 30 yields the same 1,264 perturbations —
  M&Z's 5 and 10 are both no-ops here. Raising it costs a lot: 50 → 577,
  100 → 96.
- **[Know] scPertEval's own preprocessing is fully specified**, unlike M&Z's:
  label cleaning → `normalize_total(1e4)` + `log1p` → `filter_cells(min_genes=200)`,
  `filter_genes(min_cells=3)` → trim. **No HVG selection, no scaling, no PCA,
  no batch correction.** This is the decisive advantage over the blocked track.
- **[Know] Corroboration of common ancestry:** scPertEval's measured control
  counts are *identical* to M&Z's reported controls (10,691 / 11,485 / 4,976 /
  12,013). Both derive from the same deposits with the same control definition.
  They diverge downstream — M&Z report 188,590 K562 cells / 1,383 perturbations
  against 308,646 / 1,971 here — via a filtering step M&Z never document. So our
  numbers are **not expected to match theirs**, and a mismatch is not evidence
  of a bug in either.
- Year correction: the HepG2/Jurkat source is **Nadig et al. 2025**
  (Nat. Genet. 57:1228–1237, doi 10.1038/s41588-025-02169-3), not "2024" as
  M&Z's config says.

### Storage preflight

`df -h ~` → 1.8 Ti total, 620 Gi used, **1.2 Ti available**. The 7.549 GB
download is **0.6 % of free space**. No decompression copy needed. **Safe.**

### Proposed pipeline (full detail in the spec)

Frozen 1,264 × 6,640 balanced design, committed to `data/splits/` before any
analysis. `ctrl_mean[c]` = mean over that line's `control` cells on the shared
genes; `delta[c,p]` = perturbation mean − control mean; **no further transform**
(X is already log1p(CP10K)). All 6,640 shared genes used initially; HVGs belong
in the sensitivity analysis only, since undocumented HVG selection is precisely
what made M&Z irreproducible. Decomposition via
`virtual_cell.decomposition.anova` on the uncentred `mean ||delta||²`
denominator, reporting `mu` and `alpha` both separately and merged. Split-half:
disjoint equal halves within each (line, perturbation), both referenced to the
**full** control mean, 50 resamples, fixed seed; feasible everywhere since the
smallest perturbation has 30 cells.

Leakage rule recorded now for the later transfer stage: intersections use
**identity only, never expression**, so freezing them is safe; but any step that
looks at expression (HVG, PCA, scaling) must be fitted on **source lines only**.

### Gate criteria (do not require 27.8 / 29.4 / 23.5 / 19.3)

1. decomposition mathematics passes all invariants;
2. shares stable across preprocessing choices, tolerance declared in advance;
3. beta and gamma quantified with uncertainty and non-degenerate;
4. split-half separates reproducible signal from noise (share in (0,1),
   monotone in cells per perturbation, stable across seeds);
5. no single preprocessing choice flips a qualitative conclusion.

### Commands run

```
uv run pytest -v          # 101 passed
uv run ruff check .       # All checks passed
uv run ruff format --check .
uv build                  # wheel contains both packages
```

No new code was written this session; the decomposition module and its 41 tests
from the previous entry are unchanged and are what the new track will use.

### Next step

**Awaiting download approval.** On approval: fetch the four files (7.549 GB) to
`data/raw/scperteval/`, verify against the recorded MD5s, write provenance to
`data/provenance/` in the same form as the Arc bundle, freeze the 1,264 × 6,640
split files, then run the decomposition and the §5 sensitivity battery.
Still no novel architecture until the gate passes.

---

## 2026-09-19 — Independent four-context decomposition v1 (canonical result)

### What was done

Downloaded the four approved scPertEval datasets and ran the frozen protocol
end to end. Full report: `reports/four_context_decomposition_v1.md`.

**Protocol preserved first.** Before any download,
`data/provenance/scperteval/protocol_freeze.txt` pinned SHA-256 digests of both
spec documents, `virtual_cell/decomposition/anova.py` and
`tests/test_decomposition.py`. Re-verified after the run: **all four match**, so
the decomposition protocol was not altered in response to results.

Downloads: all four files at exact audited byte sizes, **all four upstream MD5s
match** the values recorded before download; local SHA-256 computed. Provenance
in `data/provenance/scperteval/` (per-dataset JSON + manifest + checksum file).
scPertEval @ `4685f11927e887745737600170da7a655b727553`, retrieved 2026-09-19.

### [Know] The frozen design reproduces the remote audit exactly

Recomputed locally, independently of the HTTP metadata audit: **1,264 shared
perturbations, 6,640 shared genes** — identical to the audited expectation.
Committed to `data/splits/four_context_v1/` with a manifest recording counts,
ordering rules, source hashes and the explicit statement that intersections use
identifier presence only.

### [Know] Canonical decomposition, delta tensor (4, 1264, 6640)

All invariants pass before interpretation: reconstruction 4.44e-16, zero-sum
≤2.6e-12, **SS partition exact (rel err 0.0)**, max cross term 2.53e-17.
Total response energy 41.4168.

| component | uncorrected | noise-corrected (50 resamples, seed 42) |
|---|---:|---:|
| μ | 12.75 % | — |
| α | 7.54 % | — |
| template (μ+α) | 20.29 % | **20.27 %** (sd 0.010) |
| β conserved | 37.22 % | **30.07 %** (sd 0.016) |
| γ interaction | 42.49 % | **21.05 %** (sd 0.033) |
| noise | — | **28.62 %** (sd 0.044) |

Resample uncertainty is negligible (<0.12 pp range on every share).

### [Know] Per-component reproducibility is the key result

Cross-half signal ÷ raw SS: **μ 100.0 %, α 99.8 %, β 80.8 %, γ 49.5 %**
(overall 71.4 %). The template is essentially noise-free, β is largely
reproducible, and **γ is about half signal and half noise** — which is exactly
why uncorrected γ (42.5 %) halves to 21.1 %. γ carries 75 % of all measurement
noise in the decomposition. Skipping the noise correction would roughly double
the apparent interaction.

- **β is materially non-zero**: 30.07 % of energy, 80.8 % reproducible.
- **γ is materially non-zero**: 21.05 %, comparable to the template, with median
  ‖γ_{c,p}‖ (3.3–4.3) of the same order as median ‖β_p‖ (3.02).
- **γ is reproducible, but only about half of it is.** Real and substantial, yet
  the noisiest component and untrustworthy per-(context, perturbation) without
  accounting for depth.

### [Know] Simpson's paradox in the reliability-vs-depth relationship

Marginal Spearman(cells, reliability) is only **0.110**, which initially looked
like a failure of gate criterion 4. It is Simpson's paradox: ‖δ‖ correlates with
reliability at 0.756 but **negatively with cell count at −0.492**, because
`E‖δ̂‖² = ‖δ‖² + noise/n` inflates measured effect magnitude at low depth.
Stratified by effect-size quintile, Spearman(cells, reliability) is **+0.48 to
+0.92 in every stratum** and median reliability rises monotonically across
cell-count quintiles within each. **Criterion 4 is met once effect size is
controlled for.** Worth remembering: raw ‖δ‖ is depth-biased upward, so any
future filtering or ranking on effect magnitude must not use it naively.

### [Know] Other findings

- **Correction to our own spec:** §1.4 predicted `expm1(X)` row sums slightly
  below 1e4 (gene filtering after normalisation). They are **exactly 10000.0**
  in all four datasets, so filtering preceded normalisation. Harmless direction;
  protocol unaffected.
- Reliability is **strongly bimodal** (mode near 0, second near 0.75) — expected
  in an essential-gene screen where many perturbations are near-null. Overall
  median 0.313; only 25.7 % (K562) to 57.1 % (RPE1) exceed r = 0.5.
- **RPE1 is cleanest** (median reliability 0.570) yet has the **largest** median
  ‖γ‖ (4.26) — higher reliability with larger γ argues its interaction is
  genuinely biological.
- **HepG2 is thinnest** (median 57 cells/pert) and contributes disproportionate
  noise, as the spec predicted.
- **Transferability is narrowly distributed**: β fraction median 0.412, only
  21.4 % β-dominated (>0.5), 8.5 % below 0.3. Most perturbations mix conserved
  and context-specific response rather than being cleanly one or the other.

### [Don't compare] Molina & Zhang

Their template 27.8 / β 29.4 / γ 23.5 / noise 19.3 % is recorded in the report
**once, as contextual reference only**. Different cell set, different
undocumented 2,000-HVG response space, different filtering. Our β (30.1) and γ
(21.1) landing nearby is interesting but is **not** reproduction, and the
differences are **not** evidence of a bug in either analysis. The naming rule
holds: this is the *independent four-context decomposition*.

### Performance

Runtime **9.4 min** end to end (pseudobulk ~105 s for all four contexts;
split-half ~7 min). Peak RSS **21.2 GB** of 68.7 GB — dominated by the 13.4 GB
disk-backed half-mean memmap's page cache, deleted after use. No full matrix was
ever densified and only one context's cell matrix was resident at a time.

### Gate status — NOT yet passed

| # | criterion | status |
|---|---|---|
| 1 | mathematics passes all invariants | **PASS** |
| 2 | stable across preprocessing choices | **NOT TESTED** (v1 scope) |
| 3 | β and γ quantified | **PASS** |
| 4 | split-half separates signal from noise | **PASS** (monotone within effect strata) |
| 5 | no pathological preprocessing dependence | **NOT TESTED** |

Criteria 2 and 5 need the sensitivity battery, deliberately out of v1 scope.
**No novel prediction architecture until the gate passes.**

### Predeclared for the sensitivity stage (recorded, not run)

A secondary check must **independently split control cells too**. The primary
scheme reuses one full-context control mean in both halves, inducing correlated
error that can bias the noise estimate downward. Predeclared before the run;
deliberately not executed, so the primary result stands unmodified.

### Commands run

```
./scripts/download_scperteval.sh                       # 4/4 OK, exact byte sizes
uv run python scripts/scperteval_provenance.py         # 4/4 MD5 OK
uv run python scripts/build_four_context_decomposition.py   # 9.4 min, exit 0
uv run pytest -v                                       # 124 passed
uv run ruff check . && uv run ruff format --check .
uv build
```

### Next step

Stop and review the canonical result before anything else. Then the sensitivity
battery for criteria 2 and 5 (min-cells 30/50/100, gene space, response scaling,
3-context subsets, 50 vs 200 resamples, seeds, and the control-split secondary
scheme), with the tolerance declared in advance.

---

## 2026-09-19 (later) — Predeclared robustness battery; independent gate PASSES

### What was done

Froze canonical v1 (`data/provenance/scperteval/canonical_v1_freeze.txt`, 12
digests covering the report, frozen design, derived arrays, summary and the two
code files) **before** running the battery, then ran variants A–E. Full report:
`reports/four_context_decomposition_sensitivity.md`. Canonical v1 was not
modified; both freezes re-verify (12/12 and 4/4) after the run.

Battery code lives in a **separate** module (`virtual_cell.analysis.robustness`)
precisely because `scperteval.py` is frozen. At canonical settings the new,
independent code path reproduces canonical v1 to 0.01 pp — two implementations
agreeing is itself a useful check.

Runtime 57.8 min, peak RSS 23.2 GB of 68.7 GB.

### [Know] All three canonical conclusions survived every variant

Across all 21 variant x feature-space combinations:

| quantity | min | max |
|---|---:|---:|
| beta share | **27.98 %** | 30.79 % |
| gamma share (corrected) | **20.55 %** | 22.80 % |
| gamma share (uncorrected) | 38.51 % | 45.05 % |
| beta reproducibility | 77.4 % | 85.0 % |
| gamma reproducibility | 45.6 % | 57.7 % |

Invariants everywhere: reconstruction <= 4.44e-16, SS partition <= 5.15e-16,
zero-sum <= 1.07e-13. Uncorrected gamma is roughly **double** corrected gamma in
every single variant.

### [Know] A — control-estimation error can only touch the template

Independently splitting the control cells changed template by **-0.246 pp** and
noise by **+0.246 pp**, and changed **beta and gamma by exactly 0.000 pp**.
Reproducibility: mu -0.48 pp, alpha **-2.46 pp**, beta and gamma **0.000 pp**.

This is algebraically necessary, and worth remembering: the control profile
enters `delta[c,p] = pert_mean[c,p] - ctrl[c]` as a term depending on context but
**not** perturbation. In `beta_p = mean_c delta - mu` the two `mean_c ctrl` terms
cancel; in `gamma = delta - mu - alpha_c - beta_p` the `-ctrl[c]` in delta
cancels against the one inside `alpha_c`. A per-context constant lives entirely
in the mu+alpha subspace, orthogonal to beta and gamma.

**So reusing the control mean in the canonical scheme could not have inflated
beta or gamma.** The concern was real for the template and provably void for the
two components the science rests on.

### [Know] B — depth causally improves reliability; the v1 paradox is resolved

Same **643 fixed (context, perturbation) pairs** (all with >=200 cells)
re-estimated at n = 15/30/50/100 cells per half, 25 repeats each. Nothing
conditions on observed ||delta||.

| n | median reliability | 95 % CI |
|---:|---:|---|
| 15 | 0.0966 | [0.081, 0.111] |
| 30 | 0.1762 | [0.149, 0.205] |
| 50 | 0.2626 | [0.226, 0.291] |
| 100 | **0.4188** | [0.373, 0.455] |

Adjacent CIs do not overlap. Within-pair: 15->100 median +0.3010, **99.53 %
improve**; 90.8 % (584/643) increase at every step; all four contexts improve
(HepG2 +0.370 with 100 % improving). Within-repeat sd falls monotonically
(0.037 -> 0.018). The weak marginal correlation in v1 was indeed the
depth-magnitude confound; this design removes it rather than conditioning on the
biased norm.

### [Know] C — feature space

Control-derived global HVG rule (per-context variance over **control cells
only**, averaged across contexts, ranked descending; one ranking for all four
contexts; no perturbation response consulted).

| genes | template | beta | gamma | noise | beta repro | gamma repro |
|---:|---:|---:|---:|---:|---:|---:|
| 6,640 | 20.27 % | 30.07 % | 21.04 % | 28.62 % | 80.8 % | 49.5 % |
| 4,000 | 22.38 % | 29.69 % | 21.79 % | 26.13 % | 82.0 % | 52.7 % |
| 2,000 | 25.25 % | 30.79 % | 22.22 % | 21.74 % | 85.0 % | 57.7 % |

beta moves 1.1 pp and gamma 1.2 pp across a 3.3x change in feature-space size.

### [Know] D — aggregation order

`log(mean(CP10K))` is the noisier estimator: on all 6,640 genes it moves 2.1 pp
of beta and 0.5 pp of gamma into noise (beta 27.98 %, gamma 20.55 %, noise
32.39 %). Direction is what Jensen predicts (linear-scale averaging gives
high-count cells more leverage). The gap narrows to 0.7 pp of beta at 2,000 HVGs.
Both conclusions survive under either order.

### [Know] E — seed stability

Five seeds, all 6,640 genes: ranges of **0.0104 / 0.0072 / 0.0057 / 0.0199 pp**
for template / beta / gamma / noise. Seed choice is irrelevant.

### Gate: ALL SIX CRITERIA PASS

Invariants valid; beta non-degenerate (never < 27.98 %); gamma non-degenerate
after correction (never < 20.55 %); conclusions independent of any single
preprocessing choice; controlled subsampling confirms depth improves reliability;
no variant reverses the interpretation. **The independent four-context
decomposition gate is met.** It does not become a Molina & Zhang reproduction.

### [Don't understand / constraints for what comes next]

1. **Absolute per-perturbation reliability is low.** Even at n=100 median r is
   0.42; at observed median depth most perturbations are well below. Aggregate
   shares are stable to ~0.01 pp, but an individual (context, perturbation)
   response is often poorly determined. Future per-perturbation modelling must
   carry a reliability weight or filter.
2. **gamma's ~50 % reproducibility is the binding ceiling on any gamma
   predictor.** Evaluated against raw gamma, achievable correlation is bounded
   well below 1 regardless of model quality. Evaluate against the
   noise-corrected/reliable portion, or report the ceiling alongside.
3. HepG2 stays the weak context (only 43 pairs with >=200 cells).
4. **||delta|| is depth-biased upward** — never use it as an effect-size filter
   or ranking key without a depth correction.

### Commands run

```
uv run python scripts/run_four_context_sensitivity.py     # 57.8 min, exit 0
uv run python scripts/plot_four_context_sensitivity.py
uv run pytest -v ; uv run ruff check . ; uv run ruff format --check . ; uv build
```

### Next step

Stop for review. No predictive model was built — no Ridge, MLP, gamma predictor,
D predictor, transferability classifier or generative model.

---

## 2026-09-19 (later) — Zero-shot recoverability diagnostic v1

Froze all decomposition-phase artifacts first
(`decomposition_phase_freeze.txt`, 21 digests; re-verified 21/21 after). Full
report: `reports/zero_shot_recoverability_v1.md`. Runtime 7.4 min, peak RSS
18.3 GB. **No model was built.**

### [Know] Two exact identities govern source-only transfer

With `A[p] = mean_S delta[c,p]` over the three sources:

1. `delta[c*,p] - A[p] = (4/3)(alpha_{c*} + gamma[c*,p])`; centred over
   perturbations it is **exactly `(4/3) gamma[c*,p]`**. So "predict the transfer
   residual" and "predict gamma" are the same problem — and gamma must never be
   a predictor input.
2. `A[p]` carries `-gamma[c*,p]/3`, so conserved transfer is **mildly
   anti-correlated with the held-out interaction by construction**.

Both proved as tests on planted components.

### [Know] Conserved transfer: right direction, wrong magnitude

Median per-perturbation Pearson **0.306** pooled (0.277 K562, 0.357 RPE1, 0.326
HepG2, 0.283 Jurkat). Basal-weighted combinations are indistinguishable from the
plain source mean (0.310); **nearest-context is clearly worse (0.239)** —
averaging beats picking with only three sources.

**But response energy explained is NEGATIVE** (−0.138 pooled; −0.467 K562,
+0.098 RPE1, −0.065 HepG2, −0.184 Jurkat): as a point prediction of raw response
it is worse than predicting zero. The held-out template `alpha_{c*}` is
unrecoverable from source responses and the source mean substitutes
`-alpha_{c*}/3`. **Any model scored on magnitude — Arc's metrics included — must
fix template and scale, not just direction.**

### [Know] Reliability normalisation, derived and validated

`corr(h1,h2) = rho_half`; a perfect latent predictor reaches `sqrt(rho_half)`,
**not** `rho_half`; `rho_latent = r_obs / sqrt(rho_half)`;
`rho_full = 2 rho_half/(1+rho_half)`. Validated on synthetic latent+noise:
reliability recovers the planted ratio to ±0.02, the ceiling matches `sqrt(rho)`
and differs from `rho` by >0.05, and disattenuation recovers planted
correlations of 0.3/0.6/0.9 to ±0.05.

Target reliabilities — K562 rho_half 0.241, RPE1 **0.571**, HepG2 0.226, Jurkat
0.276. Reliability-normalised conserved transfer: **0.49–0.64**, i.e. about half
the attainable latent correlation.

### [Know] Gamma is recoverable zero-shot ONLY where a similar partner exists

| held out | r(gamma_true, gamma_hat) | ceiling | frac > 0 |
|---|---:|---:|---:|
| K562 | **+0.185** | 0.607 | 85.8 % |
| Jurkat | **+0.217** | 0.585 | 93.3 % |
| HepG2 | +0.033 | 0.534 | 59.9 % |
| RPE1 | +0.002 | 0.743 | 50.8 % |

Explanation: cross-context gamma correlation against the **forced null of
-1/(C-1) = -0.333** (because `sum_c gamma = 0`) shows exactly one pair above
null — **K562–Jurkat, excess +0.196**, also the most basally similar pair
(0.934). Spearman(basal similarity, excess) = +0.714 over 6 pairs, but it rests
on one point: hypothesis, not finding. K562 and Jurkat each have a
gamma-correlated partner in their source set; RPE1 and HepG2 do not.

**"gamma exists and is reproducible" != "gamma can be predicted zero-shot".**
It is context-dependent, and four contexts cannot say more.

### [Know] Source agreement is the actionable transferability feature

Spearman vs conserved-transfer success: **source agreement +0.828 raw, +0.726
after reliability normalisation** — and it is computable at inference from source
contexts alone. Target reliability drops +0.864 → +0.472 under normalisation
(most of its raw association was the noise artefact). Target-gene basal
expression +0.094, basal context distance −0.036, cells/pert −0.023: all useless
per-perturbation.

Three classes found: **transferable** (POLRMT, SMG5, TFAM; normalised r
0.86–0.92); **reliably context-specific** (BCR in K562 at rho 0.917 but r
−0.074; NSMCE2/EXOSC1/SHQ1/GRWD1 in RPE1; GAB2, IPO7) — all with near-zero
source agreement, so flaggable in advance; and **large but unreliable** (PIAS4,
PSMC1, ANAPC1... all Jurkat, 31–36 cells, ‖delta‖ ~6, reliability ~0.00–0.05) —
the depth-biased ‖delta‖ artefact made concrete.

### Bug found and fixed mid-run (recorded because it nearly passed)

The first gamma-reliability implementation replaced only the *target* row with a
half while leaving the three sources at full data. Since `gamma[c,p]` depends on
all four contexts, the source contribution was identical in both "halves" and
reliability came out at **0.90–0.94** — irreconcilable with the canonical 49.5 %
gamma reproducibility, which is what exposed it. It also averaged gamma across
repeats before correlating. Fixed to split **every** context and average
per-repeat correlations; corrected values are 0.17–0.38. **Lesson: any
split-half on a quantity derived from all contexts must split all of them.**

### Decision point

Ranked by evidence, not novelty:

1. **Transferability / D prediction — strongest.** Source agreement already
   gives +0.726 normalised association using only inference-time information,
   with no model. Well-posed, evaluable, immediately useful.
2. **Template and scale calibration — underrated, arguably first.** Negative
   energy explained is the largest failure mode found; the target's own control
   cells are available at inference; cheap and high leverage.
3. **Pathway-level gamma — plausible, untested.** Gene-level gamma sits near its
   noise floor; aggregation should raise reliability.
4. **Exact gene-level gamma — weakest.** Works in 2 of 4 folds, recovers a third
   to half of a ceiling that is itself only 0.53–0.74, with four contexts.

Recommended: (1) with (2), and (3) as a cheap diagnostic before committing to (4).

### Caveats

Four contexts / six pairs — do not fit anything of appreciable capacity to
context. K562–Jurkat are both suspension leukaemia lines *and* share a dataset of
origin (Replogle vs Nadig), so the "similar partner" effect is confounded and
cannot be separated here. Disattenuation is undefined below rho 0.05, excluding
13–33 % of pairs per fold. Reliability used 10 repeats, not the canonical 50.

### Next step

Stop for review before any modelling.

---

## 2026-09-20 — Transferability foundations v1

Froze zero-shot artifacts first (`zero_shot_v1_freeze.txt`, 14 digests). All
four freeze records re-verify after: 4/4, 12/12, 21/21, 14/14. Full report:
`reports/transferability_foundations_v1.md`. Runtime 0.8 min. **No model
trained.** Every new estimator has a leakage test that scrambles the held-out
response matrix and requires bit-identical output.

### [Disproved] Basal control expression cannot recover the context template

Three independent signals agree: gene-wise r(alpha, basal deviation) is ≈0 and
**flips sign** across folds (−0.056, −0.046, −0.263, +0.065); only **0.15–8.2 %**
of ‖alpha‖² lies in the span of the source basal deviations; and the fitted
global scalar k changes sign across folds (−0.0071, −0.0078, **+0.0025**,
−0.0134). Only a weak magnitude association survives (Spearman 0.27–0.34
between |alpha| and |basal dev|) — no usable direction.

Performance confirms it: `direct_basal` is catastrophic (energy −6 to −17),
`global_scalar` changes nothing, and `ridge_subspace` helps in three folds and
**hurts in HepG2** — a concrete demonstration of over-capacity with three
source points.

**This overturns the hypothesis I carried in from the zero-shot phase.** I
expected template offset to be the main lever for negative energy explained. It
is not.

### [Know] Scale calibration is the lever, not template offset

| estimator | K562 | RPE1 | HepG2 | Jurkat |
|---|---:|---:|---:|---:|
| source mean (zero) | −0.467 | +0.098 | −0.065 | −0.184 |
| **scale_only** (available) | **−0.211** | **+0.189** | **+0.094** | **−0.012** |
| oracle_alpha (ceiling) | −0.424 | +0.378 | −0.058 | −0.133 |
| scale + oracle_alpha (ceiling) | **+0.046** | **+0.425** | **+0.119** | **+0.117** |

**Even the oracle template leaves energy negative in 3 of 4 folds** — with true
alpha the residual is exactly `(4/3) gamma`, so the template was never the
dominant error. A single scalar shrinkage fitted by inner leave-one-source-out
(**0.455, 0.440, 0.428, 0.469** — strikingly consistent) lifts every fold, turns
HepG2 positive and Jurkat to ≈0, and is fully inference-available. Shrinkage is
≈0.44 rather than 1 because the source mean carries `beta − gamma/3` and cannot
predict gamma at all, so the variance-optimal prediction is heavily shrunk.

### [Know] Source agreement validated per fold and under confounding

Spearman vs reliability-normalised success, **each fold independently**: K562
+0.729, RPE1 +0.838, HepG2 +0.622, Jurkat +0.704, all with bootstrap CIs far
from zero. Joint partial controlling source magnitude, source reliability,
source cells and target-gene basal: **+0.525, +0.301, +0.142, +0.547**. It does
add information beyond detecting weak/noisy perturbations, but the margin is
modest and weakest in the noisiest context.

**Analysis error caught and corrected:** the first confound set included
`source_min_pair_agreement`, which drove the partial to +0.12–0.19. That is the
*minimum* of the same pairwise correlations whose *mean* is the exposure —
over-adjustment, not confound control. Reported separately and excluded.

### [Know] Pathway-level gamma is ~3x more recoverable than gene-level

MSigDB 2024.1.Hs, Hallmark + Reactome, SHA-256 recorded, used as released.
Predeclared aggregation: unweighted mean over set genes present in the frozen
space, sets with <10 genes dropped (45 Hallmark, 874 Reactome qualify).
Aggregation is linear so decompose-then-score == score-then-decompose (tested).

Median r(gamma_true, gamma_hat), basal_affine:

| held out | gene | Hallmark | ceiling | normalised |
|---|---:|---:|---:|---:|
| K562 | +0.185 | **+0.570** | 0.874 | **0.652** |
| Jurkat | +0.217 | **+0.510** | 0.852 | **0.599** |
| RPE1 | +0.002 | **+0.182** | 0.937 | 0.194 |
| HepG2 | +0.033 | +0.026 | 0.801 | 0.032 |

Partly a reliability effect (pathway rho_full 0.64–0.88 vs gene 0.28–0.55) but
**ceiling-normalised values also rise sharply** (K562 0.32 → 0.65), so it is not
purely noise suppression. **HepG2 stays at zero at every resolution** —
aggregation reveals structure where it exists rather than manufacturing it.
Hallmark (45 coarse sets) beats Reactome (874 finer) consistently.

### [Know] Best-behaved candidate D

`D_unexplained_fraction = <h1−A, h2−A> / <h1, h2>`. Unbiased under a stated
noise model; recovers planted `(1−share)²` to ±0.05; correctly **exceeds 1** for
predictions worse than zero (so it must not be clipped); immune to the depth
bias that disqualifies a raw residual norm (which grows monotonically with noise
at fixed latent signal — tested). **Edge case:** the per-pair denominator is a
noisy ‖L‖² and destabilises near zero, so `pooled_unexplained_fraction` is the
headline statistic and per-pair values need a reliability filter. Observed
pooled: 1.296 / 0.752 / 0.798 / 1.007; fraction worse than zero 40–84 %.
**Not selected.**

### Recommendation: B — pathway-level gamma, with scale calibration and source agreement

Evidence-ranked. Pathway gamma is the only place a *large* amount of previously
invisible structure appeared. Scale calibration is a cheap available fix for the
biggest point-prediction failure and should be folded in regardless.
Transferability is real but modest after adjustment — better as a confidence
output than a primary target. **Template recovery from basal is ruled out.**
Exact gene-level gamma is now clearly dominated.

Lowest-capacity model justified: scale-calibrated conserved transfer (one
scalar) + source-agreement confidence. Two scalars and one descriptive feature.

### Commands run

```
uv run python scripts/run_transferability_foundations.py
uv run python scripts/plot_transferability_foundations.py
uv run pytest -v ; uv run ruff check . ; uv run ruff format --check . ; uv build
```

### Next step

Stop for review before building anything.

---

## 2026-09-20 (later) — Pathway gamma representation falsification

Froze foundations first (`foundations_v1_freeze.txt`, 12 digests). All six freeze
records re-verify. Report: `reports/pathway_gamma_falsification_v1.md`. Runtime
5.0 min. **No model trained.**

### Design

Every representation — Hallmark, Reactome, matched random gene sets, random
projections — applied as the same linear map `R @ W.T`, so output dimensionality
and arithmetic are exactly matched and only *which genes are grouped* differs.
Primary null is **gene-label permutation**, which preserves set sizes, **all
pairwise set-set overlaps** (`(MP)(MP)^T = MM^T`) and the gene-degree multiset.
100 replicates per null type, constructions fixed in advance.

### [Know] Pathway biology beats matched random aggregation — 3 of 4 contexts

Reliability-normalised gamma recovery:

| held out | Hallmark | permuted null | gaussian null | genes | p | z |
|---|---:|---|---|---:|---:|---:|
| K562 | **0.652** | 0.305±0.057 | 0.300±0.040 | 0.305 | 0.010 | +6.1 |
| RPE1 | **0.194** | 0.015±0.040 | 0.015±0.027 | 0.002 | 0.010 | +4.5 |
| HepG2 | 0.032 | 0.065±0.037 | 0.065±0.029 | 0.062 | 0.842 | **−0.9** |
| Jurkat | **0.599** | 0.378±0.032 | 0.363±0.024 | 0.372 | 0.010 | +6.8 |

p = 0.010 is the 100-replicate floor: no null replicate reached the observed
value in those three folds.

**The decisive panel:** random 45-dim aggregation reproduces **gene-level
performance almost exactly** after reliability normalisation (K562 0.305 vs
0.305; Jurkat 0.378 vs 0.372). Raw `r_gamma` does rise under random aggregation
(0.185 → 0.194 for K562) purely because averaging raises reliability — and the
normalisation removes exactly that. **So the whole benefit of dimensionality
reduction per se is a reliability artefact; what remains is biology.**

All three null types are statistically indistinguishable, so aggregation
geometry (overlaps, gene degrees, sparsity) contributes nothing.

**HepG2 is a clean negative**: Hallmark (0.032) sits *below* its null (0.065).
Pathway aggregation does not manufacture signal where none exists.

Reactome (874 sets, independent ontology) reproduces the pattern fold-for-fold
against its own permuted null: K562 z=+7.9, RPE1 z=+7.9, Jurkat z=+5.1, HepG2
z=+1.5 (n.s.). Hallmark beats Reactome in the two strong folds.

### [Know] Dependence is pair-specific and asymmetric

Hallmark, every two-source subset (separates "lost a partner" from "fewer
sources"):

| target | all 3 | drop K562 | drop RPE1 | drop HepG2 | drop Jurkat |
|---|---:|---:|---:|---:|---:|
| K562 | +0.570 | — | +0.341 | +0.567 | **+0.337** |
| RPE1 | +0.182 | +0.171 | — | +0.011 | +0.182 |
| HepG2 | +0.026 | +0.003 | +0.116 | — | −0.078 |
| Jurkat | +0.510 | **−0.195** | +0.424 | +0.533 | — |

- **Jurkat is wholly K562-dependent**: −0.195 without it, but +0.533 when HepG2
  is dropped instead. Spread across dropped sources = 0.73.
- **K562 survives without Jurkat** (+0.337, −41 %); losing RPE1 costs the same.
  Its signal is not specifically Jurkat-dependent.
- RPE1's modest signal depends on **HepG2** — a third, different pairing.
- Two-source reduction is not the cause: RPE1 is unchanged at two sources.

### [Disproved] The dataset-ancestry confound I have been flagging since 2026-09-19

K562/RPE1 are Replogle; HepG2/Jurkat are Nadig. The two **same-dataset** pairs
rank 4th and 5th of six on gamma excess over null (−0.074, −0.077); mean excess
is **−0.076 same-dataset vs +0.034 cross-dataset**. Both informative pairings
(K562–Jurkat, RPE1–HepG2) **cross** the dataset boundary. Dataset ancestry does
not explain gamma sharing and is mildly anti-correlated with it. Weak (n=6, two
same-dataset pairs) but it points firmly away from the confound.

**The lineage confound for K562–Jurkat (both suspension leukaemia) is untouched
and remains unresolvable with four contexts.**

### Decision: pathway-level modelling is JUSTIFIED, with scope stated

Criteria 1, 2, 4 pass cleanly; criterion 3 passes in quantified partial form.
Honest framing is **not** "pathway gamma is predictable" but "pathway-level
gamma is recoverable for *some* target contexts, and which ones is itself
context-dependent and not predictable from four contexts."

Scope for any model built next: evaluate **per held-out context, never pooled**;
report the matched-random null alongside every number; treat **HepG2 as a
known-negative control** (success there is suspect); do **not** claim
independence between K562 and Jurkat results; carry scale calibration and
source-agreement confidence regardless.

### Commands run

```
uv run python scripts/run_pathway_falsification.py
uv run python scripts/plot_pathway_falsification.py
uv run pytest -v   # 233 passed
uv run ruff check . ; uv run ruff format --check . ; uv build
```

### Next step

Stop for review.

---

## 2026-09-20 (later) — Pathway residual model v1: NEGATIVE result

Froze all discovery artifacts first (`discovery_phase_freeze.txt`, 26 digests).
All eight freeze records verify **before and after** modelling. Modelling code
went into a new namespace `virtual_cell.modelling`; no frozen module changed.
Report: `reports/pathway_residual_model_v1.md`.

### Headline

**The learned pathway correction never improves the outer target.** Selected by
inner pseudo-LOCO: K562 **M0/lam=0**, RPE1 M1/lam=0.75, HepG2 M2/lam=0.75,
Jurkat **M0/lam=0**. Outer deltas: K562 0.0000, RPE1 **−0.0153**, HepG2
**−0.0034**, Jurkat 0.0000. M0 is at least as good as M1/M2 on Pearson in every
fold. M3 never evaluated (precondition not met).

### [Know] It is NOT a nested-LOCO selection failure — oracle sweep settles it

Oracle lambda sweep (evaluation-only): max attainable gain is **+0.0033**
(K562), **+0.0024** (HepG2), **0.0000** (RPE1, Jurkat). lam=0 is essentially
optimal everywhere. The obvious alternative explanation — "inner folds can't see
the K562–Jurkat pairing because the partner is the outer target" — is excluded.

### [Know] The mechanism: beta contamination cancels a real gamma gain

Residual algebra, asserted as a test:
`R = Y − B = (4/3)alpha + (1−s)beta + (1+s/3)gamma`. Centring removes alpha
exactly but **leaves (1−s)beta**, which at s≈0.55 is nearly half the conserved
effect. Forced-M2 diagnostic:

| held out | r(R_hat, gamma) | r(R_hat, beta-like) | deterministic r(gamma-hat, gamma) |
|---|---:|---:|---:|
| K562 | **+0.469** | **−0.424** | +0.570 |
| Jurkat | **+0.428** | **−0.491** | +0.510 |
| RPE1 | +0.132 | +0.838 | +0.182 |
| HepG2 | −0.054 | +0.767 | +0.026 |

**The model DOES find gamma in K562/Jurkat (r≈0.43–0.47, approaching the
deterministic 0.51–0.57) — but simultaneously predicts the beta component with
the WRONG SIGN, and the two cancel.** In RPE1/HepG2 the correction is almost
pure beta re-prediction. Third link in a chain that must never be collapsed:
*gamma exists and is reproducible* ≠ *gamma is recoverable zero-shot* ≠
**predicting gamma improves the response prediction**. The first two hold; the
third is falsified for this construction.

### [Know] Other results

- **Hallmark vs matched-random modelling: indistinguishable** (p=0.48–1.00,
  z=−1.0 to +0.2) — because neither produces a gain. Does NOT contradict the
  falsification battery, which compared gamma *recovery*, not response
  *prediction*.
- **Reactome reproduces the pattern**: one small positive (K562 +0.012), two
  negatives, one abstention.
- **Features**: two dominate by an order of magnitude — `source_mean` (0.082)
  and `weighted_source` (0.072), i.e. the raw conserved response. Every
  biological/target feature is ~0: source_agreement 0.0021 (sign-inconsistent),
  basal_pathway_deviation 0.0019 (sign-inconsistent), target_gene_basal 0.0007.
- **[Criticism of our own model] It did NOT shrink in the known-negative.**
  HepG2 got lam=0.75 and degraded (−0.018 energy); its correction is pure beta
  (r=+0.767) with no gamma (−0.054). Shrinkage worked in K562/Jurkat, not HepG2.
- **Source-agreement confidence works**: Spearman(agreement, corrected r) =
  +0.602/+0.534/+0.529/+0.572, strictly monotone by quartile in all four
  contexts. **The only component of the system that works as intended.**

### Limitation inherent to n=4

Inner pseudo-targets are built from **two** sources while the outer prediction
uses **three**, so `source_mean` has different noise/attenuation at fit vs
application time. Unavoidable here; partially absorbed by per-fold
standardisation; genuinely limits what any nested-LOCO model can learn.

### Next step (user has NOT approved)

**(a) One narrow pre-specified retest**: define the learning target against the
*unshrunk* transfer (s=1), where `centred R = (4/3)gamma` exactly and the beta
contamination vanishes by construction, then add the gamma correction to the
scale-calibrated baseline. One-line change; directly removes the identified
mechanism. **If it does not produce an outer gain in K562 or Jurkat, stop
pathway modelling.** Do not iterate beyond one attempt.

**(b) Otherwise pivot to transferability / D**, the only demonstrably working
component.

### Commands run

```
uv run python scripts/run_pathway_residual_model.py
uv run python scripts/plot_pathway_residual_model.py
uv run pytest -v   # 264 passed
uv run ruff check . ; uv run ruff format --check . ; uv build
```

---

## 2026-09-20 (final) — Pathway residual model v2 (clean gamma): NEGATIVE; pathway modelling TERMINATED

Froze v1 permanently first (`pathway_residual_v1_freeze.txt`, 15 digests, with
the retained conclusion written into the freeze file itself). All nine freeze
records verify before and after v2. Report:
`reports/pathway_residual_model_v2_clean_gamma.md`.

**Exactly one thing changed: the training target.** Enforced by
`test_features_are_identical_to_v1` (v2 design matrix must be byte-identical to
v1's) paired with `test_v2_target_differs_from_v1_target` so the claim cannot be
vacuous.

### [Know] The beta fix worked — gamma prediction improved a lot

`centre_p(Y[c] - A) = gamma[c] - mean_S gamma`; beta cancels **exactly, before
centring, for any source-set size**. With 3 sources it is `(4/3)gamma`
(verified to 1e-10 on planted data); with 2 sources (inner folds) it is
`gamma[c] - (gamma[a]+gamma[b])/2` — beta-free but a different combination.

| held out | r(R_hat, gamma) v1 | **v2** | ||R_hat||/||R_true|| v2 |
|---|---:|---:|---:|
| K562 | +0.469 | **+0.778** | 0.98 |
| Jurkat | +0.428 | **+0.623** | 0.77 |
| RPE1 | +0.132 | +0.182 | 0.22 |
| HepG2 | −0.054 | +0.125 | 0.37 |

### [Know] And it STILL did not improve response prediction — in any context

| held out | sel lam | lam_theory | Δr | bootstrap CI |
|---|---:|---:|---:|---|
| K562 | 0.25 | 0.90 | −0.0036 | [−0.0134, +0.0069] |
| RPE1 | 0.25 | 0.89 | **−0.0173** | **[−0.0243, −0.0107] excludes 0** |
| HepG2 | 0.50 | 0.88 | −0.0042 | [−0.0195, +0.0039] |
| Jurkat | 0.25 | 0.90 | −0.0060 | [−0.0133, +0.0051] |

**The decisive diagnostic: at lam_theory (0.88–0.90, the mathematically correct
amount) every context is substantially WORSE** — K562 −0.0175, RPE1 −0.1003,
HepG2 −0.0260, Jurkat −0.0467. If R_hat were an accurate `(4/3)gamma` estimate,
lam_theory would be optimal by construction. That it is strongly harmful means
the correction is **directionally right but not accurate enough per element**:
at r=0.78, ~40 % of its variance is error, injected into a response whose gamma
share is ~21 %. Inner selection correctly shrank lam to 0.25–0.5 but no positive
lam helps.

Hallmark vs matched-random: p=0.67–0.95, all z ≤ 0. Reactome reproduces the
pattern fold-for-fold (RPE1 again significantly harmful; Jurkat abstains at
lam=0).

### TERMINATED — predeclared stopping rule fired on all three conditions

(1) K562/Jurkat show no usable gain; (2) three of four CIs include zero and the
fourth is negative; (3) matched-random performs equivalently. Any one sufficed.

**No pathway model v3.** Revisiting would need new external evidence — more
contexts or a materially more accurate gamma estimator — not another iteration
on these four cell lines.

### The chain, now complete

1. gamma exists and is reproducible — **TRUE**
2. gamma is recoverable zero-shot at pathway resolution, genuinely from biology
   not aggregation — **TRUE** (p=0.010 vs 100 geometry-preserving nulls)
3. a learned gamma correction improves the actual response prediction —
   **FALSE**, tested twice, mechanism identified and removed, still false.

Step 3 fails quantitatively, not conceptually: **a correction must be far more
accurate than "substantially correlated" before it pays for itself.**

### Honest counter-evidence recorded

One inner-selection miss: for RPE1 the capacity ladder shows M2 at lam=0.25
would have scored 0.7206 vs baseline 0.7058 (**+0.0148**), but inner selection
chose M1 (0.6885). So in one of four folds a better option existed and was not
selected. It does not change the verdict — a single unreplicated +0.015 that was
not selectable is exactly what the predeclared rule exists to discount — but it
is the strongest counter-evidence available and is in the report.

### Next primary direction (stated explicitly as required)

**Transferability / confidence, NOT another gamma model.** Scale-calibrated
conserved transfer is the best available point predictor and was never beaten.
Source agreement predicts transfer success at Spearman +0.50 to +0.61 within
every outer context, survives joint confound control (+0.14 to +0.55), is
monotone by quartile, and is computable at inference from sources alone — the
only component that has worked in every experiment.

### Commands run

```
uv run python scripts/run_pathway_residual_model_v2.py
uv run pytest -v   # 281 passed
uv run ruff check . ; uv run ruff format --check . ; uv build
```

---

## 2026-09-20 (final) — Transferability/confidence model v1: THE SIMPLE STATISTIC WINS

Froze pathway modelling permanently first (`pathway_modelling_terminated_freeze.txt`,
12 digests with both retained conclusions written into the freeze file). All ten
freeze records verify before and after. Report:
`reports/transferability_confidence_model_v1.md`. **Point predictor unaltered.**

### [Know] Reliability-aware targets, validated synthetically

`signal_energy = <h1,h2>` (unbiased for ||L||^2), `residual_energy = <h1-B,h2-B>`
(unbiased for ||L-B||^2), `D = residual/signal`. Synthetic checks: naive ||h||^2
inflated >20 % while cross-half is unbiased to 5 %; D recovers planted 0.25 /
0.0625 / 1.0. **D never clipped** — the D>1 regime is real and observed.

**Stability rule derived from simulation BEFORE any outer fold**, then frozen:
using the empirical source half-noise sd (0.0600, 6,640 genes), min signal
energy = **1.9674** (within-tolerance 0.14 → 0.77 → 1.00 across signal
0.16 → 1.13 → 3.48). My first derivation returned the grid's own lower bound
because the grid never reached the unstable regime — re-derived over a wider
range. In practice it excludes almost nothing (stable fraction 0.995–0.999).

### [Know] Raw source agreement predicts transfer quality in every context

Spearman(score, −D): K562 **0.554**, RPE1 **0.786**, HepG2 **0.656**, Jurkat
**0.556**. Risk-coverage is **monotone in all four contexts** on both metrics.
Median D falls: K562 1.231 → 0.728 (**−40.9 %**) at 10 % coverage, HepG2
0.895 → 0.612 (−31.6 %), Jurkat 1.011 → 0.675 (−33.3 %), RPE1 0.805 → 0.664
(−17.5 %). Calibration strictly monotone by quintile in all four.

**K562 and Jurkat sit above D=1 at full coverage** — the transfer is on average
worse than predicting zero there — and selection brings them below 1. That is
the operational value.

**HepG2 works.** First phase in which it behaves normally (Spearman 0.656,
second highest; −23.1 % D at 50 % coverage).

### [Know] Neither M1 nor M2 beats M0 — predeclared rule fires, STOP

M0 >= learned on ranking in **4/4** contexts and on risk-coverage in **16/20**
(context, coverage) cells. Ridge coefficients are sign-inconsistent for 9 of 12
features, with obvious collinearity (`source_agreement` and
`source_min_agreement` get identical coefficients). **Final estimator = raw
source agreement.** One number, no fitting, no hyperparameters, nothing to leak.

### [Know] Trustworthiness and error magnitude are near-OPPOSITE objectives

**Targeting the least-confident perturbations captures ~HALF of random**
(M0 at 20 % budget: 0.10–0.12 vs random 0.20). Structural reason: low source
agreement tracks low *reproducible signal*, and a perturbation with no
reproducible response contributes little absolute error however badly predicted.

Ranking by **expected error magnitude** instead captures ~2x random (20 % budget:
**0.38–0.41**). But that ranking is a terrible trustworthiness score — its
risk-coverage similarity at 10 % coverage is 0.23–0.41, *below* the
all-perturbation baseline. **A practical system needs two scores for two
questions and must never use one for the other's job.**

Also worth remembering: **AURC on absolute residual energy is a degenerate
metric** — `source_magnitude` and even *random* beat source agreement on it,
because selecting tiny-response perturbations trivially minimises absolute
error. Only the normalised metrics (D, reliability-normalised similarity) are
meaningful for trustworthiness.

### Scientific framing implication

The defensible claim is **not** "we predict context-specific responses" but
**"we predict conserved responses, and can say in advance, per perturbation and
per unseen context, how much to trust each prediction."** Directly actionable
for Arc: a confidence score can gate or weight submitted perturbations with no
additional measurement.

### Commands run

```
uv run python scripts/run_transferability_confidence.py
uv run python scripts/plot_transferability_confidence.py
uv run pytest -v   # 314 passed
uv run ruff check . ; uv run ruff format --check . ; uv build
```

### Next step

Stop for review.

---

## 2026-09-20 — Arc Challenge Bridge v1

Turning the frozen findings into an Arc-compatible predictor over the official
18,533-gene raw-count output space. **No model fitted, no submission made, no
prior conclusion modified.** Gamma, pathway-residual and confidence modelling
all remain closed.

Reports: `reports/arc2026_submission_requirements.md`, `reports/arc_bridge_v1.md`
New package: `virtual_cell.arc` (`panel`, `generate`, `metrics`)

### [Blocker] Zero of the 300 Arc targets are perturbed in ANY research context

Across the four contexts the whole programme was built on — `replogle22k562`,
`replogle22rpe1`, `nadig25hepg2`, `nadig25jurkat` — **0/300 Arc targets appear
as perturbations**, though all 300 appear as measured output genes. Replogle and
Nadig screen **essential-gene** libraries; Arc's panel is largely non-essential
regulatory genes. The libraries are near-disjoint by construction.

Checked, not assumed: source-to-source perturbation intersections are large and
correct, NTC counts match published values exactly, and spot checks confirm
"perturbed in: NONE / measured in: all four".

Widening to all seven public scPertEval datasets: **86/300** targets have any
perturbation data (`kaden25rpe1` 80, `arch1` 13), and only **7** have it in ≥2
contexts — the minimum for separating `beta_p` from `gamma[c,p]`.

**This is a different generalization axis.** We characterised unseen *context*.
For 71% of Arc's panel the problem is an unseen *perturbation*, about which this
project has produced no evidence either way.

### [Know] The score's zero is the MEAN RESPONSE, not the control

From `cell-eval2` `configs/vcc2026.yaml` and
`docs/vcc2026_metrics/vcc2026-metrics-brief.md` (rule_version 3): six equally
weighted members, each rescaled `s = (u − b)/(r − b)` with `b` the context's
mean perturbation response and `r` a replicate of the experiment.

**A control-emitting submission scores −0.311**, not 0 — almost entirely via
`de_wilcoxon_direction_fidelity_yield_raw`, where calling nothing scores 0
against a baseline of 0.51. Predicting no effect is a third of a point *worse*
than predicting the panel mean.

The useful corollary: averaging `delta[c,p] = mu + alpha_c + beta_p + gamma[c,p]`
over `p` leaves `mu + alpha_c` by the centring constraints, so **the scale's zero
is exactly the context main effect** — a quantity the foundations phase already
studied. Estimating it for A/B/C is the cheapest real gain available.

### [Know] The output-space gap is not a gap

166 Arc genes are measured in no public dataset (olfactory receptors, defensins,
KIRs, testis/keratin genes). **Zero of them clear the metric's 5 CPM control
gate in any context**, so four of six members never test them, and they hold
**0.0000%** of the control profile's squared norm for the other two. Restricted
to the 11,957 genes that are tested, public coverage is **100.00%**.

### [Build] `virtual_cell.arc`

- `panel` — categorises every panel gene (`predicted` / `measured_no_response` /
  `unmeasured`); `fill` is an explicit argument so zero-filling is a recorded
  decision, never a default.
- `generate` — `G0` control resampling, `G1` control transport, `G2`
  negative-binomial count model. All validated against 3,000 real control cells
  per context: integer, non-negative, depth and density matching the controls.
- `metrics` — all six scored members reimplemented, validated against the
  values the specification states analytically.

**Defect found and fixed:** `G1` re-drawing from a cell's own empirical
composition lost **17% of detected genes**, which would have biased every
Wilcoxon member. Blending with the pooled control composition (`smoothing=0.5`)
restores density to 0.320 against the controls' 0.323.

### [Warning] Arc's contexts are FURTHER APART than any pair we measured

Basal `log1p(CPM)` correlations: A–B 0.780, A–C 0.749, B–C 0.842. The four
public contexts sit at 0.89–0.93. **Every transfer number in this project was
obtained on easier context pairs than Arc will score.** Our results are an
optimistic bound, not a neutral estimate. (Distance only; no identity inference.)

### Commands run

```
vcc --version ; vcc prep --help ; vcc sample --help
uv run python scripts/run_arc_bridge.py
uv run pytest -q   # 356 passed
uv run ruff check . ; uv run ruff format --check . ; uv build
```

### Next step

Stop for review. Two decisions are the user's:
1. Approve the 10 GB `arch1` + `kaden25rpe1` download (the only public data with
   any Arc target coverage, and an Arc-like held-out-context benchmark at 97%
   gene overlap).
2. Approve estimating `mu + alpha_c` for A/B/C from basal controls, which moves
   the floor from −0.31 to 0 across all 300 targets using permitted inputs only.

---

## 2026-09-20 — Unseen-perturbation generalization v1

The Arc bridge established that the central problem is **unseen context +
unseen perturbation**, not unseen context with seen perturbations. This phase
built that benchmark and ran it.

Report: `reports/unseen_perturbation_generalization_v1.md`
New: `virtual_cell.priors` (catalogue, features),
`virtual_cell.modelling.unseen_perturbation`,
`virtual_cell.modelling.context_main_effect`
Acquired: `arch1`, `kaden25rpe1` (approved), `wessels23` (0.2 GB, so the
coverage audit could be entirely local). All md5-verified against upstream.

### [Know] Unseen beta_p IS predictable — internally

STRING-network k-neighbour transfer, with no measurement of the knockdown
anywhere: **P1 r=+0.509, 24.7% of response energy removed.** Both axes held
out (P2): r=+0.235, 2.9% removed.

Prior ranking, consistent across regimes: **STRING > DepMap > Reactome > basal
> Hallmark.** Network/DepMap/pathway priors beat target-gene basal expression
by ~20x in energy removed. Combining families adds nothing over STRING alone.

### [Negative] It does NOT survive a real context shift

On `arch1` — 100 perturbations absent from every source, a genuinely external
context, estimator and hyperparameters frozen beforehand — the same predictor
scores **r=+0.057 and leaves MORE error than predicting zero** (1.704 vs
1.000). No estimator beats U0.

The isolating contrast, same context and gene axis: the 17 `arch1`
perturbations that ARE measured in the sources still transfer at r=+0.302,
unexplained 0.966 under the frozen scale calibration. **Direct measurement of a
perturbation survives the context shift; inferring it from gene annotations
does not.**

The internal benchmark was measuring transfer between four screens sharing a
study design, protocol and gene panel. Arc's contexts are further apart from
each other (basal r 0.75–0.84) than ours are (0.89–0.93).

### [Negative] The feasible context main effect also fails externally

`m_c = mu + alpha_c` is the VCC score's zero. Internally, E2 (basal-weighted,
one scalar fitted inner-LOO) beats the zero response in all four folds — but
only because of the scalar: E0/E1 recover direction (r 0.68–0.82) and get
magnitude badly wrong (K562 overshoots 2.1x, RPE1 undershoots to 0.36x).

On `arch1`: **with a perfect oracle scalar, 89.0% is still unexplained**
(direction r=0.33), against 32–51% internally. `arch1`'s responses are 3–5x
smaller in norm than the four screens'. **Target controls do not automatically
deliver score 0.** This is a limit, not a calibration bug.

### [Know] The Tier-0 confidence statistic is NEIGHBOUR AGREEMENT

Source agreement is undefined without source responses. Two candidates tested:

- **Support distance fails** — Spearman vs accuracy: basal **−0.25**
  (backwards), STRING −0.06, combined +0.02. Hallmark ranks highest and is the
  worst predictor, so the ranking is anti-correlated with usefulness.
- **Neighbour agreement works** — do the genes standing in for the unseen one
  agree with *each other*? DepMap +0.50, combined +0.47, Reactome +0.43,
  STRING +0.26. Holds externally: **+0.314 on arch1**.

### [Know] Ceilings that make the numbers interpretable

In P2: `oracle_beta` leaves 0.587 unexplained (the rest is gamma — permanently
closed), scale-calibrated Tier-2 transfer 0.903, best unseen prediction 0.971,
zero 1.000. So the entire budget for any beta model is 41.3%; Tier-2 captures
9.7%, Tier-0 captures 2.9%.

Also: **raw conserved transfer scores 1.043 — worse than zero.** The frozen
scale calibration is not a refinement, it is what makes transfer beat doing
nothing.

### [Decision] No nonlinear model. Tier 0 gets no prior-based beta.

Ridge already reaches the k-NN ceiling internally and collapses identically
externally, so the failure is not capacity. 30.8% of Tier-0 Arc targets have no
STRING representation at all. Adding a prior-based beta to Tier 0 would be
optimising against the benchmark that `arch1` just proved misleading.

Support tiers frozen in `data/splits/arc_target_support_v1.csv`:
**TIER 2 = 7, TIER 1 = 79, TIER 0 = 214.**

### Leakage protocol

Both axes held out in the hardest regime; projections fitted on training
perturbations only. Leakage is tested, not argued: replacing the outer
context's responses with noise, replacing held-out perturbations' responses in
every context with noise, and both at once, each leave predictions
**bit-identical** for every estimator — plus a guard test confirming that
corrupting the *readable* block does change them.

### Commands run

```
uv run python scripts/audit_public_datasets.py
uv run python scripts/run_unseen_perturbation.py
uv run python scripts/run_arch1_external.py
uv run pytest -q   # 405 passed
uv run ruff check . ; uv run ruff format --check . ; uv build
```

All 18 freeze manifests verified: 199 files, 0 failed.

### Next step

Stop for review. The main caveat on the central negative claim is that `arch1`
is a **single** external observation; a second external context would confirm
whether the collapse is general or specific to hESC biology and that study's
CRISPRi efficiency.

---

## 2026-09-21 — External unseen-perturbation validation v1

Diagnostic phase: **why** did the internal unseen-perturbation result fail on
arch1? No model built, nothing tuned.

Report: `reports/external_unseen_perturbation_validation_v1.md`
New: `virtual_cell.modelling.external_benchmark` (transcription of the frozen
arch1 protocol, plus `target_reliability`)
Audit record: `data/provenance/external_candidates.md`

### [Method] The protocol is provably the frozen one

`run_arch1_external.py` is frozen, so the protocol was transcribed and the
transcription **checked**: the script re-scores arch1 and refuses to touch
Kaden unless every number matches `arch1_external.csv`. Max difference
**9.0e-17**.

### [Negative, and uninformative] Kaden cannot answer the question

`kaden25rpe1` was meant to separate study shift from context shift — an
independent screen in RPE1, where Replogle RPE1 is already a source.

Unseen prediction fails (knn r=+0.016, ridge +0.018; both worse than zero on
energy). **But the positive control fails too**: directly measured transfer of
the same 113 perturbations in the same cell line gives **r=-0.033** (arch1's
equivalent: +0.302).

Checked before interpreting: Kaden's basal control profile correlates with
replogle22rpe1 at **r=0.974**, the highest pair in the project — cell line and
gene axis are right. Yet for the same 113 perturbations Kaden vs each source is
-0.013 to -0.016, while replogle22rpe1 vs the other three *cell lines* is +0.157
to +0.205. Two screens of the same perturbations in the same cell line agreeing
*less* than screens in different cell lines is not a plausible biological
result.

### [Know] The cause: target-side reliability

Split-half Spearman-Brown of each target context's own responses:

| context | subset | cells/pert | SB | ceiling |
|---|---|---|---|---|
| arch1 | measured (17) | 790 | **0.952** | 0.976 |
| arch1 | unseen (100) | 1,182 | **0.838** | 0.916 |
| replogle22rpe1 | — | — | 0.570 | 0.755 |
| kaden25rpe1 | measured (113) | 379 | **0.128** | 0.358 |
| kaden25rpe1 | unseen TF (300) | 400 | **0.115** | 0.339 |

Kaden is the least reproducible context in the study — 7x below arch1 — with
the smallest effects (‖delta‖ ~1.9 vs 3.6-6.1 in sources) despite healthy cell
counts. So this is weak knockdown or weak transcriptional consequence, not thin
sampling. Disattenuation changes nothing.

**Verdict: Interpretation C.** Kaden is uninformative; it supports neither A nor
B. The study-vs-context question **remains open**. The arch1 conclusion is
unchanged and strengthened — its reliability is the highest measured, so its
failure was never attenuation-limited.

### [New binding criterion] Qualify a benchmark on reliability BEFORE using it

A negative result on a target whose responses do not reproduce is
uninterpretable: the ceiling on any correlation is sqrt(rho). This criterion
did not exist before and now governs candidate selection.

### [Audit] Second external context

**Accession correction: GSE212396 is the PILOT** (50/200-gene libraries). The
~2,285-gene TeloHAEC screen is **GSE210681**.

| | TeloHAEC GSE210681 | Feng targeted |
|---|---|---|
| context | endothelial (very distinct) | iPSC (pluripotent, like arch1) |
| Arc targets | **66** | 5 |
| our gene axis | **6,613/6,640** | 5,213/6,640 |
| unseen perturbations | 1,925 nominal | 141 |
| measured-transfer controls | 250 | **170** |
| **targets with >=10 sig genes** | **93/2,345 (4.0%)**, median 0 | **444/444 (100%)**, median 189 |
| usable unseen | **48** | **141** |

TeloHAEC wins every structural criterion and fails the one that matters. Its
weakness is corroborated by the authors' own program-level MAST (12.9% with any
significant program) and is consistent with their design — they aggregate into
60 NMF programs because individual CAD-locus perturbations are weak.

**Selected: Feng.** Limitations stated up front — only 5 Arc targets, a
pluripotent context not maximally distant from arch1, 5,213 shared genes, and
median 74 cells/target/line so only the pooled estimate is usable. Its
distinctive value is 20 cell lines, which tests degradation across contexts
with study and protocol held fixed — the variable Kaden failed to isolate.

### [Unchanged] Arc Tier-0 policy

Zero perturbation-specific beta for Tier 0. The policy changes only on new
external evidence; this phase produced none in either direction.

### Commands run

```
uv run python scripts/run_external_validation.py
uv run pytest -q   # 411 passed
uv run ruff check . ; uv run ruff format --check . ; uv build
```

All 19 freeze manifests verified.

### Next step

Stop for review. One decision is the user's: approve the **1.92 GB** download of
Feng `TargetedScreen_LFC_byGene-perLine.tsv.gz`. Nothing further was downloaded.
