# Zero-Shot Virtual Cell Modeling Across Unseen Cellular Contexts

Research code for predicting transcriptional responses to CRISPRi genetic
perturbations in cellular contexts that were never seen perturbed.

The project serves two purposes at once:

1. an entry to the **2026 Arc Institute Virtual Cell Challenge**, and
2. an **independent research project** on zero-shot cross-context perturbation
   prediction.

The full specification, phased plan, and scientific background live in
[`plans.MD`](plans.MD). Literature notes are in
[`reports/literature_notes.md`](reports/literature_notes.md) and a running
research log in [`reports/research_log.md`](reports/research_log.md).

## Research question

Molina and Zhang (bioRxiv, July 2026) decompose a pseudobulk perturbation
response into a global term, a cell-line term, a **conserved perturbation
effect** shared across contexts, and a **context x perturbation interaction**.
On four CRISPRi cell lines the interaction holds about a quarter of the
reproducible variance, the conserved effect transfers to unseen lines, and
the interaction could not be predicted zero-shot from basal expression, gene
priors, or source-context responses by any model they tested.

This project therefore asks three questions, in order:

- **A.** How much context-specific response is identifiable at all under truly
  zero-shot context shift?
- **B.** Can richer context-conditioned biological priors recover any
  predictable portion of the interaction beyond existing approaches? This is
  treated as a high-risk hypothesis to test, not an expected result.
- **C.** When the interaction cannot be predicted, can we identify
  perturbation/context pairs whose conserved effect is transferable enough to
  act on, and flag the rest for experiment?

A related direction is predicting whether a perturbation is likely conserved
or strongly context-dependent in a new context, which is useful even if the
interaction itself stays unpredictable. The novelty claim is provisional and
will not be frozen until the decomposition is reproduced on our own splits.

## The Arc 2026 task in one paragraph

Arc performed CRISPRi Perturb-seq in six undisclosed cell lines from
different tissues of origin. Participants
receive only non-targeting control cells for each context (about 18,400 per
context) and a list of about 300 target genes. They must submit predicted
**raw single-cell count profiles**: exactly 400 cells per perturbation, across
all 18,533 genes, for all three contexts in the current phase, with no control
cells included. Three contexts (A/B/C) are used for validation with a live
leaderboard; three different contexts (D/E/F), released October 22, 2026, are
used for the final ranking. Final submissions are due November 5, 2026, 23:59
UTC. Scoring uses six metrics normalised so that Arc's official mean-response
baseline scores 0 and an experimental replicate anchor scores 1. A
perturbation-specific conserved-effect model is not that baseline and can
score above 0. Sources are cited in
`reports/literature_notes.md`.

## Current status (as of 2026-09-18)

Phase 1 (understand the problem and set up infrastructure). No models have been
trained. The official Arc 2026 **validation control bundle has been downloaded
and audited**; no public perturbation datasets have been downloaded yet.

Implemented:

- `virtual_cell.data.io`: robust `.h5ad` loader with integrity checks
  (non-empty, unique genes and cells, finite non-negative integer counts,
  exactly one context label per file, exact label matching).
- `virtual_cell.data.summary`: per-context statistics (cells, genes, library
  size mean/median/min/max, sparsity, genes detected per cell).
- `virtual_cell.data.arc2026`: the official Arc 2026 control bundle as a module
  — manifest/`gene_names.csv`/`pert_counts.csv` loading, memory-safe streaming
  over the CSR count matrices, a single-pass per-context audit, cross-context
  invariant checks, panel composition, and checksum helpers.
- `virtual_cell.data.synthetic`: deterministic synthetic control contexts for
  pipeline development (Poisson noise, no biology).
- `virtual_cell.preprocessing.pseudobulk`: library-size normalisation, mean
  expression profiles, shared-gene alignment, cross-context basal comparison.
- `scripts/audit_arc2026_controls.py`: reproducible read-only audit of the
  official controls; writes tables and figures to
  `outputs/arc2026_controls_audit/`.
- `scripts/make_synthetic_controls.py`: writes synthetic contexts A/B/C.
- `scripts/explore_synthetic_contexts.py`: prints AnnData structure, summary
  table, basal-mean comparison, and saves a figure to `outputs/exploration/`.
- `tests/`: 60 tests covering the data assumptions above, 29 of them pinning
  invariants of the official bundle (skipped when the git-ignored data are
  absent).

### Official validation controls, audited 2026-09-18

All 44 invariants passed. Full report:
[`reports/arc2026_controls_audit.md`](reports/arc2026_controls_audit.md).

| | A | B | C |
|---|---|---|---|
| shape (cells x genes) | 18,400 x 18,533 | 18,400 x 18,533 | 18,400 x 18,533 |
| sparsity | 0.678 | 0.702 | 0.683 |
| library size median | 20,109 | 19,946 | 20,034 |
| genes detected median | 6,147 | 5,756 | 6,006 |

Control cells only (`target_gene == 'non-targeting'`), 46 shared non-targeting
guides x 400 cells per context, identical gene order across contexts matching
`gene_names.csv` exactly, raw integer counts stored as float32 CSR. Basal
pseudobulk Pearson: A-B 0.688, A-C 0.602, B-C 0.732 — the three contexts are
far apart at baseline and separate completely under PCA.

Two findings that constrain later work: the 18,533-gene panel excludes all
ribosomal protein genes and mitochondrial rRNA, so absolute expression is not
comparable to unfiltered public data; and context B has a low-depth tail
(2.5% of cells under 2,000 UMIs) that A and C do not, so any per-cell QC
threshold will hit B alone.

Not yet implemented: public perturbation dataset download (Replogle 2022,
Nadig 2025), gene intersection across datasets, differential expression,
leave-one-context-out splits, baselines, the response decomposition, models,
and the `.vcc` submission pipeline. See the research log for the ordered next
steps.

## Setup

Requires Python 3.11 and [`uv`](https://docs.astral.sh/uv/).

```bash
uv sync                                   # create .venv and install everything
uv run python scripts/audit_arc2026_controls.py    # needs data/raw/arc2026/controls/
uv run python scripts/make_synthetic_controls.py
uv run python scripts/explore_synthetic_contexts.py
uv run pytest
uv run ruff check . && uv run ruff format --check .
```

Arc credentials are not needed for anything in the repository today. When they
are, the official tooling is the `vcc-cli` package (command `vcc`) and the
`cell-eval2` scorer; API keys must never be committed.

## Repository layout

```
plans.MD                 project specification and phased plan
README.md
pyproject.toml           uv project; packages live under src/
scripts/                 reproducible entry points (data generation, exploration)
src/virtual_cell/        research package
  data/                  io, summary statistics, synthetic data
  preprocessing/         normalisation and pseudobulk
  models/ evaluation/ visualization/   placeholders for later phases
tests/                   pytest suite for data assumptions
reports/                 literature notes, research log, data audits
data/raw, data/processed, data/external   git-ignored datasets
data/provenance/         source, checksums and manifests for downloaded data
outputs/                 git-ignored exploration outputs and figures
```

## Research principles

- **Scientific validity over leaderboard rank.** The research track exists
  independently of Arc placement.
- **No leakage.** A held-out context contributes only its control cells; its
  perturbation responses are never seen during training or model selection.
- **Frozen evaluation.** Test contexts and the public leave-one-context-out
  benchmark are fixed before modelling and never used for tuning.
- **Baselines first.** Control-only, mean-perturbation, conserved-effect,
  nearest-context, and linear baselines must be reported before any deep model.
- **Multiple metrics.** Arc's six metrics plus research metrics (conserved and
  interaction recovery, top-DE precision/recall, direction accuracy).
- **Exact labels.** Context labels are preserved byte-for-byte; the loader
  refuses files whose labels do not match expectations.
- **Reproducibility.** Core results run from scripts with fixed seeds; notebooks
  are for exploration only; data and model artefacts stay out of git.
- **Honest reporting.** Negative results, failed hypotheses, and ideas that are
  already in the literature are recorded as such.
