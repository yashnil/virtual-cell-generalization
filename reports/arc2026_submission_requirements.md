# Arc Virtual Cell Challenge 2026 — submission requirements

**Every requirement below was read off the installed tooling and the official
controls bundle during this phase. Nothing here is recalled from memory or
inferred from documentation.** Sources are named per section, and the exact
command output is quoted where it is the authority.

Date: 2026-09-20
Tool: `vcc 0.2.0` (Python 3.12.9, macOS-14.4.1), endpoint
`https://virtualcellchallenge.org`, installed at `/Users/yashnilmohanty/.local/bin/vcc`
Bundle: `data/raw/arc2026/controls/` (audited in `arc2026_controls_audit.md`,
44 invariants passed, hashes frozen in
`data/provenance/scperteval/arc2026_controls_sha256.txt`)

---

## 1. The output object

A submission is a single `.h5ad` packaged into a `.vcc` by `vcc prep`.

| property | required value | authority |
|---|---|---|
| gene dimension | **18,533** | `vcc prep --expected-gene-dim` default `18533` |
| gene identity and order | exactly `gene_names.csv`, in file order | `vcc prep -g/--genes` "Headerless CSV of expected gene symbols, **in order**" |
| contexts | **A, B, C** — all three | `vcc prep --contexts` default `A,B,C` |
| perturbations per context | **300**, exactly the official list | `--verify-targets` default on; `manifest.json` `n_constructs: 300` |
| cells per perturbation | **400** | `vcc prep --cells-per-pert` default `400`; `manifest.json` `cells_per_pert: 400` |
| total cells | **360,000** = 300 x 3 x 400 | product of the above; cap is 400,000 |
| value space | **raw integer counts** | `--require-counts` default: "Require raw integer counts (2026 scores in counts space)" |
| dtype | float32 | `vcc prep -e/--encoding` default `32` |
| perturbation column | `target_gene` | `--pert-col` default; `manifest.json` `pert_col` |
| context column | `context` | `--context-col` default; `manifest.json` `context_col` |
| control cells | **must not be submitted** | `--reject-controls` default: "scoring uses the held-out controls, never yours" |
| max counts per cell | 1,000,000 | `--max-counts-per-cell` default |
| max cells | 400,000 | `--max-cell-dim` default |
| max nonzeros | the scoring hardware's limit | `--max-nnz` default |

`vcc prep` is a "Native reimplementation of `cell-eval prep` — no network calls,
no cell-eval dependency. Catches locally the same errors the server would
reject." So local `--dry-run` is a faithful pre-check, not an approximation.

### The counts requirement is load-bearing

`--require-counts` is the **default** and its help text states the reason
directly: *2026 scores in counts space*. The legacy alternative
(`--no-require-counts`) restores log-normalization, and `--allow-discrete` is
explicitly marked "Legacy". A 2026 submission is therefore a matrix of raw
integer UMI counts, not a log-normalized matrix and **not a rounded
log-normalized matrix**. Any predictor we build must terminate in a count
generator, not in a continuous expression profile that is coerced afterwards.

---

## 2. The controls bundle

`data/raw/arc2026/controls/` — 6 files, verified this phase:

```
context_A.h5ad   224,973,316 B
context_B.h5ad   210,966,111 B
context_C.h5ad   226,056,691 B
gene_names.csv       119,295 B
manifest.json            715 B
pert_counts.csv        1,906 B
```

`manifest.json` verbatim:

```json
{ "season": "2026", "partition": "val", "panel_id": "vcc2026-val-1",
  "contexts": ["A","B","C"], "pert_col": "target_gene",
  "context_col": "context", "control_label": "non-targeting",
  "n_genes": 18533, "n_constructs": 300,
  "per_context": { "A": {"n_perturbations":300,"control_cells":18400,
                         "ground_truth_cells":138400,"n_ntc_ids":46}, ... },
  "cells_per_pert": 400 }
```

B and C carry the same per-context block as A.

- `gene_names.csv`: 18,533 rows, **all unique**, column `gene_name`, HGNC
  symbols, first three `TSPAN6, TNMD, DPM1`, last three `SLC39A4, TKT, PRH1`.
- `pert_counts.csv`: **one column only, `target_gene`**, 300 rows. There is
  **no `n_cells` column and no `context` column**. This matters: `vcc prep`
  says "An n_cells column in --perts wins over it", so with none present the
  400-cells-per-perturbation default governs, and the same 300 targets apply to
  all three contexts. Ground truth is 138,400 cells per context = 300 x 400
  perturbed + 18,400 control.
- `context_{A,B,C}.h5ad`: each **18,400 cells x 18,533 genes**, CSR, float32,
  values confirmed **integer-valued** (spot-checked min 1.0, max 567/663/762,
  `all(v == floor(v))` over the first 200k stored values). `obs` columns are
  `context`, `ntc_id`, `target_gene`. These are non-targeting controls only —
  46 distinct `ntc_id` values per context.

---

## 3. Validated end-to-end with a dummy submission

To confirm the shape requirements empirically rather than by reading help text,
a full-size dummy was generated and validated locally:

```
vcc sample -g gene_names.csv -p pert_counts.csv --full
  -> 360,000 cells x 18,533 genes, CSR float32, integer-valued,
     obs = ['target_gene','context'], density 0.016
vcc prep --dry-run
  -> "targets: verified against the official list"
     "normalization: counts-preserved"
```

This is a pipeline check only. `vcc sample`'s own help says the values are
noise and it is "NOT for scoring well". **No submission was made, and no
prediction was uploaded.**

---

## 4. Basal statistics of the official controls

These are the **only** things the standing scope rule permits us to take from
A/B/C, and they are inference inputs, never training signal. Computed from the
three control files:

| context | cells | median library | q10 | q90 | max | median nnz/cell | density | genes >=1% detected | genes all-zero |
|---|---|---|---|---|---|---|---|---|---|
| A | 18,400 | 20,110 | 9,027 | 34,330 | 52,420 | 6,147 | 0.322 | 11,664 | 2,398 |
| B | 18,400 | 19,950 | 10,550 | 29,600 | 42,670 | 5,756 | 0.298 | 11,545 | 2,853 |
| C | 18,400 | 20,030 | 9,972 | 33,880 | 50,960 | 6,006 | 0.317 | 12,171 | 2,317 |

Basal `log1p(CPM)` profile correlation between contexts:

| | A | B | C |
|---|---|---|---|
| A | 1.0000 | 0.7802 | 0.7493 |
| B | 0.7802 | 1.0000 | 0.8416 |
| C | 0.7493 | 0.8416 | 1.0000 |

**This is a warning, and it is the single most consequential number in this
report after the coverage audit.** The four public contexts used throughout
this project sit at basal correlations of 0.89-0.93 (`zero_shot_recoverability_v1.md`).
A, B and C sit at **0.75-0.84** — they are substantially further apart from each
other than our public contexts are. Every transfer result we have measured was
measured across *closer* context pairs than the ones Arc will score. Our
zero-shot numbers are therefore an **optimistic** bound on Arc performance, not
a neutral one.

Saved: `outputs/arc_bridge_v1/arc_basal_log1p_cpm.npy` (18,533 x 3, float32).

No identity inference was attempted and none of the above depends on knowing
what A, B and C are.

---

## 5. What this implies for a predictor

1. It must emit **raw integer counts**, per cell, for 360,000 cells.
2. It must cover **all 18,533 genes** in the official order — including genes
   with no public measurement, which must be handled by a stated rule, not by
   silent zero-fill.
3. It must cover **all 300 targets in all 3 contexts** — `--verify-targets` is
   on by default and a missing target is a hard reject.
4. It must not include control cells.
5. Per-cell count totals must stay under 1,000,000; the observed control median
   is ~20,000, so a generator calibrated to control library sizes is far inside
   the cap.
