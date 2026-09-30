# Final-round runbook (VCC 2026, D / E / F)

**Purpose:** produce the champion C1 prediction for the final panel **without editing
source code**.

**Champion:** C1a license-clean atlas (`current_champion.md`). No model change is
allowed on release day.

**No submit command appears in this runbook.** Submitting is a human decision, taken
after the package exists.

Pipeline code:

* `src/virtual_cell/competition_v2/{panel,final}.py`
* `scripts/competition_v2/{run_final_panel,package_final_panel,vcc_pack_panel}.py`
* registry: `configs/source_registry.yaml`
* tests: `tests/test_competition_v2_final.py`

---

> **Provenance note.**
> * **Code basis:** the final-round implementation is built on the committed C1 source code
>   (`5a28314`).
> * **Why the submitted package differs:** the validation package that was submitted was
>   emitted from an earlier, uncommitted state of the same code. After emission, four
>   source files (`fusion.py`, `generator.py`, `sources.py`, `evaluation.py`) were edited,
>   apparently a formatting pass. The source statistics are byte-identical.
> * **Size of the difference:** 2,465 of 360,000 cells (0.68 %), each by a one-count
>   redistribution between genes with unchanged library size.
> * **Reproduction:** the current implementation exactly reproduces the committed C1 builder.
>   That builder is the code every C2–C4 public-fold evaluation used; the C1 fold scores
>   match it to ≤ 1e-6.
> * **Frozen separately:** the original validation package remains frozen by SHA-256
>   (`fcc4e2508798805d93b9f1bb3a6957ca7cc8fbaa9f8318296161e4184bc43ea7`).

---

## BEFORE FINAL RELEASE

### Readiness status (2026-09-30)

| check | status | evidence |
|---|---|---|
| **Generic C1 = committed C1 code** | **exact.** The generic pipeline and the unmodified `build_c1_candidate.py --stage emit` produce array-for-array identical predictions on A/B/C: 360,000 cells, nnz 2,437,123,253, 0 differing entries | `outputs/final/regression_generic_vs_committed.json` |
| **Generic C1 vs the *submitted* package** | **not byte-identical, and not caused by the refactor.** 2,465 of 360,000 cells (0.68 %) differ, in 5,214 of 2.44e9 entries. Each difference is ±1 count moved between genes; every cell depth is identical. They sit in 106 of 900 blocks | `outputs/final/regression_generic_vs_frozen_submitted.json` |
| Mock D/E/F panel (111 shuffled targets: 4 unsupported, 25 one-source, 70 two-source, 12 three-plus) | **all 16 local invariants pass**; `vcc prep --dry-run` exit 0 (`verified_targets: true`, nothing dropped); `.vcc` built (not submitted) | `outputs/final/c1_mock_def/{panel_audit.md, provenance.json, package_manifest.json}` |
| Tests | `tests/test_competition_v2_final.py`, 10 tests: D/E/F, arbitrary count, shuffled order, unsupported / single / multi-source, blocked sources, gene order, determinism, exact cells per perturbation, no controls, raw integer counts, and block regression against the committed builder | `uv run pytest -q` |

**Why the submitted package differs.**

* **The cause.** The package was emitted at 15:17 on 2026-09-26. `fusion.py`,
  `generator.py`, `sources.py` and `evaluation.py` were edited at 15:29 the same day and
  committed as `5a28314`. The source statistics are byte-identical to the manifest.
* **What changed.** Only the bulk moment `p_bulk` of the 269 K562-supported targets
  differs, at ≤ 7.6e-5 relative.
* **The pre-edit code was never committed**, so the submitted bytes cannot be
  regenerated.
* **Final-round implication.** The final round will be built with the committed code,
  which is what every public fold (C2–C4) evaluated. This is the same finding as C2
  amendment 2.

**Panel dependence that is part of frozen C1 (not a change).**

* The K562 and CD4 statistics rows are the frozen retained set: the panel's official
  targets, then the H1 2025 training targets. Each source is centred over its rows.
* A new panel therefore **re-prepares K562 and CD4 from raw** with the unchanged C1
  preparation code. `run_final_panel.py` does this automatically: 43 s for the mock,
  and about 3 min for a 300-target panel by the C1 logs. H1 statistics do not depend on
  the panel.

### Validation-only historical builders (do not use for D/E/F)

`scripts/competition_v2/build_c1_candidate.py` (the exact record of submission #2) and
`build_c2_candidate.py` (never run; pinned by the C2 freeze) hard-code A/B/C, 300
targets and 360,000 cells. They are kept unchanged for provenance. The final round uses
`run_final_panel.py` only.

### Keep ready

1. **Raw sources on disk** (do not move them). The paths are in the registry:
   * `data/raw/competition_v2/K562_gwps_raw_singlecell_01.h5ad` (65.8 GB)
   * `GWCD4i.DE_stats.h5ad` (16.8 GB)
   * `pert_counts_Training.csv`
   * `gencode.v47.annotation.gtf.gz`
   * the frozen H1 statistics
2. **Free disk ≥ 30 GB**, for the prediction (~5 GB), compact file, `.vcc` and scratch.
3. **`vcc whoami` shows "ready to submit".** The CLI token lives in the keyring.
4. **Tooling:**
   * the vendored AtlasShift env at `third_party/atlasshift/.venv` (used by
     `compact.py`);
   * `vcc-cli` 0.2.0 (used by `vcc_pack_panel.py`, which calls `vcc.prep.run_prep`).
5. **A green test suite:** `uv run pytest -q`.
6. **Optional dry rehearsal** on the mock (about 5 min):

   ```bash
   uv run python scripts/competition_v2/make_mock_final_panel.py --out outputs/final_mock/controls
   uv run python scripts/competition_v2/run_final_panel.py --controls-dir outputs/final_mock/controls \
       --output-dir outputs/final/rehearsal
   uv run python scripts/competition_v2/package_final_panel.py --output-dir outputs/final/rehearsal \
       --controls-dir outputs/final_mock/controls
   ```

---

## WHEN D / E / F RELEASE

All commands run from the repository root. Measured wall times are on this machine
(10 cores, 64 GB).

### 0. Tooling check (1 min)

```bash
vcc --version                       # expect 0.2.0; if newer, read its changelog before step 5
vcc whoami
uv run pytest -q tests/test_competition_v2_final.py tests/test_competition_v2_c5.py
```

`run_final_panel.py` (step 0) and `package_final_panel.py` both run the compatibility
gate `competition_v2/vcc_compat.py`:

* It **records** the detected vcc version in `panel_audit.md` and the provenance files.
* It **checks** the `vcc prep` options and the `vcc.prep.run_prep` keywords the pipeline
  uses.
* On a mismatch it **stops with an actionable error** (exit code 2) and never adapts
  silently.
* An untested version with every capability present proceeds, with an explicit
  WARNING.

Control files are identified, in order of precedence, by:

1. the file the manifest names;
2. the validation-era `context_<label>.h5ad`;
3. the unique `.h5ad` whose obs `context` is the label.

Ambiguous or unknown layouts stop with a diagnostic (exit code 3). Nothing is guessed.

### 1. Download (network-bound)

```bash
vcc datasets list --json | tee data/provenance/competition_v2/final_datasets_list.json
# pick the final-round controls id from the listing (validation used `controls`)
mkdir -p data/raw/arc2026/final
vcc datasets download <FINAL_CONTROLS_ID> -d data/raw/arc2026/final   # verifies the server checksum
unzip -o data/raw/arc2026/final/<ZIP_NAME> -d data/raw/arc2026/final/controls
ls data/raw/arc2026/final/controls   # expect manifest.json, gene_names.csv, pert_counts.csv, context_<label>.h5ad
```

If the zip nests a folder, point `--controls-dir` at the folder that holds
`manifest.json`.

### 2. Freeze provenance (1 min)

```bash
( cd data/raw/arc2026/final/controls && shasum -a 256 * ) \
    > data/provenance/competition_v2/final_controls_sha256.txt
```

### 3. Coverage audit (< 5 min, including K562 / CD4 re-preparation)

```bash
uv run python scripts/competition_v2/run_final_panel.py \
    --controls-dir data/raw/arc2026/final/controls \
    --source-registry configs/source_registry.yaml \
    --output-dir outputs/final/c1 --audit-only
cat outputs/final/c1/panel_audit.md
```

Read the contexts, target count, cells per perturbation, the 0 / 1 / 2 / 3+ counts and
the unsupported targets. **Expect lower coverage than validation.** The H1 source
covers only the 2025 H1 panel; K562 and CD4 are genome-wide. That is information, not a
reason to change the model.

### 4. Prediction and local invariants (~3–4 min for a 300-target panel)

```bash
caffeinate -i uv run python scripts/competition_v2/run_final_panel.py \
    --controls-dir data/raw/arc2026/final/controls \
    --source-registry configs/source_registry.yaml \
    --output-dir outputs/final/c1
```

Step 3 already prepared the panel's sources, and they are reused here. The command
refuses to overwrite an existing `prediction.h5ad`; use a new `--output-dir` to rerun.
All checks must print PASS, and the exit code must be 0.

### 5. Dry-run validation and package (~10 min for a 300-target panel)

```bash
caffeinate -i uv run python scripts/competition_v2/package_final_panel.py \
    --output-dir outputs/final/c1 \
    --controls-dir data/raw/arc2026/final/controls
cat outputs/final/c1/package_manifest.json
```

This runs:

1. compaction (vendored `compact.py`, `--cells-per-target` from the panel);
2. `vcc prep --dry-run` with `--contexts`, `--cells-per-pert` and
   `--expected-gene-dim` taken from the panel; it must exit 0;
3. `vcc.prep` packaging through `vcc_pack_panel.py`.

It stops before packaging if the dry-run fails.

### 6. Record

* Append the `package_manifest.json` SHA-256 values and the `panel_audit.md` summary to
  `reports/research_log.md`.
* **Stop.** Submission is a separate, human decision. After any submission, capture the
  entry id immediately (`vcc submit` prints it) so the official score can be frozen with
  `vcc status <id> --json`.

**Time from controls on disk to validated package:** about 15–20 min for a 300-target
panel.

* **Measured:**
  * A/B/C prediction and invariants: 178 s, with sources reused;
  * mock 111-target prediction: 125 s, including a 43 s K562 re-preparation;
  * mock packaging: 73 s for dry-run plus pack, and about 140 s for compaction.
* **From the C1 logs for 300 targets:** compaction ~5 min, pack ~3.5 min, dry-run
  ~2–3 min.
