# Molina & Zhang reproduction — frozen specification

Status: **Task 1 complete. Tasks 2–3 BLOCKED on unavailable reference data.**
Date: 2026-09-18 (amended 2026-09-18, see §7)

> **Amendment.** An exact reproduction of this paper remains blocked and is not
> being pursued further. The project has moved to an *independent four-context
> re-derivation* on standardized public scPertEval data — a separate track with
> its own gate. See `reports/scperteval_four_context_data_spec.md` and §7 below.
> Nothing in that track may be described as reproducing this paper.

| item | value |
|---|---|
| Paper | Alexis Molina & Xinyi Zhang, *Perturbation response decomposition enables biologically aligned generalization to unseen perturbations and cellular contexts* |
| DOI | [10.64898/2026.07.24.740459](https://doi.org/10.64898/2026.07.24.740459) |
| Preprint server | bioRxiv, posted 2026-07-27, CC-BY 4.0 |
| Official code | https://github.com/xinyizhanglab/perturbation-decomposition |
| **Commit SHA audited** | **`a15214780619736d393f40240e56ba992fd416a3`** (branch `main`, the only branch) |
| Repo last pushed | 2026-07-24T06:53:20Z (4 commits total: `Init repo`, 2x `Update README`, `Clean up`) |
| **Retrieval date** | **2026-09-18** |
| Repo size | 189 KB; 75 tree entries; no releases, no tags, no forks, empty wiki |
| Org | `xinyizhanglab` — one public repo only |
| License | MIT (stated in README; no `LICENSE` file in tree) |

A read-only copy was extracted to the session scratchpad and made immutable
(`chmod -R a-w`). **The authors' repository was not modified.**

---

## 1. Repository audit (Task 1)

### 1.1 What exists

| Requested item | Present? | Location at `a152147` |
|---|---|---|
| ANOVA/decomposition implementation | **Yes (two, and they disagree)** | `decomposition/anova.py`; `figures/fig1_f.py` |
| Split-half / reliability implementation | **Yes, but only inside figure scripts** | `figures/fig1_f.py` L115–190; `figures/fig1_suppl_decomp.py`; `figures/fig1_supplementary_full.py` (S1e) |
| Zero-shot cross-cell-line benchmark code | Yes | `benchmarks/cross_cl.py` |
| Ridge and MLP benchmark code | Yes | `models/ridge.py` (both); driven by `benchmarks/{within,cross}_cl.py` |
| Evaluation metrics | **Partial — module is broken** | `evaluation/metrics.py` |
| Raw-data preparation scripts | **Partial — upstream step missing** | `data/process.py`, `data/prepare_release.py` |
| **Processed genetic-response files** | **NO — absent** | claimed `data/processed/genetic/{k562,rpe1,hepg2,jurkat}.pkl` |
| **DepMap embeddings** | **NO — absent** | claimed `data/processed/genetic/depmap_embeddings.pkl` |

### 1.2 Blocker: the released artifacts do not exist

The README states:

> "Data is included (processed pseudobulk + DepMap embeddings, ~90 MB)"
> "All predictions are included in `results/`. No retraining needed."

Both claims are false for the published repository. `.gitignore` excludes
exactly those paths:

```
data/raw/
data/processed/
*.h5ad
*.pkl
results/
```

Verified absent: no `data/processed/`, no `results/` in the git tree; 0 releases;
0 tags; 0 forks; the organisation has no second repository; no Zenodo or
figshare deposit was found. The repository is 189 KB of source only.

**Consequence:** there is no reference fixture. Tasks 2 and 3, as specified
("use the authors' released processed artifacts"), cannot be executed.

### 1.3 Blocker: 34 of ~45 Python files hardcode the authors' cluster

Every figure script — including the ones that compute the paper's headline
decomposition — reads from an absolute path on the first author's cluster:

```python
BENCH = Path("/mnt/storage01/home/amolina/benchmark_2026")
...
a = ad.read_h5ad(BENCH / "data" / cl / "replogle.h5ad")
X = np.array(a.obsm["X_hvg"])
```

`data/prepare_release.py` is explicitly documented "Run ON THE CLUSTER", and also
reads a sibling private tree `agop-perturbation/morph_pipeline/` that is not
public. The DepMap embeddings ship from that private pipeline, not from the
public `build_depmap_embeddings` path.

### 1.4 Blocker: the upstream preprocessing step is absent

Both `data/process.py` and `data/prepare_release.py` **consume** a pre-existing
`obsm["X_hvg"]`; neither creates it. Nothing in the repository performs
normalisation, log-transformation, or HVG selection. `configs/datasets.yaml`
says `n_hvgs: 2000` and nothing more. So even with the raw data in hand, the
authors' processed artifacts cannot be regenerated exactly: the normalisation
scheme, HVG flavour, batch handling, and whether `X_hvg` is log-scale or
z-scored are all unspecified.

### 1.5 Defects found in the released code

| # | File | Defect |
|---|---|---|
| D1 | `evaluation/metrics.py` | `evaluate_predictions` calls `centroid_accuracy(...)` and `median_rank(...)`, **neither of which is defined anywhere in the repository**. The function raises `NameError`. It is the exact snippet in the README's Quick Start step 4. |
| D2 | `configs/datasets.yaml` | The `nadig_hepg2:` block is mis-indented, so it parses to `None` and its 8 keys leak as siblings of the dataset entries. Verified: `datasets['nadig_hepg2'] is None`. HepG2 is unconfigurable as shipped. |
| D3 | `configs/datasets.yaml` | `raw_url` for **both** Replogle K562 and RPE1 is `plus.figshare.com/articles/dataset/.../21999546`. Figshare article 21999546 is *"Characteristics of the Greifer and the Axon-Hook"*, an unrelated PLoS ONE figure. The correct record is **20029387** (see §4). |
| D4 | repo-wide | `decomposition/anova.py` reads `{cl}.pkl`; `data/process.py` writes `{cl}_pseudobulk.pkl`. Filenames do not match. |
| D5 | repo-wide | Minimum cells per perturbation is **5** in `process.py`/`prepare_release.py`/`configs`, but **10** in `figures/fig1_f.py`, which is the script that produces the paper's numbers. |
| D6 | `decomposition/anova.py` | The documented decomposition entry point does **not** implement the paper's decomposition (see §2.4). |
| D7 | `configs/datasets.yaml` | `models.state.repo: https://https://github.com/...` (doubled scheme). |

---

## 2. The decomposition, as actually implemented (Task 2 specification)

**Authority: `figures/fig1_f.py`, not `decomposition/anova.py`.** The module
advertised in the README produces neither the paper's components nor its
denominator, and has no noise correction. The figure script is what generates
the reported percentages.

### 2.1 Response space and preprocessing

| property | value | source |
|---|---|---|
| Input | `data/{cl}/replogle.h5ad`, one per cell line | `fig1_f.py` L49 |
| Response matrix | `adata.obsm["X_hvg"]` | `fig1_f.py` L50 |
| Response dimensionality | 2,000 HVGs (`n_hvgs: 2000`) — **not verifiable from code; the HVG step is absent** | `configs/datasets.yaml` |
| Perturbation label column | `obs["gene"]` | `fig1_f.py` L51 |
| Control label | `"non-targeting"` (exact string) | `fig1_f.py` L52 |
| Control reference | `ctrl_mean = X[genes == "non-targeting"].mean(0)` — a single per-cell-line vector | `fig1_f.py` L52 |
| Minimum cells per perturbation | **≥ 10** | `fig1_f.py` L57 |
| Response definition | `delta[c,p] = mean(cells of p in c) - ctrl_mean[c]` | `fig1_f.py` L75 |
| Scale | difference of means in `X_hvg` space. **Not** a log-fold-change, **not** re-centred, **not** re-scaled | `fig1_f.py` L75 |
| Cell-line order | `["k562", "rpe1", "hepg2", "jurkat"]` (declaration order, *not* sorted) | `fig1_f.py` L26 |
| Shared perturbation set | `sorted(set.intersection(...))` over the four per-cell-line perturbation sets | `fig1_f.py` L64 |
| Shared perturbation count | **not stated in the repo; printed at runtime only** | — |

### 2.2 Decomposition

```python
mu    = D.mean(axis=(0, 1))                                  # (G,)
alpha = D.mean(axis=1) - mu                                  # (C, G)
beta  = D.mean(axis=0) - mu                                  # (P, G)
gamma = D - mu - alpha[:, None, :] - beta[None, :, :]        # (C, P, G)
```

Plain balanced two-way ANOVA by cell means. Consequences, all of which we test:

* **Exact reconstruction:** `D == mu + alpha + beta + gamma` identically.
* **Zero-sum side conditions:** `sum_c alpha = 0`, `sum_p beta = 0`,
  `sum_c gamma = sum_p gamma = 0`.
* **Orthogonality:** the four components are mutually orthogonal under §2.3.
* **No constraint is imposed or optimised** — there is no fitting, no
  regularisation, no iterative solve.

### 2.3 Sums-of-squares denominator — the convention that matters

```python
ss_total = np.mean([np.sum(D[c, p] ** 2) for c in range(n_cl) for p in range(n_p)])
ss_mu    = np.sum(mu ** 2)
ss_alpha = np.mean([np.sum(alpha[c] ** 2) for c in range(n_cl)])
ss_beta  = np.mean([np.sum(beta[p] ** 2)  for p in range(n_p)])
ss_gamma = np.mean([np.sum(gamma[c, p] ** 2) for c, p in ...])
```

* The total is an **uncentred mean squared norm per observation**, *not*
  `np.var`. The `mu` component therefore carries a real share.
* Each component is averaged over the axes it does not span.
* In a balanced design these four sum **exactly** to `ss_total`.

### 2.4 Why `decomposition/anova.py` is not the paper's decomposition

The README's documented entry point differs materially:

| | `figures/fig1_f.py` (paper) | `decomposition/anova.py` (README) |
|---|---|---|
| Denominator | `mean_{c,p} \|\|D[c,p]\|\|²` (uncentred) | `np.var(D)` (centred on a **scalar** grand mean) |
| `mu` share | reported | **never reported** |
| Noise correction | split-half | **none** |
| Reported components | template(`mu`+`alpha`), `beta`, `gamma`, noise | `alpha`, `beta`, `gamma` only |
| Cell-line order | declaration order | `sorted()` → `hepg2, jurkat, k562, rpe1` |
| Fractions sum to 1 | yes | **no** (the `mu` share is silently dropped) |

Running `python -m decomposition.anova` **cannot** produce 27.8 / 29.4 / 23.5 /
19.3 even with the data. Anyone reproducing from the README alone would fail.

### 2.5 Split-half noise correction

```python
N_SPLITS = 50
rng = np.random.RandomState(42)   # one RNG threaded through all 50 resamples
```

Per resample, for every (cell line, perturbation):

1. `perm = rng.permutation(len(cells))`, `half = len(cells) // 2`
2. half 1 = `cells[perm[:half]]`, half 2 = `cells[perm[half:2*half]]` — an odd
   cell is **dropped** so the halves are exactly equal in size
3. each half's delta is referenced against the **full-data** `ctrl_mean`;
   **the controls are not split**
4. each half is decomposed independently (§2.2)
5. matching components are combined by **cross-half dot product**:
   `dot(mu1, mu2)`, `mean_c dot(a1[c], a2[c])`, `mean_p dot(b1[p], b2[p])`,
   `mean_{c,p} dot(g1[c,p], g2[c,p])`

Independent measurement noise has zero expected cross-product, so these estimate
the **reproducible** SS of each component on the same scale as §2.3. Averaged
over the 50 resamples, then:

```python
template_pct = (sig["mu"] + sig["alpha"]) / ss_total * 100
beta_pct     =  sig["beta"]               / ss_total * 100
gamma_pct    =  sig["gamma"]              / ss_total * 100
noise_pct    = (ss_total - sum(sig.values())) / ss_total * 100
```

**Noise is a residual, not a direct estimate.** "Template" is `mu + alpha`
combined — this is why the paper's "cell-line-specific" share is a single number.

### 2.6 Template-removed residual (second bar)

Projective, not subtractive:

```
eps[c,p] = D[c,p] - (D[c,p] · That_c) That_c,   That_c = mean_p D[c,p] / ||mean_p D[c,p]||
```

Applied **independently to each half** before the split-half decomposition, then
the same four shares are recomputed against `ss_resid_total`.

### 2.7 Reported values to reproduce (targets — never to be tuned toward)

| component | full response | template-removed residual |
|---|---|---|
| template (`mu`+`alpha`) | **27.8 %** | ~0 by construction |
| conserved `beta` | **29.4 %** | — |
| interaction `gamma` | **23.5 %** | — |
| noise | **19.3 %** | **~37.4 %** |

Source: paper Figure 1c/1f. Sum = 100.0 %.

---

## 3. Zero-shot transfer benchmark (Task 3 specification)

`benchmarks/cross_cl.py` with `--n_sources 3` is exactly the four
leave-one-cell-line-out settings requested:

| setting | sources | target |
|---|---|---|
| 1 | rpe1 + hepg2 + jurkat | k562 |
| 2 | k562 + hepg2 + jurkat | rpe1 |
| 3 | k562 + rpe1 + jurkat | hepg2 |
| 4 | k562 + rpe1 + hepg2 | jurkat |

* **Features** `V`: per-perturbation DepMap dependency embedding, looked up as
  `emb[p]` or `emb[p + "+ctrl"]`. Perturbations without an embedding are dropped
  from **both** train and test.
* **Targets** `D`: the pseudobulk delta vectors of §2.1.
* **Training set**: all (pert, delta) pairs from the three source cell lines,
  **pooled**, with each perturbation appearing once per source cell line.
* **Ridge**: `PCA(n_components=min(50, n-1, d))` → `Ridge(alpha=10)`.
* **MLP**: same PCA → `Linear→ReLU` ×2 with hidden `(256, 256)`, output `n_genes`;
  Adam `lr=5e-4`, `weight_decay=1e-5`, 400 epochs, batch 64, cosine schedule,
  grad-norm clip 1.0, MSE loss; **averaged over 3 seeds** (`n_mlp_seeds=3`;
  note `configs` says `n_seeds: 3`, README table says 5).
* **Metrics**: `pref` (primary) and `std_r`.
  * `std_r` = mean over test perturbations of `pearson(delta, delta_hat)`.
  * `pref` = mean over test perturbations of
    `pearson(delta - T, delta_hat - T)`.
  * **`T = tD.mean(0)`, the *target* cell line's own mean delta** — *not* the
    training centroid, despite the docstring in `evaluation/metrics.py` saying
    "mean delta of TRAINING perturbations". Record this deviation: it means the
    cross-CL metric removes a target-derived template.

### 3.1 Not present in the repository

* **beta-only transfer** — no implementation. Not in `benchmarks/`, not in
  `models/`. (Trivial to define from §2.2, but it is *our* construction, not
  theirs, and must be labelled as such.)
* **beta recovery / gamma recovery metrics** — only a passing reference in
  `figures/fig4_c.py` ("interaction recovery corr(delta_hat, gamma)"); no
  reusable implementation.
* `centroid_accuracy`, `median_rank` — referenced, never defined (D1).

---

## 4. Raw data acquisition plan (Task 4)

**No downloads have been started. Approval is required before any of these.**

### 4.1 Verified sources

| dataset | cell line | source | accession | exact file | size | local destination | preprocessing role |
|---|---|---|---|---|---|---|---|
| Replogle 2022 | K562 | figshare+ | [10.25452/figshare.plus.20029387.v1](https://doi.org/10.25452/figshare.plus.20029387) | `K562_essential_raw_singlecell_01.h5ad` | **10.66 GB** | `data/raw/replogle2022/` | pseudobulk delta source |
| Replogle 2022 | RPE1 | figshare+ | same record | `rpe1_raw_singlecell_01.h5ad` | **8.70 GB** | `data/raw/replogle2022/` | pseudobulk delta source |
| Nadig 2024 | HepG2 | GEO | GSE264667 | `GSE264667_hepg2_raw_singlecell_01.h5ad` | **5.2 GB** | `data/raw/nadig2024/` | pseudobulk delta source |
| Nadig 2024 | Jurkat | GEO | GSE264667 | `GSE264667_jurkat_raw_singlecell_01.h5ad` | **8.7 GB** | `data/raw/nadig2024/` | pseudobulk delta source |
| DepMap 2024Q2 | 1,208 lines | DepMap portal | 2024Q2 | `CRISPRGeneEffect.csv` | ~0.1 GB | `data/raw/depmap2024q2/` | perturbation features `V` |

**Total ≈ 33.4 GB.** All four h5ad files are the `raw_singlecell` variants;
the `normalized_*` and `*_bulk_*` variants are **not** what the pipeline reads.

Direct figshare file endpoints (stable, verified 2026-09-18):
`K562_essential_raw_singlecell_01.h5ad` → `https://ndownloader.figshare.com/files/35773219`;
`rpe1_raw_singlecell_01.h5ad` → `https://ndownloader.figshare.com/files/35775606`.

Corrections to the authors' config: the figshare URL in `configs/datasets.yaml`
is wrong for both Replogle lines (D3); the correct record is 20029387. The GEO
filenames in the config **are** correct and were confirmed against the GEO FTP
listing.

### 4.2 Verified preprocessing parameters

| parameter | value | confidence |
|---|---|---|
| assay | CRISPRi Perturb-seq, 10x 3' | high (paper) |
| perturbation type | single-gene CRISPRi knockdown | high |
| `obs` perturbation column | `gene` | high (code) |
| control definition | `gene == "non-targeting"`, pooled into one mean vector per cell line | high (code) |
| min cells per perturbation | 10 for the decomposition; 5 for the released pseudobulk | high (code, but inconsistent — D5) |
| pseudobulk construction | arithmetic mean over cells of `X_hvg`, minus control mean | high (code) |
| HVG count | 2,000 | medium (config only) |
| **normalization** | **UNKNOWN** | **none — step absent from repo** |
| **HVG selection method** | **UNKNOWN** (flavour, batch handling, per-line vs shared) | **none — step absent from repo** |
| **DepMap embedding construction** | **UNKNOWN** — released embeddings come from the private `agop-perturbation/morph_pipeline`; the public `build_depmap_embeddings` yields raw 1,208-dim profiles instead | **none** |

### 4.3 Unresolved count discrepancy

The paper reports per-line cell and perturbation counts:

| line | cells | controls | perturbation genes |
|---|---|---|---|
| K562 | 188,590 | 10,691 | 1,383 |
| RPE1 | 173,737 | 11,485 | 1,499 |
| HepG2 | 96,616 | 4,976 | 1,340 |
| Jurkat | 184,470 | 12,013 | 1,537 |

These are substantially **smaller** than the raw deposits (the Replogle essential
screens alone carry ~2,000 targeted genes and far more cells). An unspecified
filtering/subsampling step sits between the raw files and `replogle.h5ad`.
Note also that `prepare_release.py` names the per-cell-line input `replogle.h5ad`
for **all four** lines including HepG2 and Jurkat, which suggests a single
harmonised re-processing whose provenance is not documented.

**We cannot claim a faithful raw-data reproduction until this is resolved.**

---

## 5. What we implemented (Task 2, partial)

`src/virtual_cell/decomposition/anova.py` — our own implementation, written to
§2.2/§2.3/§2.5/§2.6 above. It is verified by 41 tests in
`tests/test_decomposition.py`, including:

* exact reconstruction (general and degenerate shapes)
* all zero-sum side conditions
* the balanced-design SS partition and pairwise component orthogonality
* **planted-component recovery** — components and their SS budget are planted in
  a synthetic tensor and recovered exactly, so correctness is established
  without the authors' data
* permutation invariance in both the perturbation and cell-line axes, and
  covariance of `beta` with a perturbation permutation
* split-half behaviour: unbiasedness on identical halves, monotone increase of
  the noise share as cells per perturbation fall, seed reproducibility
* frozen perturbation set, frozen cell-line order, frozen response-gene set
* rejection of unbalanced designs, missing cell lines, empty intersections,
  malformed shapes, and non-finite values
* **verbatim equivalence with the reference implementation's algebra**, by
  transcribing `fig1_f.py` L78–90, L94–104 and L131–144 into the test file and
  asserting agreement to 1e-12

No number here is tuned toward 27.8 / 29.4 / 23.5 / 19.3. Those remain
unverified targets.

---

## 6. Gate status

| gate | status |
|---|---|
| Task 1 — repository audit | **PASS (complete)** |
| Task 2 — our decomposition implemented and mathematically verified | **PASS** |
| Task 2 — numeric reproduction against the authors' processed data | **BLOCKED** — no such data exists |
| Task 3 — zero-shot baselines reproduced | **BLOCKED** — same |
| Task 4 — raw acquisition plan | **PASS (complete, awaiting download approval)** |

**The Molina & Zhang reproduction gate is NOT passed, and is now closed as
unachievable from public materials.** The hard rule stands in amended form: no
interaction model or other novel architecture until the *independent
four-context gate* passes (§7).

### Arc separation

Honoured throughout. Arc A/B/C were not touched by this task: not trained on,
not used to choose any preprocessing, not intersected with the paper response
space, and no identity inference was attempted. The Arc-alignment stage remains
a separate, later step.


---

## 7. Amendment (2026-09-18): the reproduction track is closed, a re-derivation track replaces it

| track | status | rationale |
|---|---|---|
| **Exact Molina & Zhang reproduction** | **BLOCKED — closed** | The processed artifacts their README calls "included" are excluded by their own `.gitignore`, with no release, tag, fork or external deposit, and the upstream `X_hvg` preprocessing is absent from the repository. Not resolvable from our side. |
| **Independent four-context re-derivation** | **READY** | Standardized public K562/RPE1/HepG2/Jurkat data audited at `reports/scperteval_four_context_data_spec.md`. |

Per the user's direction on 2026-09-18: **do not contact the authors and do not
download the ~33.4 GB of raw deposits** listed in §4. That acquisition plan is
retained below for the record but is **not** the active path.

### What survives from this document

Still authoritative and reused by the new track:

* §2 — the decomposition, denominator, split-half and projective-template
  conventions, which our `virtual_cell.decomposition.anova` implements and which
  are verified by 41 tests. These remain the *methodological* reference even
  though the numeric comparison is off the table.
* §1.5 — the released-code defects, which stand as findings.
* §4.1 — the corrected figshare record (20029387) and confirmed GEO filenames,
  now superseded in practice by the much smaller scPertEval processed files
  (7.549 GB vs 33.4 GB) but still the correct raw provenance.

### What must never be claimed

The reported values in §2.7 (27.8 / 29.4 / 23.5 / 19.3) are **not a target** for
the new track. They were computed on a different cell set, a different and
undocumented 2,000-gene response space, and a different filtering regime. Our
independent numbers are expected to differ, and a difference is **not** evidence
of a bug in either analysis.
