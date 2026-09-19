# scPertEval four-context dataset — source audit and pipeline specification

Status: **audit complete; nothing downloaded; awaiting download approval.**
Date: 2026-09-18

**This is not a Molina & Zhang reproduction.** It is an *independent
four-context re-derivation* on standardized public data. See §6.

---

## 1. Source audit (Task 1)

### 1.1 Provenance

| item | value |
|---|---|
| Repository | https://github.com/Virtual-Cell-Research-Community/scPertEval |
| **Commit SHA audited** | **`4685f11927e887745737600170da7a655b727553`** (`main`, release v0.2.0) |
| Commit date | 2026-09-09T16:39:25Z |
| **Retrieval date** | **2026-09-18** |
| License | MIT |
| Docs read | `docs/user-guide/datasets.md`, `docs/references.bib`, `src/scperteval/dataset.py` (all pinned at the SHA above) |
| Data bucket | `gs://scperteval/processed/` — public, read-only |
| HTTPS mirror | `https://storage.googleapis.com/scperteval/processed/<dataset>_processed_complete.h5ad` |

All four files return HTTP 200 over plain HTTPS; `gsutil` is not required.

### 1.2 Files, sizes, checksums

Sizes and MD5s from the GCS JSON API on 2026-09-18.

| dataset | cell line | exact filename | bytes | GB | GiB | MD5 (base64) |
|---|---|---|---:|---:|---:|---|
| `replogle22k562` | K562 | `replogle22k562_processed_complete.h5ad` | 2,430,512,332 | 2.431 | 2.264 | `zzGVsPPynk6Inru+fcMVZQ==` |
| `replogle22rpe1` | RPE1 | `replogle22rpe1_processed_complete.h5ad` | 1,877,364,555 | 1.877 | 1.748 | `S2VTPbbJ8p75C8i3aD1l3g==` |
| `nadig25hepg2` | HepG2 | `nadig25hepg2_processed_complete.h5ad` | 1,236,448,196 | 1.236 | 1.152 | `AmgRZbNFVwWbOYlZII8yLw==` |
| `nadig25jurkat` | Jurkat | `nadig25jurkat_processed_complete.h5ad` | 2,004,474,709 | 2.004 | 1.867 | `2lpruoqyqbt+9hARAQcAbg==` |
| **total (4 files)** | | | **7,548,799,792** | **7.549** | **7.030** | |

The bucket holds 7 datasets / 17.763 GB in total; we need only these four.

### 1.3 Contents — verified directly, not taken from the docs

Read from the remote files themselves via **HTTP range requests, 24.2 MB fetched
of 7,549 MB** (HDF5 metadata only; no expression data transferred).

| property | `replogle22k562` | `replogle22rpe1` | `nadig25hepg2` | `nadig25jurkat` |
|---|---|---|---|---|
| cells (`n_obs`) | 308,646 | 240,774 | 133,757 | 258,202 |
| genes (`n_vars`) | **8,563** | **8,749** | **9,623** | **8,881** |
| control cells | 10,691 | 11,485 | 4,976 | 12,013 |
| perturbed cells | 297,955 | 229,289 | 128,781 | 246,189 |
| unique perturbations (excl. control) | 1,971 | 2,016 | 1,818 | 2,137 |
| min cells per perturbation | **30** | **30** | **30** | **30** |
| median cells per perturbation | 125 | 82 | 55 | 89 |
| max cells per perturbation | 1,996 | 3,580 | 1,213 | 2,555 |

Cell, control and perturbation counts agree **exactly** with the published
summary table in `docs/user-guide/datasets.md`.

| structural property | value (identical in all four) |
|---|---|
| `X` representation | `csr_matrix`, **`float32`**, log-normalised |
| `obs` columns | **`perturbation` only** |
| `var` columns | none — gene names are the index |
| `layers` / `obsm` / `varm` / `obsp` / `varp` / `uns` | **all empty** |
| **raw counts retained?** | **NO — nowhere in the file** |
| perturbation label format | bare HGNC gene symbol, e.g. `AAAS`, `AARS2`, `ABCB10` |
| combination perturbations | none (`+` never appears in any of the four) |
| control label | the single literal string `control` |
| gene identifiers | **HGNC symbols** (e.g. `LINC01409`, `NOC2L`, `ISG15`) — not Ensembl IDs |

**Raw counts are not recoverable from these files.** `X` is log-normalised and
the raw layer was deliberately dropped during trimming. Any analysis needing
counts (e.g. a negative-binomial noise model) would have to go back to the
original deposits.

### 1.4 Preprocessing — A (inherited) vs B (performed by scPertEval)

**A. Inherited from the original experiments/deposits**

* Assay, chemistry and screen design: CRISPRi Perturb-seq, 10x 3'.
  Replogle 2022 contributes its **essential-gene panel** (~2,057 targets) for
  K562 and RPE1 — explicitly *not* the genome-wide K562 screen (`K562_gwps`,
  ~9,866 targets), which scPertEval does not host. Nadig 2025 contributes the
  companion HepG2 and Jurkat screens on the same essential-gene protocol.
* Guide-to-cell assignment, doublet/quality calls, and the upstream cell
  filtering that produced the deposited matrices.
* **The ≥30 cells-per-perturbation floor is inherited** — the observed minimum
  is exactly 30 in all four datasets, so it was applied upstream, not by
  scPertEval (whose own default is `min_cells`, applied at runtime).
* Gene detection per dataset: each dataset carries a *different* gene set
  (8,563–9,623), reflecting per-dataset expression filtering upstream.

**B. Performed by scPertEval when building the hosted files**

1. **Label cleaning** — one perturbation label per cell; all non-targeting
   controls collapsed to the single label `control`; guide/plasmid suffixes
   stripped; combinations joined with `+`; cells with a missing label dropped.
   Per the docs, `replogle22k562` needed *no* cleanup — its labels were already
   clean gene symbols.
2. **Log-normalisation** — `sc.pp.normalize_total(target_sum=1e4)` then
   `sc.pp.log1p`. So `X` is `log1p(CP10K)`.
3. **Light QC** — `sc.pp.filter_cells(min_genes=200)`,
   `sc.pp.filter_genes(min_cells=3)`.
4. **Trimming** — everything except `X`, `obs["perturbation"]` and the gene
   names removed: raw-count layer, `obsm` embeddings, `uns`, surplus columns.
   Written sparse `float32`, gzip-compressed.

**No HVG selection, no scaling, no batch correction, no PCA** is applied by
scPertEval. This is important: unlike Molina & Zhang's `obsm["X_hvg"]` — whose
construction is undocumented and unavailable — every preprocessing step here is
specified and reproducible.

### 1.5 Source publications and accessions

| cell line | publication | DOI | primary deposit |
|---|---|---|---|
| K562, RPE1 | Replogle et al., *Mapping information-rich genotype–phenotype landscapes with genome-scale Perturb-seq*, **Cell** 185(14):2559–2575.e28 (2022) | [10.1016/j.cell.2022.05.013](https://doi.org/10.1016/j.cell.2022.05.013) | figshare+ [10.25452/figshare.plus.20029387](https://doi.org/10.25452/figshare.plus.20029387) |
| HepG2, Jurkat | Nadig et al., *Transcriptome-wide analysis of differential expression in perturbation atlases*, **Nat. Genet.** 57(5):1228–1237 (2025) | [10.1038/s41588-025-02169-3](https://doi.org/10.1038/s41588-025-02169-3) | GEO `GSE264667` |

Note the year correction: Molina & Zhang's config calls this "Nadig et al. 2024"
(the 2024 preprint/deposit); the published article is **2025**, which is what
scPertEval cites.

### 1.6 Relationship to Molina & Zhang — corroboration, not equivalence

The scPertEval control-cell counts are **identical** to the control counts
reported in Molina & Zhang's paper:

| | K562 | RPE1 | HepG2 | Jurkat |
|---|---|---|---|---|
| Molina & Zhang, reported controls | 10,691 | 11,485 | 4,976 | 12,013 |
| scPertEval, measured controls | **10,691** | **11,485** | **4,976** | **12,013** |

So both derive from the same upstream deposits with the same control definition.
They diverge on everything downstream: M&Z report far fewer total cells and
perturbations (e.g. K562 188,590 cells / 1,383 perturbations vs 308,646 /
1,971 here), via a filtering step their repository does not document, and their
response space is an undocumented 2,000-gene `X_hvg` rather than the full
log-normalised matrix.

**Therefore results computed here are not expected to equal theirs, and a
mismatch is not evidence of a bug in either.**

---

## 2. Storage preflight (Task 2)

```
$ df -h ~
Filesystem      Size    Used   Avail Capacity  Mounted on
/dev/disk3s5   1.8Ti   620Gi   1.2Ti    35%    /System/Volumes/Data
```

| quantity | value |
|---|---|
| Download size, four files | **7.549 GB (7.030 GiB)** |
| Free space | 1.2 TiB |
| Download as a share of free space | **0.6 %** |
| Current `data/raw/` usage | 641 MB (Arc controls) |
| Peak requirement (download only, no decompression needed) | ~7.6 GB |

Files are gzip-compressed HDF5 and are read in compressed form; no
decompressed copy is required on disk. Working memory, not disk, is the real
constraint: see §4.9.

**Verdict: storage is safe with very large margin.**

---

## 3. Gate status amendment (Task 3)

| track | status |
|---|---|
| **Exact Molina & Zhang reproduction** | **BLOCKED** — their claimed processed artifacts are not publicly present (see `reports/molina_zhang_reproduction_spec.md` §1.2). Not resolvable from our side. |
| **Independent four-context re-derivation** | **READY** — standardized public K562/RPE1/HepG2/Jurkat data located, audited and costed. |

Naming rule, to be applied in every report, log entry, figure caption and commit
message from here on: this work is the **"independent four-context
decomposition"**. It must never be described as a Molina & Zhang reproduction,
and its numbers must never be presented as reproducing theirs.

---

## 4. Proposed independent decomposition pipeline (Task 4)

### 4.1 Response space — shared gene set

The four datasets do **not** share a gene space (8,563 / 8,749 / 9,623 / 8,881).
Measured intersections:

| quantity | value |
|---|---|
| **genes shared by all four** | **6,640** |
| union across the four | 11,909 |
| pairwise shared | 7,226–7,733 |

**Decision: use the 6,640-gene intersection as the primary response space, with
all 6,640 genes retained initially.** Rationale: the intersection is a function
of gene *identity* only — it never touches expression values — so it is
deterministic, leakage-free, and independent of which cell line is later held
out. Gene order frozen as the sorted symbol list.

Cost of the intersection: each dataset loses 22–31 % of its genes. That is the
price of a balanced design and must be reported.

### 4.2 Shared perturbation set

| quantity | value |
|---|---|
| **perturbations shared by all four** | **1,264** |
| union | 2,374 |
| pairwise shared | 1,511–1,872 |

Because the upstream floor is already ≥30 cells everywhere, **every min-cells
threshold from 1 to 30 yields the same 1,264 perturbations** — the threshold is
not a free parameter in the range that matters (M&Z used 5 and 10; both are
no-ops here). Raising it further costs a lot: 50 → 577; 100 → 96.

Of the 1,264 shared perturbations, **1,072 target a gene that is itself in the
6,640-gene response space** and 192 do not. Those 192 can still be decomposed;
only their on-target knockdown check is unavailable.

**Decision: freeze the 1,264 × 6,640 design and commit both lists as files
under `data/splits/` before any analysis is run.**

### 4.3 Whether to use all genes initially — yes

Use all 6,640. HVG selection is exactly the undocumented step that makes M&Z
irreproducible; re-introducing it as a primary choice would repeat their
mistake. HVGs belong in the sensitivity analysis (§4.8), not the main result.

### 4.4 Control pseudobulk

Per cell line `c`:

```
ctrl_mean[c] = mean over all cells with perturbation == "control", restricted to the 6,640 shared genes
```

One vector per cell line (10,691 / 11,485 / 4,976 / 12,013 cells respectively).
Controls are pooled, never split by guide — the files carry no guide identity.

### 4.5 Perturbation pseudobulk and the response

```
pert_mean[c, p] = mean over all cells with perturbation == p in cell line c   (shared genes)
delta[c, p]     = pert_mean[c, p] - ctrl_mean[c]
```

### 4.6 Response transformation — none

`X` is already `log1p(CP10K)`, so `delta` is a difference of mean log-expression,
i.e. a log-fold-change-like quantity on a per-cell-line basis. **No further
centring, scaling, z-scoring or re-normalisation.** Adding one would be an
untracked degree of freedom; if any is wanted it goes in §4.8.

One caveat to state in the write-up: the mean of `log1p` is not the `log1p` of
the mean, so this is a "mean of log" pseudobulk, not a library-size-weighted
one. This matches how scPertEval and M&Z both work, and it is applied
identically to every cell line, so it cannot bias the cross-context comparison.

### 4.7 Decomposition and denominator

Use `virtual_cell.decomposition.anova`, already implemented and verified by 41
tests:

* `delta = mu + alpha_c + beta_p + gamma_{c,p}` by balanced two-way cell means.
* Denominator: the **uncentred** `mean_{c,p} ||delta[c,p]||²`, with all four
  shares (`mu`, `alpha`, `beta`, `gamma`) reported and summing to 1.
* Report `template = mu + alpha` as one share **and** `mu` and `alpha`
  separately — M&Z merge them, and keeping both lets us compare either way.
* Cell-line axis order frozen as `(k562, rpe1, hepg2, jurkat)`.

### 4.8 Split-half reliability

Per resample, for each (cell line, perturbation): shuffle that perturbation's
cells, take two **disjoint equal-size** halves (drop the odd cell), average each
half over the shared genes, and subtract the **full-data** control mean — the
controls are *not* split, so the two halves share a reference and the estimate
reflects perturbation-side noise only. Decompose each half independently and
combine matching components by **cross-half dot product**; average over
resamples; noise is the residual share.

Feasible everywhere: the smallest perturbation has 30 cells → 15 per half.
HepG2 is the thinnest (median 57 cells among shared perturbations) and should be
reported separately, since its noise share will be the largest.

Settings: 50 resamples, fixed seed, both recorded. 200 resamples as a stability
check.

### 4.9 Leakage avoidance

* The perturbation and gene intersections use **identity only, never expression
  values** — so freezing them cannot leak a held-out cell line's responses.
* The decomposition itself is descriptive, not predictive: no train/test split
  exists to leak across, and nothing is fitted or tuned.
* **For the later transfer experiments** the rule is stricter and must be stated
  now: any step that *looks at expression* — HVG selection, PCA, scaling,
  feature standardisation — must be fitted on **source cell lines only** and
  applied to the held-out line. Selecting HVGs on all four pooled would leak.
* A held-out cell line contributes only its control cells; its perturbation
  responses never enter training or model selection.
* Freeze `data/splits/shared_perturbations.txt` and
  `data/splits/shared_genes.txt` before analysis; never regenerate them
  conditional on a result.
* **No tuning toward 27.8 / 29.4 / 23.5 / 19.3.** Those numbers come from a
  different cell set, gene space and filtering, and are not a target.

### 4.10 Memory note

The largest file is 308,646 × 8,563 sparse float32. Pseudobulk means are
computed by streaming row chunks and accumulating per-perturbation sums, reusing
the pattern already proven on the Arc controls in
`virtual_cell.data.arc2026.stream_row_chunks`. Split-half needs per-cell access
within a perturbation; those subsets are small (≤3,580 cells) and can be
materialised one perturbation at a time. **The full matrix is never densified.**

---

## 5. Gate criteria (Task 5)

The independent four-context re-derivation **passes** when all five hold.
Matching the paper's percentages is explicitly **not** required.

| # | criterion | how it is judged |
|---|---|---|
| 1 | **Decomposition mathematics passes all invariants** | The 41 existing tests plus the same invariants asserted on the real tensor at runtime: exact reconstruction to float tolerance, all zero-sum side conditions, the four sums of squares partitioning the total, pairwise component orthogonality. |
| 2 | **Results are stable across reasonable preprocessing choices** | Component shares move by less than a pre-declared tolerance across: min-cells 30 vs 50 vs 100; gene space 6,640-intersection vs source-selected 2,000 HVGs vs union-with-missing-dropped; response as raw log-difference vs per-gene standardised; cell subsets (all four vs each 3-subset); 50 vs 200 split-half resamples; and ≥3 seeds. Declare the tolerance **before** running. |
| 3 | **Conserved beta and interaction gamma can be quantified** | Both have a well-defined, non-degenerate share with a reported uncertainty across seeds/resamples, and per-perturbation `beta_fraction` has a sensible spread rather than collapsing to 0 or 1. |
| 4 | **Split-half reliability distinguishes reproducible signal from noise** | The noise share is strictly between 0 and 1, decreases monotonically as cells per perturbation increase (already verified as a property of the estimator on synthetic data), and the cross-half estimate is stable across seeds. A negative or >1 share fails the gate. |
| 5 | **Results do not depend on a pathological preprocessing choice** | No single filtering/normalisation decision flips a qualitative conclusion. In particular the ordering of the component shares must be stable; if it is not, report that as the finding rather than picking a favourable setting. |

If criteria 1–5 pass, the gate passes **as an independent re-derivation** and
modelling may proceed. It never becomes a Molina & Zhang reproduction.

**No novel prediction architecture is to be implemented until the gate passes.**

---

## 6. Arc separation

Unchanged and honoured. Arc A/B/C were not touched by this task. The paper
response space has **not** been intersected with the Arc 18,533-gene panel; that
remains a separate, later stage to be run only after this gate passes.

## 7. Download decision

**Nothing has been downloaded.** Only 24.2 MB of HDF5 metadata was read over
range requests to produce §1.3 and §4.

Awaiting approval to fetch the four files (7.549 GB) to `data/raw/scperteval/`,
with MD5 verification against §1.2 and a provenance record written to
`data/provenance/` in the same form as the Arc bundle.
