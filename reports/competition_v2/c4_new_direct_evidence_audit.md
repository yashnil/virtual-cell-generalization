# C4: new direct evidence, acquisition and qualification

**Competition track, 2026-09-29/30. Nothing was submitted. No bundle, no model change.**

## Result

**No new source qualifies. KEEP C1 AND WAIT FOR X-ATLAS.**

KOLF2.1J is license-GREEN, covers 282 / 300 Arc targets and carries real signal. As an
extra equal-weight donor, however, it helps **only** the pluripotent held-out context
(H1) and **hurts** K562 and CD4. It fails the predeclared source rule on PDS
preservation, cosine, both-folds and anomalous-fold. Jurkat is license-UNKNOWN.
VIPerturb-seq is GREEN but would not materially improve K562 reliability, so it was not
downloaded.

## Where things are

| artifact | path |
|---|---|
| C3 freeze | `data/provenance/competition_v2/c3_freeze_sha256.txt` (45 files, OK) |
| predeclaration | `reports/competition_v2/c4_predeclaration.md`, SHA-256 `bcedc1d0…` (unchanged) |
| license evidence | `data/provenance/competition_v2/c4/{kolf_figshare_article.json, viperturb_zenodo_record.json, jurkat_guides_reference.fa.gz, viperturb_genome_wide_manifest.txt}`; register updated (`data_license_register.md`, C4 additions); `licensing.STATUS` gains `KOLF2.1J_iPSC` GREEN, `VIPERTURB_K562` GREEN, `JURKAT_GSE249595` UNKNOWN |
| pre-download coverage | `data/provenance/competition_v2/c4/{kolf_pan_genome_obs.json, kolf_strong_obs.json, kolf_arc_coverage.csv}`, read from the **original** Figshare h5ad through HTTP range requests (75 MB fetched) |
| KOLF download | `data/raw/competition_v2/kolf/KOLF_Pan_Genome_QC_Filtered.h5ad`, 189,393,177,972 bytes, MD5 `afd30fde1e6ad32969c29868394385d1` = Figshare `supplied_md5`; logged in `download_log.jsonl` |
| code | `src/virtual_cell/competition_v2/sources_c4.py` (streaming CSC → C1-format statistics; tests `tests/test_competition_v2_c4.py`, including exact equality with the C3 split-half); `scripts/competition_v2/{c4_remote_h5ad, c4_download, run_c4_kolf}.py` |
| outputs | `outputs/competition_v2/c4_kolf/{KOLF_full_statistics.npz, KOLF_split_halves.npz, reliability_*.csv/json, transfer_mean_level.csv, kolf_single_source_and_agreement.csv, folds/*, vcc_scaled.csv, c4_decision.json}` |
| figures | `reports/competition_v2/figures/c4_P{1..4}_*.{png,svg}` + `figures/sources/c4_P*.csv` |

## Answers

### 1. KOLF license status

**GREEN: CC BY 4.0.** Evidence:

* **Record:** the original Figshare+ API for article 27261219
  (`10.25452/figshare.plus.27261219.v1`, published 2026-05-18) returns
  `license: CC BY 4.0`. Attribution is required.
* **Assay:** single-gene CRISPRi. The obs field `perturbation_type` is `CRISPRi` for all
  2,659,209 cells. Controls are the 146,747 `NTC` cells, in 3 batches (ALPHA, BETA,
  GAMMA).
* **Data:** raw integer counts in `X` (CSC, float32-typed integers), with the same counts
  duplicated in `layers/counts`. The gene axis has 37,567 features. Target-level
  pseudobulk was reconstructed; obs `total_counts` equals the X row sums.

The six files (sizes and MD5) are in the saved API response. The file sizes differ from
the brief: the pan-genome file is 189.4 GB, not 176 GB. The 46.7 GB
"Strong Perturbations" file holds only 1,655 strong perturbations and covers just 55
Arc targets, so it is insufficient. The pan-genome file was the smallest sufficient
artifact.

### 2. KOLF Arc-target overlap (decision gate D, from the original file's obs)

**282 / 300 Arc targets, all with ≥ 20 cells** (median 236.5 cells; 279 with ≥ 50).
17,971 of the 18,533 official genes are on its axis.

Not measured: CHKB, DTNBP1, EPHB2, IFNGR2, IL10RB, LENG1, MTM1, PBLD, PSMB9, RSRC1,
SNRK, STK3, TBC1D19, TBCK, TESK2, TRIP4, TTBK2, ZC2HC1A.

### 3. KOLF unique target contribution

| usable sources per Arc target | C1 | C1 + KOLF |
|---|---|---|
| 0 | 13 | **2** |
| 1 | 53 | 18 |
| 2 | 210 | 55 |
| 3+ | 24 | **225** |

* KOLF uniquely adds **11 of the 13** targets with no license-clean evidence: ABCD1,
  ANKRD52, CAPRIN2, KIF21B, NICN1, PARP3, PHF19, SEMA4F, SLC44A1, SNN, TAF4. EPHB2 and
  PBLD stay uncovered.
* 46 targets gain a second source, and 225 a third or more.

### 4. KOLF response reliability (same C3 protocol for every row)

| source | protocol | targets | **median reliability** | IQR | null p95 | above null | cells per target | on-target log2FC (median) |
|---|---|---|---|---|---|---|---|---|
| H1 2025 | C3 split-half SB | 149 | **0.787** | 0.71–0.89 | 0.061 | 100 % | 1,071 | −3.90 |
| **KOLF2.1J** | C3 split-half SB | 411 | **0.160** | 0.13–0.24 | 0.035 | **100 %** | 234 | −1.27 |
| K562 GWPS | C3 split-half SB | 384 | **0.088** | 0.04–0.15 | 0.040 | 76 % | 177 | −2.30 |
| CD4 DE | C3 SE-based | 346 | 0.000 | 0–0 | — | — | — | — |
| Kaden RPE1 | *quoted, research-track protocol* | 1,836 | 0.167 | | | | 400 | |
| arch1 (H1) | *quoted, research-track protocol* | 150 | 0.906 | | | | 1,045 | |

* **Signal:** every KOLF target clears the NTC pseudo-target null. The main-effect
  (panel-mean) reliability is **0.94**.
* **Per-target reliability** is about twice K562's and similar to Kaden's.
* **Knockdown** is real but milder than in H1 or K562.

### 5. KOLF held-out transfer value (§K: same C1 equal-fusion rule, no reweighting)

The C1 set plus KOLF, compared with the C1 set alone. C1a reproduces C2's `G0_a1.00`
exactly in both folds.

| fold | arm | PDS | MSE | NMAE | FID | REACH | JAC | **Overall** |
|---|---|---|---|---|---|---|---|---|
| H1 | C1 set | 0.626 | −0.036 | −0.120 | −0.122 | −0.012 | 0.048 | 0.070 |
| H1 | + KOLF | **0.778** | **+0.012** | **−0.077** | **−0.107** | −0.020 | 0.050 | **0.106** |
| K562 | C1 set | 0.437 | −0.037 | 0.117 | −0.158 | 0.136 | 0.001 | 0.089 |
| K562 | + KOLF | **0.393** | −0.067 | 0.096 | −0.225 | 0.087 | −0.001 | **0.059** |

Mean-level metrics, on all predicted cells:

| held out | cosine | top-200 sign acc. | effect PDS |
|---|---|---|---|
| H1 | 0.065 → 0.064 (0.065 → **0.088** on the C1-covered cells) | 0.564 → **0.588** | — |
| K562 | 0.045 → **0.026** | 0.530 → 0.517 | — |
| CD4 | 0.031 → **0.016** | 0.520 → 0.513 | 0.699 → **0.646** |

**Rule L** (predeclared):

| criterion | result |
|---|---|
| 1. mean Overall up | ✓ 0.0793 → 0.0823 |
| 2. PDS preserved | ✗ K562 90 %, CD4 92 % |
| 3. cosine up | ✗ 0.047 → 0.035 |
| 4. sign accuracy up | ✓ +0.0013 (H1 only) |
| 5. both folds up | ✗ H1 +0.036, K562 −0.030 |
| 6. no anomalous fold | ✗ 100 % of the gain is from H1 |
| 7. GREEN | ✓ |

**KOLF FAILS.**

**§M, redundant or new signal? Neither, cleanly: it is pluripotency-specific.**

* **It adds real signal beyond H1-less C1 for a pluripotent truth.** H1 PDS rises from
  0.626 to 0.778, MSE turns positive, and sign accuracy rises by 3.4 points on the same
  cells.
* **The signal does not transfer.** KOLF alone reaches K562 or CD4 at cosine ≈ **0.007**,
  against 0.029–0.082 for the existing donors.
* **It agrees with H1 and nothing else.** KOLF~H1 source agreement is 0.066–0.078,
  against 0.009–0.014 with K562 and CD4.
* **Fusion then dilutes the transferable donors** in non-pluripotent contexts.

The Arc contexts are six non-pluripotent cell lines (FAQ) and are basally unlike every
source (Fig. P4). So on the Arc panel KOLF should behave like the K562 and CD4 folds,
not like H1.

### 6. Jurkat license status

**UNKNOWN. Not downloaded; fails criterion 1.**

* **Repository:** GEO GSE249595 (a SubSeries of GSE247601; BioProject PRJNA1049794) has
  no dataset license.
* **NCBI's disclaimer:** it "cannot provide comment or unrestricted permission
  concerning the use, copying, or distribution" of GEO data whose submitters may claim
  IP.
* **The paper** (Nat Cell Biol 2025, PMC11906366) is CC BY-NC-ND 4.0. That covers the
  article, not the data, and the authors include Myllia Biotechnology employees.

### 7. Jurkat Arc-target overlap

The library (`GSE249595_guides_reference.fa.gz`, 83,401 guides, 18,601 targets plus
non-targeting) contains **299 / 300** Arc targets, including all 13 targets with no C1
evidence. These are library guides, not measured cells: cells per perturbation and the
unstimulated/stimulated split require the 3.5 GB RAW.tar, which was deliberately not
downloaded.

### 8. VIPerturb-seq license and overlap

**GREEN: CC BY 4.0.**

* **Record:** Zenodo `10.5281/zenodo.18460279`, published 2026-02-02.
* **Manifest (MD5-verified):** 18,888 perturbations. It covers **all 300** Arc targets,
  284 with ≥ 20 cells, at a median of **47 cells per perturbation**. The filtered object
  has 6,727 perturbations and includes 124 Arc targets.
* **Format:** the data are Seurat `.rds` files, 3 × ~3.5 GB plus a 3.6 GB filtered
  object. R is not installed here.

### 9. Does VIPerturb add independent K562 value?

**Not materially.** It is the same context as K562 GWPS (median 177 cells, reliability
0.088). Adding about 47 cells per target (×1.27) lifts the Spearman–Brown reliability to
about 0.11, and that is before any cross-protocol batch variance. It fails gate G, so it
was **not downloaded**.

It would, however, give K562-context evidence for all 13 C1-unsupported targets. If
coverage rather than reliability becomes the goal, that is a separate, predeclared
ablation.

### 10. X-Atlas permission status

**PENDING.** No written reply is recorded (`xatlas_permission_status.md`, log entry
2026-09-29). X-Atlas was not used.

### 11. Which new sources pass qualification?

**None.**

* **KOLF:** GREEN, measurable signal, fails the transfer rule.
* **Jurkat:** UNKNOWN license.
* **VIPerturb-seq:** GREEN, fails the value gate, not downloaded.

### 12. Exact source universe for next

**Unchanged: K562 GWPS + VCC 2025 H1 + CD4 DE (+ GENCODE for the promoter cap)**, the C1
universe.

### 13. Is a new competition model justified?

**No.** None of C2, C3 or C4 found a change that passes its predeclared rule. The only
candidate with evidence of broad transfer value is X-Atlas (+0.018 public Overall in C1,
two non-pluripotent contexts), and it is BLOCKED pending permission.

## Caveats

* **Different reliability protocols.** Kaden and arch1 are quoted under the research
  track's protocol. All other rows share the C3 protocol, including H1 and K562,
  recomputed.
* **Basal similarity (P4) is platform-dominated.** The Arc contexts are 10x Flex probe
  data, and KOLF sits nearer K562 than H1 on centred basal profiles. It is shown for
  completeness and was not used for any decision.
* **Response norms depend on the gene axis** (KOLF's is 37k features), so they are not
  compared across sources.
* **A context-gated use of KOLF was not tested**, e.g. as a donor only for pluripotent
  targets. That would be a new model, which this phase excludes.

## Validation

* `uv run pytest -q`: **604 passed, 1 skipped** (2 new C4 tests).
* `ruff check` and `ruff format --check`: clean.
* `uv build`: OK.
* `verify_c1_state.py --full` (`outputs/competition_v2/c4_kolf/state_full_end.json`):
  **0 failures**.
* C2 and C3 freeze manifests: 0 failures.
* C1 package, statistics and folds untouched.
