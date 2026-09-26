# AtlasShift public data sources — audit (competition track)

Every source listed by AtlasShift's pinned `sources.json` (commit `d24ce4f`), plus the
X-Atlas/Orion scan that `prepare.py` performs outside `sources.json`. Downloads:
`scripts/competition_v2/download_sources.py` → `data/raw/competition_v2/` (git-ignored).
Download log with every verified checksum:
`data/provenance/competition_v2/download_log.jsonl`. Prepared per-source statistics
(upstream `prepare.py`, unmodified): `outputs/competition_v2/atlasshift_c0/data/`.
Coverage matrix: `data/splits/arc_target_support_competition_v2.csv`
(`scripts/competition_v2/build_target_coverage.py`).

## 1. Source table

"Arc ∩" counts official VCC 2026 validation targets with ≥ 20 perturbed cells (AtlasShift's
`minimum_cells`); for CD4 it additionally requires the publisher quality flags. See
§2 for per-source numbers from the prepared statistics.

| source | accession / URL | bytes | SHA-256 (verified) | license | context | technology | perturbations | genes | representation | directly measured? |
|---|---|---|---|---|---|---|---|---|---|---|
| K562 GWPS (Replogle et al. 2022, *Cell*) | figshare+ `10.25452/figshare.plus.20029387.v1`, file 35775507 `K562_gwps_raw_singlecell_01.h5ad` | 65,830,941,948 | `b697ef7f…a9de2c` | **CC BY 4.0** (figshare API) | K562 (CML) | CRISPRi Perturb-seq, 10x 3′, genome-wide | 9,866 targets (1,989,578 cells, 75,328 non-targeting, 267 gem groups) | 8,248 (measured-gene subset) | raw UMI counts, single cell | yes |
| X-Atlas/Orion HCT116 (Xaira Therapeutics 2025, bioRxiv 10.1101/2025.06.11.659105) | HF `slaf-project/X-Atlas-Orion` @ `598aa544…`, `data/HCT116` (SLAF/Lance re-release) | 50,544,250,046 | per-file LFS SHA-256 verified by `huggingface_hub` (no single-file checksum is published) | **CC BY-NC-SA 4.0** | HCT116 (colorectal) | CRISPRi FiCS Perturb-seq (fixed cells), genome-wide | 18,294 | see §2 | raw counts, single cell | yes |
| X-Atlas/Orion HEK293T | same repo, `data/HEK293T` | 82,965,837,799 | as above | **CC BY-NC-SA 4.0** | HEK293T (embryonic kidney) | as above | 18,312 | see §2 | raw counts, single cell | yes |
| VCC 2025 H1 train | `gs://arc-institute-virtual-cell-atlas/virtual-cell-challenge/2025/train/adata_Training.h5ad` | 15,482,497,461 | `a0997710…19c02b` | Arc public release; **no license file in the bucket** (terms unverified) | H1 hESC | CRISPRi Perturb-seq (Arc) | 150 | 18,080 | raw counts | yes |
| VCC 2025 H1 validation | `…/2025/validation/adata_Validation.h5ad` | 6,928,967,541 | `376f0bab…77bfe1` | as above | H1 hESC | as above | 50 | 18,080 | raw counts | yes |
| VCC 2025 H1 test | `…/2025/test/adata_Test.h5ad` | 11,950,739,168 | `ba5ce667…831f03` | as above | H1 hESC | as above | 100 | 18,080 | raw counts | yes |
| CD4 GW DE statistics (Marson lab 2025) | `s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad` | 16,786,240,107 | `c355f535…cfbb62` | **no license file in the bucket** (`data_sharing_readme.md` states none; terms unverified) | primary human CD4⁺ T cells, Rest / Stim8hr / Stim48hr | CRISPRi Perturb-seq, genome-wide | 11,526 contrasts | 10,282 | publisher DE table (log2FC, adj. p, lfcSE) | yes (as DE statistics, not cells) |
| GENCODE v47 annotation | EBI FTP `gencode.v47.annotation.gtf.gz` | 58,953,903 | `df11938c…a3db24` | open | — | genome annotation | — | — | TSS table | n/a |

Full 64-character hashes: `data/provenance/competition_v2/download_log.jsonl`. Every
URL source was downloaded resumably and verified against `sources.json` by our downloader
**and again** by upstream `prepare.py`, which re-hashes before use.

## 2. Arc 300-target coverage per source

From `outputs/competition_v2/target_coverage_summary.json`:

| source | measured (any cells) | usable (≥ 20 cells; CD4 + quality flags) | V1 Tier-0 targets rescued |
|---|---|---|---|
| arch1 | 13 | 13 | 0 |
| kaden25rpe1 | 80 | 80 | 0 |
| nadig25hepg2 | 0 | 0 | 0 |
| nadig25jurkat | 0 | 0 | 0 |
| replogle22k562 | 0 | 0 | 0 |
| replogle22rpe1 | 0 | 0 | 0 |
| wessels23 | 0 | 0 | 0 |
| K562_GWPS | 272 | 269 | 183 |
| XAtlas_HCT116 | 300 | 290 | 204 |
| XAtlas_HEK293T | 300 | 294 | 208 |
| H1_2025_full | 25 | 25 | 6 |
| CD4_DE | 291 | 251 | 179 |

Full VCC 2025 H1 contributes 25 Arc targets (train 13, validation 4, test 8); `arch1` = train only (13), so the complete 2025 release adds **12**.

## 3. Promoter-neighbour pairs (GENCODE v47, upstream rule)

`outputs/competition_v2/atlasshift_c0/data/official_pairs.csv`: **75 pairs, 74 Arc targets,
75 distinct neighbour genes**; TSS distance median 442 bp (5 – 4,821 bp); 40 pairs within
500 bp (capped to 15 % of control). Leakage risk: none (annotation + official lists only).

## 4. Leakage review

| question | answer |
|---|---|
| Does any source contain VCC 2026 perturbed cells? | **No.** 2026 contexts A–F are new Arc screens released only as controls. The 2025 H1 data is Arc's *previous* challenge, fully public since 2025-12. |
| Is any hidden context one of the source cell lines? | **No evidence.** Basal log-CPM correlation of A/B/C controls with each source's control profile is 0.48–0.62, never above the A/B/C cross-correlations (0.48–0.64); HCT116–HEK293T correlate at 0.74 for comparison (`outputs/competition_v2/basal_identity.json`). |
| Were AtlasShift's constants chosen on the A/B/C leaderboard? | **Unknown** — not recorded upstream. This is a *selection* risk (over-fit to A/B/C), not data leakage. |
| Rejected sources | **None** on leakage grounds. |

## 5. License summary for Challenge use

* K562 GWPS: CC BY 4.0 — unrestricted with attribution.
* X-Atlas/Orion: **CC BY-NC-SA 4.0**. Non-commercial use only; ShareAlike attaches to
  adaptations of the data. A Challenge entry by a non-commercial participant is
  consistent with NC; whether a submitted prediction is an "adaptation" under SA is a
  legal question we do not resolve. **Flagged for the human.**
* VCC 2025 H1 and CD4: no license file found. Arc itself states the 2025 data is relevant
  to 2026. **Terms unverified — flagged.**
* The Challenge FAQ allows any data the entrant has rights to, and finalists must disclose
  datasets and public models used.
