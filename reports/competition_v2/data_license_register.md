# Competition-v2 data license register

**Date accessed: 2026-09-26.** Competition track. Permission is **never** inferred from
public accessibility. A source is GREEN only when a named license, or an explicit
rights-holder statement, permits use in a prize-bearing Challenge entry. The
machine-readable mirror is `src/virtual_cell/competition_v2/licensing.py`. C1 code refuses
any non-GREEN source (`assert_sources_allowed`).

Method. A research agent fetched every page listed below. The key quotes were then
re-checked in this session:

* CC0 on the Arc Virtual Cell Atlas page;
* "MIT License" on the CZI Virtual Cells Platform page;
* the 2026 datasets-page permission for H1;
* the EMBL-EBI and Ensembl statements;
* the figshare API `license` field;
* the Hugging Face license tag of X-Atlas/Orion.

The local X-Atlas copy also ships `LICENSE.md` (CC BY-NC-SA 4.0 full text).

This register is an engineering record, not legal advice.

## Summary

| source | licensor | license | Challenge-use status | used in C1 |
|---|---|---|---|---|
| VCC 2025 H1 hESC (train / validation / test) | Arc Research Institute | CC0 1.0 | **GREEN** | yes |
| K562 genome-wide Perturb-seq (Replogle 2022) | Replogle & Weissman (figshare+) | CC BY 4.0 | **GREEN** (attribution) | yes |
| CD4⁺ T-cell GWCD4i DE statistics (Marson lab 2025) | Zhu, Dann et al.; listed by CZ Biohub Virtual Cells Platform | MIT License (platform listing) | **GREEN** (caveat below) | yes |
| GENCODE v47 annotation (promoter-neighbour TSS table) | GENCODE / EMBL-EBI / Ensembl | no named license; EMBL-EBI Terms of Use and Ensembl "no restrictions" | **GREEN** (caveat below) | yes (promoter cap) |
| X-Atlas/Orion HCT116 | Xaira Therapeutics | CC BY-NC-SA 4.0 | **BLOCKED_PENDING_PERMISSION** | **no** |
| X-Atlas/Orion HEK293T | Xaira Therapeutics | CC BY-NC-SA 4.0 | **BLOCKED_PENDING_PERMISSION** | **no** |
| Kaden 2025 RPE1 (Southard et al., Nat Genet 2025; scPertEval `kaden25rpe1`) | Southard, Norman lab (MSKCC) | CC BY 4.0 (Zenodo) | **GREEN** (attribution) | no (scientific exclusion) |
| VCC 2026 A/B/C control bundle | Arc Research Institute | Challenge Terms of Use (no open license) | Challenge use only; do not redistribute | yes (as the Challenge input) |

**UNKNOWN: none.** Every source AtlasShift uses was resolved. The only other direct source
in the repository, the v1 scPertEval datasets, contributes no Arc coverage beyond `arch1`
(= H1 2025 train, CC0) and Kaden.

## Governing Challenge terms (VCC 2026)

The site is client-rendered. Its text was read from the page bundles. The rules are "as
revised on September 16, 2026".

> "You may not use any property, know-how, code, information, or data in the Challenge
> that you do not have rights to use for that purpose." — https://virtualcellchallenge.org/rules

> "As long as you have the appropriate permissions to use the data, you may use any data to
> improve the performance of your predictions." — https://virtualcellchallenge.org/faq

Finalists must publish "the datasets you used in training" (rules). The Challenge carries
a cash prize. We therefore treat a NonCommercial license as incompatible with an entry
until the licensor says otherwise.

## Per-source records

### VCC 2025 H1 hESC — GREEN

* URLs: `https://storage.googleapis.com/arc-institute-virtual-cell-atlas/virtual-cell-challenge/2025/{train,validation,test}/…`;
  evidence at https://arcinstitute.org/tools/virtualcellatlas and https://virtualcellchallenge.org/datasets.
* License evidence (Virtual Cell Challenge 2025 block, arcinstitute.org/tools/virtualcellatlas):
  > "© 2025. This work is openly licensed via CC0 1.0." (links to creativecommons.org/publicdomain/zero/1.0)
* Challenge organisers, 2026 datasets page:
  > "You are free to train on the H1 hESC data released for the 2025 Challenge, and on any
  > other public or proprietary data you have the right to use"
* Attribution: not legally required (CC0). We cite the Arc Cell commentary.
* Redistribution: unrestricted.
* Gaps: the bucket and the `ArcInstitute/arc-virtual-cell-atlas` GitHub repo carry no
  LICENSE file. These are gaps, not contradictions.
* Status: **GREEN_FOR_VCC_USE** (the entering policy, confirmed; no contradictory evidence).

### K562 GWPS (Replogle et al. 2022, *Cell*) — GREEN

* URL: figshare+ `10.25452/figshare.plus.20029387.v1`, file 35775507.
* License evidence (https://api.figshare.com/v2/articles/20029387, `license` field,
  re-fetched 2026-09-26):
  > `{'value': 1, 'name': 'CC BY 4.0', 'url': 'https://creativecommons.org/licenses/by/4.0/'}`
* Attribution: **required**. Cite Replogle J., Weissman J. (2022), figshare+ DOI above, and
  Replogle et al. 2022 *Cell*.
* Redistribution: allowed with attribution and an indication of changes. No ShareAlike.
* Supporting: the VCC 2026 datasets page lists this record as a recommended public dataset.

### CD4⁺ T-cell GWCD4i DE statistics (Marson lab 2025) — GREEN with caveat

* URL: `https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/GWCD4i.DE_stats.h5ad`.
* License evidence (https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq,
  dataset metadata panel):
  > "Version v1.0 , processed released 22 Dec 2025 License MIT License Repository
  > https://github.com/emdann/GWT_perturbseq_analysis_2025/…"

  and the same page:
  > "Do not use the dataset for the following purposes: Discriminatory or biased analyses
  > Any use that is not in accordance with the Acceptable Use Policy . Any use prohibited by
  > the MIT License ."
* The associated GitHub repository (`emdann/GWT_perturbseq_analysis_2025`) has an MIT
  LICENSE ("Copyright (c) 2025 Emma Dann"). The README asks users of the data or code to
  cite Zhu R., Dann E. et al. (2026) *Cell*.
* Caveat:
  * The MIT text is written for "Software".
  * The authors' own `data_sharing_readme.md` in the bucket names no license.
  * The MIT label *for the dataset* comes from the official hosting platform the preprint
    points to. The bioRxiv data-availability statement links to that CZI page.
  * The entering policy was "GREEN if the official source continues to show MIT licensing".
    It does, so the status is GREEN.
  * Recommended: record the provenance in the finalist disclosure, and optionally ask the
    authors to confirm.
* Attribution: keep the MIT notice and cite Zhu, Dann et al.
* Redistribution: allowed. No ShareAlike. The CZI Acceptable Use Policy applies.

### GENCODE v47 (promoter-neighbour correction) — GREEN with caveat

* URL: `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz`.
* GENCODE names no license of its own. Its "Terms of use" link resolves to the EMBL-EBI
  terms:
  > "EMBL-EBI itself places no additional restrictions on the use or redistribution of the
  > data available via its Data Resources and Tools other than those provided by the
  > original data owners" — https://www.ebi.ac.uk/about/terms-of-use/
* The GENCODE human annotation is the Ensembl/HAVANA annotation:
  > "Ensembl imposes no restrictions on access to, or use of, the data provided and the
  > software used to analyse and present it." — Ensembl legal disclaimer (read from
  > grch37.ensembl.org; www.ensembl.org returned HTTP 403)
* Only derived facts are used: TSS coordinates and distances between genes. Attribution is
  expected by good practice.
* Status: **GREEN**. It is the explicit "no restrictions" statement of the producing
  project and the hosting institute, not an inference from accessibility. The promoter
  correction is therefore retained exactly as frozen (§14 of the C1 report).
* AtlasShift's MIT code license was **not** relied on for this data.

### X-Atlas/Orion HCT116 and HEK293T — BLOCKED_PENDING_PERMISSION

* URLs: https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion (original);
  `slaf-project/X-Atlas-Orion` @ `598aa544…` (the SLAF re-release we downloaded).
* License evidence:
  * Hugging Face API tag `license:cc-by-nc-sa-4.0` (original repository, re-fetched
    2026-09-26, lastModified 2025-09-19).
  * The re-release README:
    > "- **License**: CC-BY-NC-SA-4.0 (Creative Commons Attribution-NonCommercial-ShareAlike 4.0)"
  * Local `data/raw/competition_v2/xatlas_orion/LICENSE.md`: the full CC BY-NC-SA 4.0 text.
* Restrictions:
  * NonCommercial. A prize-bearing entry is at least arguably "directed towards …
    monetary compensation".
  * ShareAlike applies to adapted material.
  * Attribution is required.
* Status: **BLOCKED_PENDING_PERMISSION** until Xaira Therapeutics grants written permission.
  See `xatlas_permission_status.md`.
* Use in this phase:
  * No X-Atlas response or statistic enters C1.
  * X-Atlas statistics are read only by the predeclared ablation arms (`C0_withX`), which
    measure what permission would be worth, and by the implementation-equality tests.
  * Both are non-commercial research uses, and neither produces a submission.

### Kaden 2025 RPE1 — GREEN license, excluded from C1

* Source: Southard K.M. … Norman T.M., "Comprehensive transcription factor perturbations
  recapitulate fibroblast transcriptional states", *Nat Genet* 57:2323 (2025).
* Deposits: Zenodo records 15213597, 15200179, 15211972, 15213619, 15215389, 15215414,
  15215154 and 15215216. The Zenodo API `metadata.license` is `cc-by-4.0` for all eight.
* Not verified: which record scPertEval built `kaden25rpe1` from. All candidates are
  CC BY 4.0.
* Scientific exclusion (not a licensing one):
  * It is not part of the C0 backbone being reproduced.
  * Our frozen results found its responses barely reproduce (median split-half
    reliability 0.17).
  * The deposit names an RPE1 **CRISPRa** TF screen, and its modality relative to the
    CRISPRi Arc screens was not verified.

### VCC 2026 control bundle

Governed by the Challenge Terms of Use (https://virtualcellchallenge.org/terms):
> "Your use of our Content is subject to these Terms"

It is used only as the Challenge input and never redistributed.

## C4 additions (2026-09-29)

### KOLF2.1J iPSC genome-scale CRISPRi atlas: GREEN

* **Status:** GREEN.
* **Source:** "A genome-scale CRISPRi perturbation atlas of human induced pluripotent
  stem cells", *Nat Biotechnol* (2026), `10.1038/s41587-026-03199-w`.
* **Data record:** Figshare+ `10.25452/figshare.plus.27261219.v1`, article 27261219,
  published 2026-05-18.
* **License evidence:** the original Figshare API returns
  `license = {"name": "CC BY 4.0", "url": "https://creativecommons.org/licenses/by/4.0/"}`.
  The saved response is `data/provenance/competition_v2/c4/kolf_figshare_article.json`.
* **Files, sizes and MD5s** are recorded from the same response.
* **Used:** `KOLF_Pan_Genome_QC_Filtered.h5ad`, 189,393,177,972 bytes,
  MD5 `afd30fde1e6ad32969c29868394385d1`.
* **Attribution required.** `licensing.STATUS["KOLF2.1J_iPSC"] = GREEN`.

### Jurkat genome-scale CRISPRi (GSE249595 / PRJNA1049794): UNKNOWN

* **Status:** UNKNOWN.
* **Repository:** GEO series GSE249595, a SubSeries of GSE247601. It carries **no
  dataset license**.
* **NCBI GEO disclaimer:** "NCBI places no restrictions on the use or distribution of the
  GEO data. However, some submitters may claim patent, copyright, or other intellectual
  property rights in all or a portion of the data they have submitted. NCBI … cannot
  provide comment or unrestricted permission concerning the use, copying, or
  distribution of the information contained in GEO."
* **Publication:** Nat Cell Biol 2025 (PMC11906366) is licensed CC BY-NC-ND 4.0. That
  covers the article, not the dataset. Its data statement only points to GEO.
* **Commercial authors:** several are Myllia Biotechnology employees.
* **Consequence:** without an explicit permissive dataset license it is UNKNOWN, and
  `licensing.STATUS["JURKAT_GSE249595"] = UNKNOWN`. Not downloaded.

### VIPerturb-seq (K562 genome-wide): GREEN

* **Status:** GREEN. **Not downloaded** (see the C4 report, gate G).
* **Record:** Zenodo `10.5281/zenodo.18460279` (concept `…18460278`), published
  2026-02-02.
* **License evidence:** Zenodo API `metadata.license.id = "cc-by-4.0"`. The saved
  response is `data/provenance/competition_v2/c4/viperturb_zenodo_record.json`.
* `licensing.STATUS["VIPERTURB_K562"] = GREEN`.

## Change control

A status may change only with new written evidence, and only by a dated addition to this
file. `licensing.STATUS` must change in the same commit. X-Atlas may move to APPROVED only
through `xatlas_permission_status.md`.
