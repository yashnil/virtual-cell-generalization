# External benchmark candidates — audit record, 2026-09-21

Audited from primary sources (GEO, Figshare API, ENA, publisher pages) using
metadata, supplementary statistics tables and HTTP range requests. **Neither
dataset itself was downloaded.** Incidental audit files only, listed below.

## Candidate 1 — TeloHAEC endothelial CRISPRi Perturb-seq

**Accession correction.** The phase plan named GSE212396. That series is
*"[Pilot scRNA-seq]"* — 50-gene and 200-gene libraries used to optimise
conditions for, and cross-validate, the larger study. The ~2,285-gene screen is
**GSE210681**.

| | |
|---|---|
| accession | GSE210681 (pilot: GSE212396) |
| context | CRISPRi TeloHAEC, telomerase-immortalised human aortic endothelial |
| design | 36,880 promoter-targeting guides, CROP-seq, 20 10x lanes |
| targets | 2,345 (released log2FC matrix); CAD GWAS loci ±500 kb |
| genes | 25,054 (p-value matrix) / 17,472 (NMF spectra) |
| cells | 32,991 in one audited lane × 20 ≈ 660,000 (~280/target) |
| control | non-targeting guides; log2FC already referenced to control |
| representation | log2 fold change; raw counts only inside `RAW.tar` |
| publication | PMID 38326615 |
| assets | log2fcs 475 MB · Pvalues 415 MB · RAW.tar 6.4 GB · RDS 1.4 GB |

Arc/overlap: 66 Arc targets · 15,321 Arc panel genes (82.7%) · 6,613 of our
6,640-gene axis · 250 measured-transfer controls · 1,925 unseen.

**Rejected on signal.** Two independent released statistics agree that
per-perturbation transcriptional effects are weak:

| statistic | result |
|---|---|
| program-level MAST, experiment-wide FDR < 0.05 | 304/2,357 (12.9%) with ≥1 significant program |
| gene-level, BH per perturbation, FDR < 0.05 | 93/2,345 (4.0%) with ≥10 significant genes; median **0** |

Usable subset: 48 unseen, 20 measured-transfer controls, 3 Arc targets. This is
consistent with the authors' own design — they aggregate into 60 NMF programs
because individual CAD-locus perturbations are weak in endothelium.

## Candidate 2 — Feng et al., iPSC CRISPRi across cell lines

| | |
|---|---|
| publication | Cell Genomics 2026; PMC12903452; bioRxiv 2024.11.28.625833 |
| context | human iPSC, CRISPRi |
| targeted screen | 444 genes, 1,355 guides (3/gene), 20 non-targeting |
| cell lines | 20 lines / 10 donors (genome-scale: 34 / 26) |
| cells per target per line | median 74 (pooled over lines ≈ 1,480) |
| genes | 6,520 (targeted screen) |
| raw data | ENA ERP165335 |
| processed | Figshare 10.6084/m9.figshare.26819743 (9.12 GB, 11 files) |
| counts | Figshare 10.6084/m9.figshare.27989294 (7.93 GB, 9 files) |
| licence | MIT |

Arc/overlap: 5 Arc targets · 6,194 Arc panel genes (33.4%) · 5,213 of our
6,640-gene axis · 170 measured-transfer controls · 141 unseen.

**Signal: 444/444 targets (100%) have ≥10 significant genes at adjusted
p < 0.05; median 189 of 6,520 (2.9%).**

## Selection

**Feng targeted screen**, on the criterion Stage A added: a benchmark whose own
responses do not reproduce cannot distinguish a failed method from a failed
measurement. `kaden25rpe1` passed every structural criterion and produced an
uninterpretable result at reliability 0.115; TeloHAEC would repeat that.

Requested download (pending approval): `TargetedScreen_LFC_byGene-perLine.tsv.gz`,
**1.92 GB** — per-line deltas for 444 targets × 20 lines, giving a cross-line
reliability estimate, the pooled estimate, and the degradation-across-contexts
test in one file. Preferred over the 4.48 GB count matrix.

## Incidental audit files fetched

Held in the session scratchpad, not the repository:

| file | size | purpose |
|---|---|---|
| Feng `TargetedScreen_LFC_byGene.tsv.gz` | 77.6 MB | 444-target list, gene axis, signal statistics |
| GSE210681 `ALL_Pvalues_dup4_s4n3.99x.txt.gz` | 415 MB | gene-level signal, BH per perturbation |
| GSE210681 `Table.S2.X.AllMASTStatisticalTestResults.txt.gz` | 5.4 MB | program-level signal |
| GSE210681 `2kG.library.gene_spectra_score...txt.gz` | 10 MB | measured-gene axis |
| GSE210681 log2FC header (ranged GET) | 4 MB | 2,345-target list |
| GSM6435403 barcodes + guide calls | 0.6 MB | cells per lane |
