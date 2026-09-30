# C4 predeclaration: new direct evidence, source qualification

Written 2026-09-29, **before any statistic, reliability or transfer number was computed
from a new source**. Its SHA-256 is stored in
`data/provenance/competition_v2/c4_predeclaration_digest.txt`, and the C4 scripts refuse
to run if the file has changed.

## Already known when this was written (metadata only)

| source | original record | license | Arc coverage | action |
|---|---|---|---|---|
| **KOLF2.1J** CRISPRi | Figshare+ `10.25452/figshare.plus.27261219.v1` (API JSON saved) | **CC BY 4.0** | 282 / 300 with ≥ 20 cells (read from the original h5ad by range requests) | full-file download (the only file with Arc coverage), MD5-verified |
| **Jurkat** CRISPRi (GSE249595) | GEO; no dataset license. The NCBI disclaimer cannot grant unrestricted permission. The paper is CC BY-NC-ND 4.0 | **UNKNOWN** | 299 / 300 in the guide library (cells not counted) | **not downloaded**: fails criterion 1 |
| **VIPerturb-seq** (K562) | Zenodo `10.5281/zenodo.18460279` | **CC BY 4.0** | 284 / 300 with ≥ 20 cells; median 47 cells vs K562 GWPS 175 | **not downloaded**: gate G, see below |

**VIPerturb-seq gate G.** It is the same context as K562 GWPS, with about 27 % extra cells
per target. By Spearman–Brown (`k = 222 / 175`), K562's median split-half reliability of
0.088 would rise to about 0.11. That is not the "material improvement in K562 response
reliability" the phase requires, so it is not downloaded. It can be revisited as a
K562-replication ablation later.

## KOLF source statistics

* **Same definition as the K562 GWPS source** (`sources._matched_statistics`):
  * the frozen 437-target retained set (300 official targets, then the H1 training
    targets), so centering runs over the same rows as every other source;
  * NTC cells are the controls; KOLF `batch` (ALPHA / BETA / GAMMA) is the matching batch;
  * duplicated symbols are summed;
  * a target is usable with ≥ 20 cells, and the pseudocount shrinkage is the same.
* **Name:** `KOLF2.1J_iPSC`, license GREEN (added to `licensing.STATUS`).
* **Construction:** one streaming pass over the CSC `X` (raw counts) computes per-target
  count sums, mean CPM, per-batch counts and libraries, and the split-half means.

## Reliability audit (§J), same protocol as C3

* **Split-half Spearman–Brown reliability** per target on centred log2fc, genes with
  control mean CPM ≥ 5, targets with ≥ 40 cells. The halves use the exact
  `fusion_c3.split_half` draw: `rng(SEED + 31)`, a permutation per target in
  retained-target order. H1 and K562 are recomputed by the same code.
* **Null:** 100 pseudo-targets of 218 cells (KOLF's median) drawn from NTC cells and
  excluded from the control mean, through the same split-half. **Detectable signal**
  means reliability above the null's 95th percentile.
* **Also reported:** cells per target, response norm, on-target knockdown (the target's
  own-gene log2fc where measured), and main-effect reliability (the SB of the
  half-level Pearson over all targets).
* **Kaden and arch1** are quoted from `reports/kaden_source_reliability_diagnostic_v1.md`
  under **its** protocol and labelled as such. They are not recomputed.

## Held-out transfer value (§K), no reweighting

For each held-out atlas t ∈ {H1, K562, CD4}, the arms are:

* **C1 set** = the other two GREEN atlases with C1's equal-weight fusion (C1a);
* **C1 set + KOLF** = the same rule with KOLF as one more equal-weight source (C1a+KOLF).

Everything else is the C1/C2/C3 harness: the unchanged G0 emitter, promoter cap and
frozen amplitude; the frozen local ruler; and C1a reproducing C2's `G0_a1.00` to
≤ 1e-12. Mean-level metrics use the C3 definitions (cosine, sign accuracy on the top-200
truth genes, norm ratio, energy explained, and CD4 effect PDS). VCC metrics are computed
on H1 and K562.

## Source pass rule (§L; not to be lowered)

KOLF may enter the next competition candidate only if **all** hold:

1. The mean held-out Overall (H1, K562) is greater than C1a's.
2. PDS is preserved: the mean over H1, K562 and CD4-effect is ≥ C1a's, and no single
   fold is below 99 % of C1a.
3. The mean-response cosine (mean of the three fold means) is greater than C1a's.
4. Top-200 sign accuracy (mean of the three fold means) is greater than C1a's.
5. Overall improves in **both** H1 and K562.
6. Not one anomalous fold: neither fold contributes more than 75 % of the summed Overall
   gain.
7. License GREEN.

Plus the §I qualification: direct CRISPRi evidence, harmonisable gene axis, and
measurable signal. The median reliability must exceed the null's 95th percentile.

## §M scientific question (reported, not a gate)

On the H1 fold, KOLF-alone transfer to H1 is compared with K562-alone. On all folds the
KOLF-vs-H1 source agreement is compared with the K562-vs-H1 agreement. This shows
whether KOLF adds redundant pluripotent signal or new transferable signal: does it help
the K562 and CD4 folds, where no pluripotent truth is involved?

No Arc bundle, no new weights, no submission in C4.
