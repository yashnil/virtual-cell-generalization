# C5 predeclaration: C1 + X-Atlas (prepared, NOT run)

Written 2026-09-30, **before any C5 result exists and while X-Atlas permission is
PENDING**. Its SHA-256 is stored in
`data/provenance/competition_v2/c5_predeclaration_digest.txt`.
`scripts/competition_v2/run_c5_xatlas.py` checks that digest and refuses to run unless
the permission gate below is open.

## Permission gate (checked in code before anything is read)

C5 may run only if **both** hold:

1. `reports/competition_v2/xatlas_permission_status.md` shows `status` = **APPROVED**,
   and its log records:
   * the sender and organisation (Xaira Therapeutics);
   * the date;
   * the **exact permission text**;
   * whether it explicitly covers the Virtual Cell Challenge 2026;
   * whether it explicitly covers prize eligibility and the ShareAlike question for
     submitted predictions.
2. `licensing.STATUS["HCT116"]` and `licensing.STATUS["HEK293T"]` are `GREEN`. Per the
   register's change control, they change in the same commit as the status file.

The outcomes of the gate are:

* **PENDING:** do not run. Keep C1, and continue paper consolidation.
* **DENIED:** close X-Atlas permanently. Keep C1. The data-expansion line is closed.

## The only change

**C5 = the C1 source universe + X-Atlas HCT116 + X-Atlas HEK293T.**

Everything else is exactly C1:

* **Fusion:** C1a equal fusion, where every usable source (K562, H1, HCT116, HEK293T and
  CD4) has weight 1. This is **not** C0's upstream weights (2/1/1/2/0.5).
* **Unchanged:** the source centering, preprocessing, statistics and minimum-cells rule;
  the promoter correction (GENCODE v47, same pairs); the target panel; the generator,
  amplitude and clip.
* **Nothing new:** no source weights, no shrinkage, no generator or amplitude change.
* **X-Atlas inputs:** the X-Atlas statistics already prepared and frozen for C0
  (`outputs/competition_v2/atlasshift_c0/data/{HCT116,HEK293T}_full_statistics.npz`,
  pinned by `c0_manifest.json`). They are not re-derived.

## Evaluation (the same public folds as C1–C4)

| fold | truth | C1 predictors | C5 predictors | scored with |
|---|---|---|---|---|
| H1 | H1 2025 cells | K562 + CD4 | K562 + CD4 + HCT116 + HEK293T | VCC (unchanged G0 emitter, frozen local ruler) + mean level |
| K562 | K562 GWPS cells | H1 + CD4 | H1 + CD4 + HCT116 + HEK293T | VCC + mean level |
| CD4 | CD4 DE | K562 + H1 | K562 + H1 + HCT116 + HEK293T | effect PDS + mean level |

* **Baseline:** C1 is this run's C1a, which must reproduce C2's `G0_a1.00` raw members
  to ≤ 1e-12.
* **Mean-level metrics:** the C3 definitions (per-target cosine and top-200 sign
  accuracy on the log2fc centred effect, over all predicted cells).

## C5 pass rule (not to be changed after any C5 number exists)

C5 is submission-worthy only if **all** hold, on these folds:

1. The mean Overall (H1, K562; primary ruler, MSE clamped) is greater than C1's.
2. PDS improves or stays ≥ 98 % of C1 in **each** of H1, K562 (scaled PDS) and CD4
   (effect PDS).
3. The mean-response cosine, averaged over the three fold means, does not decrease.
4. At least 2 evaluable folds improve. Improvement means Overall up for H1 and K562,
   and effect PDS up for CD4.
5. No single fold contributes more than 75 % of the summed Overall gain (H1 + K562).

Nothing is tuned after C5 is seen. If C5 passes, the Arc bundle is built with the C1
builder and C5's sources, the C1 checks and `vcc prep --dry-run` are run, and the
package is made. **Submission remains a human decision.**

## Prior expectation (for honesty, not a gate)

C1's X-Atlas ablation used C0's weights, not equal weights. It gave +0.018 mean local
Overall and PDS flat (H1), up (K562) and down (CD4, −0.026 effect PDS). Criterion 2
could therefore fail on CD4.
