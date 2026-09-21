# Feng multi-context external validation — v1

**Feng confirms the arch1 failure, and this time the positive control is alive.**

Across 19 iPSC lines screened under one protocol, the frozen unseen-perturbation
predictor scores `r = +0.007` (median over lines) and `+0.018` pooled —
indistinguishable from zero in every line. On the *same* lines, the *same*
targets' measured counterparts transfer at `r = +0.128` per line and **`+0.320`
pooled**, positive in **19 of 19 lines**.

That contrast is the result. Kaden could not deliver it because its positive
control was dead; here both arms are measured on the same material, so the gap
between them cannot be blamed on the target measurement.

Nor is it a property of the target sets. Matched on the dataset's own
per-target signal strength, **measured transfer improves steadily with signal
(+0.084 → +0.452 across quartiles) while prior-based prediction stays flat at
+0.01 to +0.03.** Better measurement helps direct transfer and does nothing for
prior interpolation. They are qualitatively different strategies.

**This is predeclared Interpretation 2.** One honest qualification: Feng cannot
test Interpretation 1, because its 19 lines are all iPSC and all sit within
0.019 of each other in basal similarity to the training contexts. **Context
distance barely varies, so context-distance dependence remains untested** — and
what *does* explain the between-line variation is measurement reliability
(`r = +0.956` against the measured-transfer control).

**No model was fitted and no constant changed.** The scoring path reproduces the
frozen arch1 result to **9.0e-17** before any Feng line is scored. Tier-0 policy
is unchanged.

Date: 2026-09-21
Reproduce: `uv run python scripts/run_feng_multicontext.py`
Outputs: `outputs/feng_multicontext_v1/`
Provenance: `data/provenance/feng/feng.md`

---

## 0. The protocol is provably the frozen one

Feng ships log fold changes, not cells, so the frozen cell-level entry point
does not apply. `virtual_cell.modelling.multicontext` scores a delta matrix
instead — **importing** `PRIOR_FAMILY`, `KNN_K`, `RIDGE_ALPHA`,
`PCA_COMPONENTS` and `prepare_matrix` from the frozen modules rather than
redefining them (a test asserts no local copy exists).

The equivalence is checked, not claimed: arch1's pseudobulk is pushed through
the new path first.

```
max |difference| against frozen arch1 beta-level result: 9.021e-17
REPRODUCED. The scoring path is the frozen one.
```

## 1. How many cell lines are usable?

**19 lines, all scored.** 53,490,080 rows; 42,767,452 on our gene axis; **zero
duplicated** `(line, target, gene)` rows; 8,204 of 8,436 `(line, target)` cells
present (97.2%); 2.75% missing values overall.

| | value |
|---|---|
| cell lines | **19** (the paper describes 20 for the targeted screen) |
| targets | **444** |
| output genes on our 6,640-gene axis | **5,213** |
| globally unseen targets | **141** |
| directly measured in the sources | **170** |
| STRING coverage of the unseen set | 0.766 |

All counts were recomputed locally and match the pre-download audit.

Two lines are usable only nominally: `tolg_4` has **0** targets with ≥10
significantly changed genes and `fiaj_3` has 6. They are reported rather than
dropped, and flagged.

## 2. How reliable are they?

**There are no split halves and no replicate guides in this file, so true
split-half reliability cannot be computed.** Two things can, and they are kept
strictly distinct from it.

**(a) A lower bound from cross-line agreement.** For two lines measuring the
same perturbation, `r_ij ≤ sqrt(rho_i · rho_j)`, so `rho_i·rho_j ≥ r_ij²`. The
bound is loose by exactly the amount context genuinely changes the response.
Median cross-line `r` runs **0.009–0.046**, giving reliability floors of
0.0001–0.0022 — a floor so low it constrains nothing. It is reported for
completeness, not as an estimate.

**(b) The dataset's own signal strength.** Per line, the number of targets with
≥10 significantly changed genes (`pval_adj < 0.05`):

| line | cross-line `r` | targets ≥10 sig | unseen ≥10 sig | median sig genes |
|---|---|---|---|---|
| `oikd_2` | 0.016 | **211** | **65** | 9 |
| `oikd_5` | 0.045 | 187 | 28 | 7 |
| `fiaj_1` | 0.046 | 181 | 28 | 7 |
| `zapk_3` | 0.036 | 178 | 26 | 7 |
| `jejf_3` | 0.046 | 172 | 25 | 5 |
| … | | | | |
| `pipw_4` | 0.030 | 56 | 5 | 2 |
| `pipw_5` | 0.014 | 11 | 1 | 2 |
| `fiaj_3` | 0.009 | 6 | 0 | 0 |
| `tolg_4` | 0.014 | **0** | **0** | **0** |

Full table in `outputs/feng_multicontext_v1/reliability.csv`.

**The per-line data is severely underpowered.** Median significant genes per
target is 0–9 out of 5,213, against a median of 189 out of 6,520 in the
publisher's *pooled* table. That is the expected consequence of **median 74
cells per (line, target)** — q10 = 22, only 25.1% of pairs reach 100 cells — and
of **8,241 non-targeting control cells in total (0.7%)**, median 488 per line and
**minimum 86**, so each line's LFCs share a thin, correlated control estimate.

Pooled across lines, 83 of the 141 unseen targets carry ≥10 significant genes
in at least one line.

## 3. Frozen unseen prediction, per line

| line | n | k-NN `r` | k-NN unexplained | ridge `r` |
|---|---|---|---|---|
| `jejf_3` | 141 | **+0.0112** | 1.193 | — |
| `eipl_1` | 141 | +0.0100 | 1.202 | — |
| `oikd_5` | 141 | +0.0098 | 1.221 | — |
| `eipl_3` | 141 | +0.0099 | 1.184 | — |
| `kolf_2` | 141 | +0.0095 | 1.197 | — |
| … median over 19 lines | | **+0.0073** | **1.193** | |
| `tolg_4` | 137 | **−0.0003** | 1.113 | — |
| **pooled** | 141 | **+0.0179** | 2.794 | +0.0118 |

**Zero in every line.** The best single line reaches `r = +0.011`.

## 4. Direct measured transfer, per line

| line | n | measured `r` | unexplained |
|---|---|---|---|
| `oikd_5` | 170 | **+0.1792** | 1.179 |
| `jejf_3` | 169 | +0.1757 | 1.151 |
| `fiaj_1` | 169 | +0.1701 | 1.189 |
| `zapk_3` | 170 | +0.1615 | 1.246 |
| … median over 19 lines | | **+0.1284** | 1.197 |
| `tolg_4` | 167 | +0.0313 | 1.199 |
| **pooled** | 170 | **+0.3201** | 2.500 |

**Positive in 19 of 19 lines**, and larger than the unseen arm in 19 of 19.

> **On the energy metric.** Unexplained fractions exceed 1 throughout, in both
> arms. This is a scale mismatch, not a failure of both methods: Feng's LFCs come
> from their own shrunken model on ~74 cells, while our source deltas are
> differences of log-normalised pseudobulk means with much larger magnitude. The
> ratio `‖truth − pred‖²/‖truth‖²` is dominated by that difference. **Correlation
> is the comparable statistic across the two datasets and is what this report
> reads.** The limitation is recorded in `data/provenance/feng/feng.md` and is
> not removable without the count matrix.

### Is the arm gap just a difference in effect strength?

The measured arm is our essential-gene core; the unseen arm is transcription
factors. Matched on the dataset's own per-target signal count:

| signal quartile (significant genes) | measured `r` (n) | unseen `r` (n) |
|---|---|---|
| 16–53 | **+0.084** (25) | +0.031 (55) |
| 53–78 | **+0.148** (35) | +0.014 (42) |
| 78–188 | **+0.335** (49) | −0.018 (28) |
| 188–3363 | **+0.452** (61) | +0.031 (16) |

**Measured transfer beats prior prediction in every bin, and the gap widens as
the measurement improves.** Direct transfer rises fivefold from the weakest to
the strongest quartile; prior prediction is flat and never leaves the noise.
The gap is a property of the two strategies, not of the two target sets.

## 5. Does performance degrade with context distance?

**Unanswerable from Feng, and that is a finding about the dataset.**

Basal profiles were available after all: the file's `wt_expr` column is exactly
constant within `(line, gene)` across all targets — the line's unperturbed
expression — and it reads no perturbation outcome. But all 19 lines sit at
**0.876–0.895** mean basal similarity to the four training contexts, a span of
**0.019**. They are iPSC lines from 10 donors; as contexts they are nearly
identical.

Spearman across the 19 lines, with bootstrap 95% CIs:

| x | y | Spearman | 95% CI |
|---|---|---|---|
| basal similarity to sources | k-NN `r` | +0.009 | [−0.49, +0.51] |
| basal similarity to sources | measured `r` | −0.239 | [−0.61, +0.28] |
| **reliability floor** | **measured `r`** | **+0.956** | **[+0.83, +0.98]** |
| reliability floor | k-NN `r` | +0.661 | [+0.22, +0.87] |
| n unseen ≥10 sig | measured `r` | +0.677 | [+0.23, +0.98] |

**Between-line variation is explained by measurement reliability, not by
context** — which is exactly predeclared Interpretation 4, *for the between-line
comparison*. Context distance has no detectable relationship, and with a span of
0.019 it could not have.

**This does not undermine the between-arm conclusion.** Both arms are scored on
the same lines with the same reliability, so reliability cannot produce a gap
between them; and the signal-matched analysis in §4 shows the gap survives
conditioning on signal strength.

## 6. Does response conservation across lines explain model success?

**No, and the test is underpowered.**

Median cross-line agreement is **+0.026** for unseen targets and **+0.083** for
measured ones; median line-specific energy is **0.929**. Taken at face value
that says 93% of a response is context-specific — but at 74 cells per
(line, target), "line-specific" and "noise" are not separable, so this number is
an upper bound on genuine context specificity, not a measurement of it.

Spearman(cross-line conservation, frozen-predictor accuracy) = **+0.052** over
141 unseen perturbations. The hypothesis that biological priors capture a
generic conserved program and fail where context-specific response dominates
gets **no support here** — but with per-line reliability this low, the test
could not have detected a moderate effect either.

## 7. Does neighbour agreement remain a useful confidence signal?

**No, on this benchmark.**

| | Spearman vs accuracy |
|---|---|
| arch1 (frozen, previous phase) | **+0.314** |
| Feng, median over 19 lines | **−0.057** |
| Feng, pooled (best-powered) | **−0.064** |
| Feng, support distance (pooled) | +0.056 |

The signal does not transfer. The honest qualification: the quantity it is being
asked to predict — per-perturbation k-NN accuracy — is itself at `r ≈ 0.01` here,
so this is partly a correlation against noise. Neighbour agreement is not
refuted as a concept; it is **unsupported outside arch1**, and it should not be
relied on.

## 8. Does Feng support or refute the arch1 external failure?

**Supports it, and strengthens it.**

arch1 was one external observation, and its weakness was exactly that. Feng is a
second, independent one: different lab, different cell type, different protocol,
19 replicated contexts, 141 unseen perturbations. The frozen predictor fails
there too.

Feng adds what arch1 could not: a **live positive control across 19 contexts**.
On arch1 measured transfer gave `r = +0.302` on 17 perturbations; here it gives
`+0.320` pooled on 170, positive in every line. The two external benchmarks now
agree on both arms.

The standing conclusion is unchanged and is now better supported:

> *STRING-based unseen-perturbation prediction succeeds internally but fails
> externally.*

## 9. Is unseen beta prediction viable under any defensible condition?

**Not on current evidence, and the signal-matched result is why.**

If prior-based prediction were merely *weak*, it should improve when the target
is better measured — as direct transfer does, fivefold across signal quartiles.
It does not move at all. Two external datasets, 241 unseen perturbations
between them, 20 contexts, and no condition has yet been found where it works
outside the internal benchmark.

The one condition **not** yet tested is context proximity: whether a target
context *close* to the training contexts would support it. arch1 and Feng are
both distant, Kaden was unreliable, and Feng's own lines are too uniform to
provide the contrast. That gap is stated in §11 rather than filled by assumption.

## 10. Should the Arc Tier-0 zero-effect policy change?

**No.** The rule was that Tier 0 changes only on external evidence clearly
supporting unseen-perturbation prediction. Feng is external evidence pointing
the other way, in all 19 lines and in the pooled estimate.

Unchanged:

| tier | n | prediction | confidence |
|---|---|---|---|
| **Tier 2** | 7 | direct multi-context transfer, scale-calibrated | raw source agreement |
| **Tier 1** | 79 | direct single-context evidence, strongly shrunk | raw source agreement (weak) |
| **Tier 0** | 214 | **zero perturbation-specific beta** | none — neighbour agreement did not transfer |

One tightening follows from §7: neighbour agreement was carried as a Tier-0
confidence signal on arch1 evidence alone. It did not replicate, so it should be
reported as diagnostic only and not used to gate or weight anything.

No Arc submission has been created.

## 11. What remaining uncertainty matters most?

1. **Context-distance dependence is still untested.** Three external attempts
   and none could measure it: arch1 is a single distant context, Kaden was too
   unreliable, and Feng's 19 lines span 0.019 in basal similarity. Interpretation
   1 — that unseen prediction might be viable *conditional on context support* —
   has not been ruled out and cannot be, on this evidence. Settling it needs a
   dataset with **several contexts at genuinely different distances from the
   training set**, with per-target reliability above roughly 0.3.
2. **The estimator mismatch.** Feng's model-based LFCs and our pseudobulk deltas
   are not the same quantity. It biases the energy metrics badly and correlations
   mildly, in an unknown direction. Resolving it needs the count matrix
   (4.48 GB, deliberately not downloaded).
3. **Per-line power.** At 74 cells per (line, target) the conservation test in
   §6 and the confidence test in §7 are both underpowered; neither null is
   strong.
4. **Whether the measured-transfer arm itself is good enough.** It works
   consistently, but `r = +0.32` pooled is modest in absolute terms, and it is
   the ceiling for Arc Tier 2 — the only tier with real supervision.

## Validation

`uv run pytest` — 423 passed. `ruff check`, `ruff format --check`, `uv build`
clean. All 20 freeze manifests verified.
