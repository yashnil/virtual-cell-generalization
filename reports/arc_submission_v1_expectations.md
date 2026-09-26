# Arc submission candidate v1: pre-result expectations and interpretation rules

Written 2026-09-25, **before any hidden score exists**, and committed before
submission. Its purpose is to let the hidden A/B/C result *test* our reasoning
rather than let us rewrite the reasoning afterwards. These are qualitative
hypotheses, not score predictions.

Candidate: `arc_count_space_baseline_v1`
(`reports/arc_submission_v1_manifest.md`), validation partition
`vcc2026-val-1`.

## What we can and cannot say numerically

The official score rescales each member as `s = (u − b) / (r − b)`. Here `b`
is the context mean-response baseline (0) and `r` is a replicate anchor (1).
The published `b` and `r` were measured on Arc's own data. Our public
count-space benchmark (held-out K562) produced **raw** values only, and
applying Arc's anchors to them would be an invention
(`arc_count_space_baseline_v1.md` §I). So **no normalized score is predicted
here.**

Two reference points *are* frozen and can be used as-is:

- **Control-emitting submission** (`outputs/arc_bridge_v1/score_accounting.csv`):
  overall **−0.311**. Per member:

  | member | control-emitting score |
  |---|---|
  | pds | 0.000 |
  | expression | 0.000 |
  | DE direction fidelity | **−1.712** |
  | DE direction reach | −0.080 |
  | DE significance overlap | −0.078 |
  | DE log-FC | +0.002 |

- **Mixed-panel discrimination on public data:** `pds_cosine` 0.535 at Arc's
  7 / 79 / 214 tier prevalence, against the constant-prediction floor of 0.500
  (pseudobulk, public held-out contexts).

The structural fact that dominates every expectation below: **214 of 300
targets (71%) receive `m_hat` alone.** Its norm is about 0.07 because the
estimator shrank to 0.119, so those targets are **close to control-emitting**.

## K. Expectations, per member

| member | expectation | reasoning |
|---|---|---|
| **PDS** (perturbation discrimination) | The **strongest relative** member, but a **small** absolute gain above baseline. | Tier 1/2 direct effects and G1 transport produced the whole public discrimination gain. The 214 identical Tier-0 predictions tie with each other, so on a mixed panel most of that gain disappears (public 0.535 vs 0.500 floor). |
| **Expression accuracy** | **Among the weakest; near or below 0.** | On public data no generator beat control resampling on expression error (G1 1.062 vs 1.036). Tier 0 contributes nothing above control. |
| **DE direction fidelity** | Some signal for supported perturbations, but the member as a whole is expected **below baseline, possibly well below**. | It is a *yield* member: calling nothing scores −1.71. The 71% of near-control Tier-0 targets should call few DE genes. |
| **DE direction reach** | **Limited.** | Shrinkage (w = 0.50 / 0.25, m_hat scalar 0.119) makes predicted effects small, so few genes reach significance. |
| **DE significance overlap** | **Uncertain.** | It depends jointly on predicted effect magnitude and on the cell-to-cell variation G1 generates. The public result (0.122 raw vs 0.053 for control resampling) does not transfer via anchors. |
| **DE log-FC accuracy** | **Limited, likely near 0.** | Bounded by how well the context main effect and effect magnitudes are calibrated, which is the weakest part of the model on Arc's sources. |
| **Overall** | Not predicted numerically. The dominant structural risk is that direction fidelity pulls the overall **toward or below the control-emitting reference (−0.311)**. A positive overall would mean the supported targets carry more than the public evidence implies. | — |

## L. Interpretation rules, predeclared

Scores below are the official **normalized** per-member values.

**Bands** (fixed now):

| band | normalized score |
|---|---|
| clearly positive | s ≥ +0.05 |
| near baseline | −0.05 < s < +0.05 |
| below baseline | s ≤ −0.05 |

The control-emitting reference for each member is read from
`score_accounting.csv`. It is not re-derived after the result.

| pattern | observed | interpretation | what it licenses next |
|---|---|---|---|
| **1** | PDS clearly positive; expression near or below baseline | Direction and discrimination exist; magnitude and the context baseline are the main problem. | Work on the context main effect / magnitude. Not on perturbation representation. |
| **2** | PDS and DE direction fidelity both near or below baseline | Supported perturbation effects themselves fail to transfer to Arc's contexts. | Question direct transfer across distant contexts before anything else. No added capacity. |
| **3** | Direction fidelity above its control reference, but reach below baseline | Directions are often right; effects are too weak or over-shrunk. | Magnitude calibration becomes a legitimate predeclared question for v2. |
| **4** | Reach and fidelity reasonable (≥ baseline), significance overlap below baseline | Generator variance or per-gene effect allocation is wrong. | Examine G1's cell-level variance, not the mean model. |
| **5** | All six near or below baseline | Public validation materially failed to predict hidden transfer. | **Do not** increase model complexity. Re-examine the public-to-Arc evidence chain first. |
| **6** | Discrimination-related members (PDS, and fidelity/reach relative to their control references) show tier-supported signal while expression stays weak | Tier-supported biology is useful; the unsolved part is the context baseline. | Future work prioritises the context baseline / magnitude, not perturbation representation. |

**Additional predeclared readings:**

- If direction fidelity lands near **−1.7**, the near-control Tier-0 majority
  dominates it. That is an expected consequence of the frozen Tier-0 policy,
  **not** evidence that the supported predictions are wrong. Judge supported
  targets by PDS and the other members instead.
- A single validation score is **one** external measurement on three contexts.
  It is not used to tune weights, `m_hat`, or G1 on A/B/C. Any change it
  motivates must be predeclared and selected on public held-out contexts,
  exactly as v1 was.
- The final D/E/F contexts are different cell lines. A good or bad A/B/C score
  does not transfer automatically.
