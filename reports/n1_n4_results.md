# N1 / N4 results

Date: 2026-10-06/07.

* Protocol: `reports/n1_n4_protocol.md`, SHA-256 `cc2e0f0f…`. Recorded before any N1/N4 data existed
  (`data/provenance/research_v3/n1_n4_protocol_digest.txt`) and **re-hashed identical after the run**
  (`n1_n4_protocol_check.txt`).
* No metric, estimator, grid or pass rule was changed after results existed.
* One implementation deviation (D1, §0) was declared before any N4 output existed.

## Verdicts

| experiment | preregistered verdict | one-line meaning |
|---|---|---|
| **N1** | **FAIL-A**: γ⊥ learnable, but only with a large budget | γ⊥ is learnable, but the k = 20 gain is small: 12–17 % of the 885-anchor gain in 3/4 contexts (39 % in RPE1) |
| **N4** | **WEAK**: incremental information in 4/4 contexts, dominance in only 2/4 | agreement adds information beyond signal strength and reliability, but it does not beat them everywhere |

A further result was not anticipated. The project's frozen v1 estimate of how strongly agreement predicts
transfer (`transferability_confidence_model_v1.md`, Spearman 0.55–0.79) is **inflated by a biased quality
estimator**. With unbiased energies it is 0.29–0.59 (§3.4). The frozen report was not modified.

Neither result is oversold below. N1's small-k γ⊥ gain is statistically real but its lower draw-interval
bounds sit barely above zero (0.001–0.004) in three contexts.

---

## 0. What was run, exactly

**Data object** (`scripts/build_n1_n4_splits.py`):

* For each of the 4 contexts and each of 5 repeats, perturbed cells were split 50/25/25 into F/E1/E2.
* Control cells were split 50/25/25 into CF/CE1/CE2 **before any delta was formed**.
* Disjointness was asserted in code (F4 ✓). The minimum was 15 / 7 / 7 cells per part.
* Control parts:

| context | CF / CE1 / CE2 |
|---|---|
| K562 | 5,345 / 2,672 / 2,672 |
| RPE1 | 5,742 / 2,871 / 2,871 |
| HepG2 | 2,488 / 1,244 / 1,244 |
| Jurkat | 6,006 / 3,003 / 3,003 |

**N1** (`scripts/run_n1_calibration_budget.py`, `src/virtual_cell/analysis/calibration_budget.py`):

* The 4 frozen LOCO folds.
* A fixed test set of 379 perturbations per target; an anchor pool of 885.
* k ∈ {0, 1, 2, 5, 10, 20, 50, 100, 200}, plus the reference k_ref = 885.
* 200 draws per k (5 cell-split repeats × 40 anchor sets).
* Estimators E0, E0s, E1, E2 (template/scale only) and E3, E4 (γ-capable), exactly as in protocol §2.3.
* Anchors entered only as fit-part deltas (F cells vs CF controls). Truth came only from E1/E2 cells vs
  CE1/CE2 controls.
* Runtime: 21.5 min on 4 processes.

**N4** (`scripts/run_n4_agreement_null.py`, `src/virtual_cell/analysis/agreement_null.py`):

* The frozen agreement score.
* Quality q = −D of the frozen predictor B, from the disjoint halves.
* Tests T1 (marginal, vs b1–b6), T2 (partial Spearman) and T3 (exchangeable random-effects null,
  200 simulations).
* Runtime: 2 min.

**Exploratory, not preregistered, labelled as such:**

* `scripts/explore_n1_rank.py`: the rank structure of E4's learned γ⊥.
* `scripts/explore_n4_halves_bias.py`: the v1 quality-estimator bias.

No decision depends on them.

**Deviation D1 (declared before any N4 output).**

* *Protocol §3.2 text:* N4 says to average the target halves over the 5 repeats before forming ⟨h1, h2⟩.
* *Problem:* across repeats the same cells fall into both halves. Averaged halves therefore share noise, and
  ⟨h̄1, h̄2⟩ is biased upward. That contradicts the protocol's own stated purpose (§1: independent noise).
* *Implemented instead:* per-repeat energies, then averaged over repeats. Every other N4 definition is
  unchanged.
* §3.4 shows this choice matters.

**Integrity.**

* All 17 freeze manifests (290 digests) and 11 raw-checksum files verified after the run, with 0 failures
  (`outputs/n1_n4/freeze_verify_after.json`).
* The full test suite passes, including 11 new tests (`tests/test_calibration_budget.py`,
  `tests/test_agreement_null.py`).
* The output manifest is `outputs/n1_n4/manifest_sha256.txt`.
* Code hashes are in `outputs/n1_n4/n1/n1_decision.json`.

---

## 1. N1 failure checks

| check | result |
|---|---|
| F1: non-anchor and canonical target rows replaced by noise → predictions bit-identical | **pass** (test) |
| F2: E0/E0s source-only | **pass** (test) |
| F3a: anchor-permutation null (k = 20, 100; 50 draws) | **no flags.** Permuted γ⊥ M3 medians lie between −0.15 and −0.0003 in all 16 cells, against unpermuted values of +0.02 to +0.30 (`n1_f3a_summary.csv`). The γ⊥ gain needs the correct perturbation ↔ anchor pairing |
| F3b: shared-control diagnostic (k = 5, 20, 100; 50 paired draws) | Sharing one control mean inflates **raw** M0 by +0.004 to +0.024 (e.g. K562 E4 k = 5: −0.027 shared vs −0.051 split). It leaves **M3 identical**, because a control error is constant over perturbations and centring removes it. The γ⊥ endpoint is immune by construction; template-level results would have been mildly inflated without the split |
| F4: disjointness | **pass** (asserted at build) |
| F5: synthetic recovery | **pass.** Planted γ⊥: E4 M3 = 0.92 at k = 100. No γ: E4 −0.001 to −0.011; E3 −0.02 to −0.06 (E3 injects noise when there is nothing to learn) |
| F6: M2 + M3 = M1 identity | **pass** (asserted in every draw) |

---

## 2. N1 results

Notation:

* All metrics are unbiased energy explained on the 379 held-out perturbations, per context.
* Cells give the median over 200 draws, with the 2.5–97.5 % draw interval in brackets for M3.
* k_ref has only 5 draws, so its draw interval is narrow by construction. Use the perturbation bootstrap for
  k_ref.
* Figures: `outputs/n1_n4/figures/n1_M{0,1,2,3,4}_*.png`. Table sources: `outputs/n1_n4/n1/n1_summary.csv`.

### 2.1 Full response, template included (M0): the template is cheap

| context | est | k=0 | k=1 | k=2 | k=5 | k=10 | k=20 | k=50 | k=100 | k=200 | k_ref 885 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| K562 | E0s | +0.046 | | | | | | | | | |
| K562 | E1 | — | −1.088 | −0.457 | −0.006 | +0.180 | +0.261 | +0.313 | +0.330 | +0.339 | +0.347 |
| K562 | E4 | — | −1.088 | −0.520 | −0.032 | +0.186 | +0.293 | +0.413 | +0.474 | +0.517 | +0.589 |
| RPE1 | E0s | +0.289 | | | | | | | | | |
| RPE1 | E1 | — | −0.186 | +0.119 | +0.365 | +0.430 | +0.468 | +0.490 | +0.498 | +0.502 | +0.505 |
| RPE1 | E4 | — | −0.186 | +0.096 | +0.407 | +0.494 | +0.568 | +0.622 | +0.653 | +0.678 | +0.720 |
| HepG2 | E0s | +0.368 | | | | | | | | | |
| HepG2 | E1 | — | −1.176 | −0.433 | +0.090 | +0.250 | +0.332 | +0.382 | +0.402 | +0.411 | +0.421 |
| HepG2 | E4 | — | −1.176 | −0.473 | +0.078 | +0.269 | +0.385 | +0.483 | +0.533 | +0.578 | +0.651 |
| Jurkat | E0s | +0.245 | | | | | | | | | |
| Jurkat | E1 | — | −1.021 | −0.460 | +0.065 | +0.210 | +0.302 | +0.352 | +0.369 | +0.377 | +0.384 |
| Jurkat | E4 | — | −1.021 | −0.530 | +0.024 | +0.207 | +0.326 | +0.419 | +0.471 | +0.518 | +0.599 |

* A target template estimated from anchors (E1) is **harmful below about 5 anchors**. With one anchor, that
  perturbation's own response is mistaken for the template (M0 ≈ −1).
* It overtakes the zero-shot E0s between k = 5 and 10 in K562 and RPE1, and by about k = 20 in HepG2 and
  Jurkat.
* It plateaus by k ≈ 50, at +0.30 (K562) to +0.22 (RPE1) above E0s.
* E2 (also refitting the scale) adds ≤ 0.03.
* So **the template and scale are cheap: about 10–20 random anchors, saturating by about 50.**
* Everything E1/E2 gains here is invisible to M1–M4 by construction (M1 for E1 equals E0s at every k).

### 2.2 γ⊥, the primary endpoint (M3)

γ⊥ is the target response orthogonal to the source consensus. Template, global-scale and per-perturbation
scale estimators score **exactly 0** (E0, E0s, E1, E2 = 0.000 at every k, verified).

| context | est | k=5 | k=10 | k=20 | k=50 | k=100 | k=200 | k_ref 885 |
|---|---|---|---|---|---|---|---|---|
| K562 | E3 | +0.002 [−0.05, +0.03] | +0.017 [−0.09, +0.04] | +0.022 [−0.05, +0.04] | +0.048 [+0.02, +0.07] | +0.077 [+0.06, +0.09] | +0.084 [+0.07, +0.09] | +0.115 |
| K562 | E4 | +0.000 [−0.22, +0.04] | +0.014 [−0.13, +0.09] | +0.051 [+0.001, +0.14] | +0.152 [+0.06, +0.21] | +0.225 [+0.18, +0.25] | +0.283 [+0.25, +0.31] | +0.394 |
| RPE1 | E3 | +0.016 [−0.28, +0.04] | +0.019 [−0.10, +0.04] | +0.033 [+0.01, +0.06] | +0.086 [+0.07, +0.10] | +0.101 [+0.08, +0.11] | +0.114 [+0.11, +0.12] | +0.136 |
| RPE1 | E4 | +0.065 [−0.02, +0.15] | +0.102 [+0.00, +0.18] | +0.171 [+0.04, +0.22] | +0.247 [+0.19, +0.28] | +0.300 [+0.26, +0.33] | +0.347 [+0.32, +0.37] | +0.442 |
| HepG2 | E3 | +0.003 [−0.04, +0.02] | +0.020 [−0.22, +0.04] | +0.020 [−0.08, +0.04] | +0.050 [+0.03, +0.07] | +0.081 [+0.01, +0.09] | +0.076 [+0.06, +0.09] | +0.110 |
| HepG2 | E4 | +0.001 [−0.18, +0.06] | +0.029 [−0.10, +0.09] | +0.071 [+0.004, +0.15] | +0.170 [+0.09, +0.21] | +0.237 [+0.19, +0.27] | +0.304 [+0.27, +0.33] | +0.430 |
| Jurkat | E3 | +0.002 [−0.04, +0.02] | +0.017 [−0.14, +0.03] | +0.025 [−0.06, +0.04] | +0.044 [+0.02, +0.06] | +0.076 [+0.06, +0.09] | +0.084 [+0.07, +0.09] | +0.123 |
| Jurkat | E4 | +0.000 [−0.22, +0.03] | +0.014 [−0.24, +0.06] | +0.046 [+0.003, +0.11] | +0.110 [+0.07, +0.15] | +0.167 [+0.12, +0.20] | +0.229 [+0.19, +0.25] | +0.349 |

k = 1 and 2 are 0.000 for both estimators: with one anchor the residual is absorbed entirely by the template.

**Perturbation bootstrap of M3** (draw-averaged terms, 2,000 resamples of test perturbations):

| context | est | k=10 | k=20 | k=50 | k_ref |
|---|---|---|---|---|---|
| K562 | E3 | +0.011 [+0.007, +0.014] | +0.016 [+0.011, +0.022] | +0.047 [+0.034, +0.059] | +0.116 [+0.096, +0.136] |
| K562 | E4 | +0.011 [−0.002, +0.023] | +0.058 [+0.039, +0.077] | +0.144 [+0.108, +0.175] | +0.394 [+0.344, +0.444] |
| RPE1 | E3 | +0.002 [−0.006, +0.010] | +0.036 [+0.024, +0.048] | +0.086 [+0.070, +0.101] | +0.135 [+0.113, +0.158] |
| RPE1 | E4 | +0.092 [+0.077, +0.107] | +0.164 [+0.142, +0.188] | +0.243 [+0.210, +0.276] | +0.444 [+0.400, +0.484] |
| HepG2 | E3 | −0.000 [−0.007, +0.006] | −0.001 [−0.013, +0.010] | +0.049 [+0.032, +0.065] | +0.109 [+0.083, +0.133] |
| HepG2 | E4 | +0.027 [+0.012, +0.039] | +0.073 [+0.053, +0.092] | +0.164 [+0.129, +0.197] | +0.433 [+0.385, +0.476] |
| Jurkat | E3 | +0.009 [+0.005, +0.012] | +0.017 [+0.012, +0.023] | +0.043 [+0.030, +0.058] | +0.123 [+0.102, +0.146] |
| Jurkat | E4 | −0.000 [−0.009, +0.008] | +0.047 [+0.035, +0.060] | +0.108 [+0.084, +0.132] | +0.349 [+0.303, +0.395] |

**Pass-rule evaluation** (protocol §2.7; `outputs/n1_n4/n1/n1_decision.json`):

| context | E3 C1 / C2 / C3 | E4 C1 / C2 / C3 | E4 M3 at k=20 ÷ at k_ref | γ⊥ onset k (E3, E4) |
|---|---|---|---|---|
| K562 | ✗ / ✓ / ✗ | ✓ / ✓ / ✗ | 0.129 | 50, 20 |
| RPE1 | ✓ / ✓ / ✗ | ✓ / ✓ / **✓** | 0.387 | 20, 10 |
| HepG2 | ✗ / ✓ / ✗ | ✓ / ✓ / ✗ | 0.164 | 50, 20 |
| Jurkat | ✗ / ✓ / ✗ | ✓ / ✓ / ✗ | 0.131 | 50, 20 |

* E4 passes in 1/4 contexts and E3 in 0/4. **N1 does not pass.**
* For both estimators, C2 holds and the k_ref bootstrap lower bound is > 0 in 4/4 contexts. **FAIL-A.**
* Multiplicity: two γ-capable estimators were tested without correction, as predeclared. Neither passes, so
  this does not affect the verdict.

### 2.3 Other components

**Template-removed full response (M1).**

* E0s gives 0.195–0.285.
* E4 at k = 20 gives +0.027 to +0.168 over that, and +0.28 to +0.37 at k_ref.
* E3's gains are about a third of E4's.

**Conserved part (M2, along the source consensus).**

* E0s already explains 0.67–0.77 of the reliable energy along the source direction.
* Refitting the scale (E2) helps RPE1 (+0.13 by k = 5) and HepG2 (+0.05).
* It slightly hurts K562 (−0.01 at k ≤ 10) and Jurkat (−0.03 at k ≤ 10).
* E4 adds +0.01 to +0.03 at k = 20 and +0.08 to +0.13 at k_ref (per-perturbation amplitude learned from
  similar anchors).
* **The conserved component is largely solved zero-shot**; anchors buy amplitude refinements.

**γ_ANOVA (M4).** M4 is the reliable residual of plain source-mean transfer that is removed.

* The global shrink (E0s) already removes 0.14–0.44 of it in K562, HepG2 and Jurkat, but −0.03 in RPE1.
* The γ-specific gain over E2 at the same k (paired by draw, `n1_m4_gain_over_e2.csv`):

| | k = 20 | k_ref |
|---|---|---|
| E4 | +0.028 to +0.163 (lower bounds +0.001 to +0.042) | +0.22 to +0.42 |
| E3 | +0.016 to +0.031 (CIs include 0) | +0.07 to +0.13 |

### 2.4 Exploratory: what E4 learns is low-rank

Not preregistered. One repeat; 10 draws per k (`outputs/n1_n4/n1/exploratory_rank.csv`). E4's γ⊥ prediction
was truncated to its top singular components and M3 recomputed.

| context | k | M3, full | rank 1 | rank 3 | rank 10 | prediction energy in top 1 |
|---|---|---|---|---|---|---|
| K562 | 20 | 0.051 | 0.029 | 0.047 | 0.051 | 0.82 |
| K562 | 885 | 0.392 | 0.189 | 0.318 | 0.386 | 0.39 |
| RPE1 | 20 | 0.165 | 0.130 | 0.148 | 0.151 | 0.63 |
| RPE1 | 885 | 0.451 | 0.270 | 0.347 | 0.405 | 0.54 |
| HepG2 | 20 | 0.089 | 0.067 | 0.087 | 0.081 | 0.62 |
| HepG2 | 885 | 0.428 | 0.182 | 0.323 | 0.413 | 0.32 |
| Jurkat | 20 | 0.054 | 0.025 | 0.051 | 0.050 | 0.58 |
| Jurkat | 885 | 0.346 | 0.118 | 0.251 | 0.334 | 0.26 |


* 90–100 % of E4's γ⊥ gain lives in **≤ 10 target-specific response axes**.
* At k = 20, **one axis** carries 46–79 % of the gain, and three axes carry 90–98 %.
* Small budgets therefore mostly learn *a target-specific direction of response, scaled per perturbation*.
  This is a genuine context × perturbation interaction: neither template nor scale can score here. But it is
  closer to "a target-wide response program" than to many perturbation-specific patterns.
* Larger budgets add further axes.
* This is a hypothesis-generating observation, not a result.

---

## 3. N4 results

Per held-out context; all 1,264 perturbations; quality q = −D with the frozen stability rule.

| | K562 | RPE1 | HepG2 | Jurkat |
|---|---|---|---|---|
| stable fraction | 0.650 | 0.843 | 0.722 | 0.762 |

Tables: `outputs/n1_n4/n4/n4_t{1,2,3}.csv`.

### 3.1 T1: marginal Spearman with q, and the paired-bootstrap difference

The difference is ρ(agreement) − |ρ(b)|, with its 95 % interval.

| statistic | K562 | RPE1 | HepG2 | Jurkat |
|---|---|---|---|---|
| **agreement** | **0.472** | **0.589** | **0.285** | **0.438** |
| b1 source magnitude | 0.399 · diff +0.073 [+0.037, +0.111] | 0.563 · +0.026 [−0.002, +0.053] | 0.301 · −0.016 [−0.059, +0.026] | 0.337 · +0.101 [+0.069, +0.130] |
| b2 source reliability | 0.380 · +0.092 [+0.057, +0.130] | 0.559 · +0.030 [+0.010, +0.050] | **0.344** · −0.059 [−0.092, −0.024] | 0.366 · +0.073 [+0.050, +0.095] |
| b3 source reliable energy | 0.374 · +0.099 [+0.059, +0.139] | **0.575** · +0.014 [−0.011, +0.040] | 0.325 · −0.040 [−0.087, +0.012] | 0.308 · +0.130 [+0.097, +0.165] |
| b4 min source cells | 0.095 · +0.378 | 0.110 · +0.479 | −0.044 · +0.241 | 0.064 · +0.375 |
| b5 target-gene basal | 0.056 · +0.417 | 0.062 · +0.511 | −0.060 · +0.197 | 0.046 · +0.385 |
| b6 noise-implied agreement ceiling | 0.368 · +0.104 [+0.068, +0.141] | 0.545 · +0.044 [+0.027, +0.062] | 0.326 · −0.041 [−0.073, −0.009] | 0.390 · +0.048 [+0.027, +0.071] |

* Agreement dominates b1, b2, b3 and b6 in **K562 and Jurkat**.
* In **RPE1** it ties b1 (magnitude) and b3 (source reliable energy).
* In **HepG2** it is significantly **worse** than source reliability (b2) and the noise ceiling (b6).

### 3.2 T2: partial Spearman, controlling for b1–b5

| | K562 | RPE1 | HepG2 | Jurkat |
|---|---|---|---|---|
| partial ρ | +0.138 [+0.066, +0.218] | +0.217 [+0.148, +0.288] | +0.123 [+0.050, +0.191] | +0.174 [+0.113, +0.245] |
| n | 776 | 969 | 851 | 899 |

Incremental information is positive in **4/4** contexts, but modest: 0.12–0.22.

### 3.3 T3: exchangeable random-effects null (centred responses; 200 simulations)

| context | observed | null median [2.5, 97.5] | class |
|---|---|---|---|
| K562 | 0.651 | 0.666 [0.660, 0.671] | below null |
| RPE1 | 0.624 | 0.322 [0.310, 0.333] | **above** null |
| HepG2 | 0.536 | 0.717 [0.708, 0.725] | below null |
| Jurkat | 0.551 | 0.659 [0.649, 0.668] | below null |

Figure: `outputs/n1_n4/figures/n4_t3_exchangeable_null.png`.

* In 3/4 contexts, agreement predicts transfer **less** well than it would if the contexts were exchangeable
  (K562 only slightly). Real contexts carry non-exchangeable interaction that agreement does not anticipate.
* RPE1 lies far above its null. The plug-in null for RPE1 (median 0.32, against 0.66–0.72 elsewhere) behaves
  differently from the other three. I do not trust this single outlier as evidence that agreement captures
  "transport" information: the null's plug-in mean m_p double-counts source noise and heterogeneity, as
  predeclared.
* **T3 is inconclusive about novelty. It does not support "agreement goes beyond random-effects
  heterogeneity".**

### 3.4 Exploratory: the frozen v1 agreement result is inflated

Not preregistered. `outputs/n1_n4/n4/exploratory_halves_bias.csv`. On the **same** disjoint splits:

| | K562 | RPE1 | HepG2 | Jurkat |
|---|---|---|---|---|
| ρ(agreement, −D), per-repeat unbiased energies (N4) | 0.472 | 0.589 | 0.285 | 0.438 |
| stable fraction | 0.650 | 0.843 | 0.722 | 0.762 |
| ρ, halves averaged over repeats first (v1 construction) | 0.558 | 0.782 | 0.652 | 0.561 |
| stable fraction | 0.996 | 0.994 | 0.998 | 0.996 |
| frozen v1 report (M0, `transferability_confidence_model_v1.md` §3) | 0.554 | 0.786 | 0.656 | 0.556 |

**What this shows.**

* The averaged-halves construction reproduces the frozen v1 numbers to within 0.005.
* Averaging halves over repeats makes ⟨h̄1, h̄2⟩ approximately ‖full mean‖², which includes the noise energy.
  So nearly every perturbation passes the stability rule, including noise-dominated ones, and D is computed
  against a noise-inflated denominator.
* Agreement easily ranks noise-dominated perturbations as poor, which inflates the association.
* The v1 analysis also used a shared control reference. Here the averaging alone reproduces the inflation.

**Consequences for frozen claims** (stated, not edited into frozen reports):

* "Source agreement predicts transferability, Spearman 0.55–0.79" overstates it. The unbiased figure is
  **0.29–0.59**.
* HepG2 changes most: agreement goes from second-best to worse than plain source reliability.
* The v1 risk-coverage and calibration curves were built on the same quality estimator. They should be
  presumed inflated until recomputed.

---

## 4. Scientific interpretation

### The main question: do a few target perturbations improve held-out γ, beyond template and scale?

**Yes, but only modestly at small budgets.** The test separates the cases cleanly:

* Template and scale estimators score exactly zero on γ⊥.
* The anchor-permutation null scores ≤ 0.
* The γ⊥ gain cannot come from shared control noise: F3b leaves it unchanged.

**What the budget buys:**

| budget | E4 (kernel transfer of anchor residuals across perturbations similar in the sources) |
|---|---|
| ≤ 5 random anchors | nothing reliable (RPE1 is a partial exception) |
| 20 random anchors | 5–17 % of the reliable γ⊥ energy; intervals exclude 0 in 4/4 contexts, barely in 3 |
| 885 anchors | 35–44 % |

**Shape of the curve.**

* From about k = 10 the curve rises roughly linearly in log k, with no plateau by 885 anchors (70 % of the
  panel).
* So **γ is learnable from target measurements, and the return per anchor is steady but slow.**
* The preregistered bar (k = 20 reaching ≥ 25 % of the k_ref gain) is met only in RPE1, the deepest and
  cleanest context. Hence FAIL-A.

**Template and scale.**

* They are cheap: 10–20 random anchors capture most of their value, saturating by about 50.
* A single anchor is actively harmful.
* That is the practical floor for any "pilot screen" design.

**What is learned at small k.**

* Exploratory: mostly 1–3 target-specific response axes.
* The cheapest γ is a context-specific program, not perturbation-specific idiosyncrasy.

**What this does *not* show.**

* It does not show that γ is unlearnable at small k by better estimators. E4 is a deliberately simple
  linear smoother on half-depth anchors.
* It does not show that a selection rule could not do better than random. That was not tested and is not
  justified by the preregistered rule.
* It does not generalise beyond these 4 essential-gene screens from 2 labs.

**Relation to the literature.** Molina & Zhang and State report γ learnable at about 30 % of target
perturbations. Our k = 200–885 regime (16–70 % of the panel) agrees. Below about 50 anchors we add the
information that the template and scale are already solved, while γ gains are real but small.

### N4

* Source agreement carries **some** information beyond signal strength and reliability: partial ρ 0.12–0.22
  in 4/4 contexts.
* It does **not** dominate the simple statistics: it beats all four required baselines in only 2/4 contexts,
  and it is worse than source reliability in HepG2.
* It does not exceed exchangeable-heterogeneity expectations in 3/4 contexts.
* Combined with the correction in §3.4, the defensible framing is now:
  **"Agreement is a reasonable, mostly signal-strength-driven trust baseline with a small incremental
  component."**
* It should not be presented as a contribution. It should not be presented as stronger than "source
  reliability", which is equally simple.

---

## 5. Recommendation for the next experiment

**Not N2.** By the predeclared mapping (protocol §4, FAIL-A), acquisition or selection development is not
justified. I have not started it.

### Recommended: N3, an unchanged-protocol replication of N1 on the 6-context design

Add X-Atlas HCT116 and HEK293T (about 1,062 shared perturbations), with the N1 protocol frozen as-is: same
estimators, metrics, grid and pass rule.

**Why this, first.**

1. Every N1 statement rests on 4 essential-gene screens from 2 labs. The log-linear γ⊥ curve, "template
   saturates at about 50" and the FAIL-A boundary all need replication before they are findings.
2. Three more source contexts per fold (5 instead of 3) directly tests whether more source diversity raises
   the small-k γ⊥ gain. That is the cheapest lever not yet pulled, and it decides whether FAIL-A is a property
   of γ or of having only 3 sources.
3. It also re-runs N4 on 6 folds, which decides whether "WEAK" is stable.

**Prerequisites.**

* A written research-use license memo for X-Atlas (CC BY-NC-SA). Non-commercial academic use only, separate
  from all competition code paths (already enforced in code).
* A new predeclaration that copies the N1/N4 protocol verbatim, changing only the context list.

**Cost.** CPU only. Streaming about 2 × several-million FiCS cells is probably 1–3 h, plus about 40 min of
compute.

**What would change the plan:**

* **Small-k γ⊥ gain clearly larger with 5 sources** (C3 met in ≥ 4/6): revisit P1, and only then consider a
  predeclared N2.
* **Same FAIL-A pattern in 6/6:** P1 becomes a limits paper. "Template and scale need about 20 anchors; γ⊥
  grows log-linearly and slowly; few-shot γ is dominated by a handful of target programs."
* **The X-Atlas folds behave differently:** lab/assay shift is entangled. The backup direction P2 (N5 ladder)
  moves up.

### Separate housekeeping (not an experiment)

The frozen v1 transferability claims (README finding 5, `project_synthesis.md` §5, the paper outline §7)
quote the inflated 0.55–0.79. Correct them with a dated erratum that points to §3.4. Do not edit the frozen
reports. I have not made this change; it is your call.
