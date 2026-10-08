# N3 preregistered protocol: six-context replication, source-count ablation, rank structure

Written 2026-10-07. Status at the time of writing:

* No X-Atlas expression value has been read.
* No six-context data object exists.
* No N3 result exists.

Only identifier metadata was read (perturbation labels, gene symbols, counts of cells, batches,
`pass_guide_filter`). The SHA-256 of this file is recorded in
`data/provenance/research_v3/n3_protocol_digest.txt` before any data object is built.

After that, nothing below may change. Any deviation found while running is reported as a deviation, with its
consequence.

Parents:

* `reports/n1_n4_protocol.md` (N1/N4; incorporated by reference where stated, with deviation D1 adopted);
* `reports/n1_n4_results.md` (N1 FAIL-A, N4 WEAK);
* `reports/xatlas_license_memo.md` (research-only use; no competition use).

**Hard constraints:**

* No model tuning to rescue N1/N4.
* No acquisition or selection model.
* No source subset selected using target outcomes.
* No frozen result modified.

---

## 0. Six-context design

### Contexts (index order)

| index | context | source |
|---|---|---|
| 0 | K562 | scPertEval `replogle22k562` |
| 1 | RPE1 | scPertEval `replogle22rpe1` |
| 2 | HepG2 | scPertEval `nadig25hepg2` |
| 3 | Jurkat | scPertEval `nadig25jurkat` |
| 4 | HCT116 | X-Atlas/Orion |
| 5 | HEK293T | X-Atlas/Orion |

### Axes

Both are identifier-only and written to `data/splits/six_context_n3/` with hashes.

* **Perturbations: 1,062.** Members of `four_context_v1/shared_perturbations.txt` with ≥ 30 cells labelled with
  that `gene_target` in **both** X-Atlas lines. Order follows the four-context file.
* **Genes: 6,499.** Members of `four_context_v1/shared_genes.txt` present as a unique `gene_id` symbol in both
  X-Atlas lines. Order follows the four-context file.

### scPertEval contexts (0–3)

* Canonical deltas: the frozen `delta_tensor.npy`, subset to the N3 axes.
* Split-part means: the N1 object `outputs/n1_n4/splits/` (same seeds, same cells), subset to the N3 axes.
* Subsetting a mean of per-gene values is exact.

### X-Atlas contexts (4, 5)

* **Cells:** every cell whose `gene_target` is one of the 1,062 perturbations, plus every `Non-Targeting` cell,
  across all batches (all have `pass_guide_filter = True`).
* **Per-cell value:** `log1p(1e4 · count / total_counts)`, with `total_counts` from the cells table (all 38,606
  features), restricted to the 6,499 genes. Pseudobulk is the mean over cells (the canonical `mean_log`
  aggregation).
* **Canonical delta:** mean over all of the perturbation's cells − mean over all `Non-Targeting` cells.
* **Disjoint splits:** the N1 algorithm (`calibration_budget.disjoint_split_means` logic).
  * Within each perturbation, cells are ordered by `cell_integer_id`, then shuffled.
  * Parts are F = ⌊n/2⌋, E1 = ⌊n/4⌋, E2 = ⌊n/4⌋.
  * Controls are split CF/CE1/CE2 the same way, before any delta exists.
  * 5 repeats, seed `[20261006, 4 or 5, r]`.

### Data integrity gates (stop if any fails)

* **G1:** non-zero entries read for the selected cells ≥ 99 % of Σ `n_genes_by_counts` over those cells.
* **G2:** the canonical delta equals the cell-count-weighted combination of the part means, to 1e-6 relative.
* **G3:** on-target knockdown. Among perturbations whose target gene is on the 6,499-gene axis, the median
  canonical delta of the target gene is < 0 in each X-Atlas line.
* **G4:** part disjointness asserted.

### Assay incompatibilities (stated; never silently pooled)

The X-Atlas lines differ from contexts 0–3 in:

* lab;
* chemistry (fixed-cell FiCS, 10x Flex-type);
* sgRNA library (genome-wide vs essential-wide);
* sequencing depth (median 17–21k UMI per cell);
* normalisation basis (all genes vs scPertEval's filtered set);
* batch structure (controls pooled over 109/223 batches);
* knockdown strength (published median 75 % HCT116, 52 % HEK293T).

**Every metric is computed and reported per held-out context. No metric is averaged across contexts for any
decision.** A context's responses enter other folds only as a source.

### Assay-compatibility diagnostics (descriptive, §6)

Per context, on the N3 axes:

* median ‖δ‖;
* median per-perturbation split-half reliability;
* median on-target delta;
* the 6 × 6 basal-similarity matrix (control means);
* the 6 × 6 matrix of median per-perturbation Pearson between centred canonical responses.

---

## 1. N3-A: six-context replication of N1

This is N1 protocol §2, verbatim, except for the substitutions below.

* **Folds:** 6 LOCO folds. The target is each context; the sources are the other **5**.
* **Test set:** 319 perturbations (= round(0.3 × 1,062)), seed `[20261006, 11, t]`. The pool is the other 743.
* **Budgets:** k ∈ {0, 1, 2, 5, 10, 20, 50, 100, 200}, plus k_ref = 743 (the full pool). All are supported.
* **Draws:** 200 per k (5 repeats × 40), anchors seeded `[20261006, 12, t, k, r, d]`. k = 0 and k_ref have 1 per
  repeat.
* **Estimators:** E0, E0s, E1, E2, E3, E4 exactly as in N1 §2.3, with the 5 sources.
  * E3's prior is `s/5` per source.
  * E4's features are `concat_s(S̃_s/‖S̃_s‖)/√5`.
  * `s = fit_scale(D, sources)` over the 5 sources.
  * λ grids, the k < 5 rule and LOO selection are unchanged.
* **Metrics:** M0, M1, M1r, M2, M3 (γ⊥ relative to the 5-source consensus), M3r and M4, exactly as in N1 §2.4–2.5.
* **Uncertainty:** draw percentiles; perturbation bootstrap with 2,000 resamples, seed `[20261006, 13, t, k]`.
* **Failure checks:**
  * F1/F2 as tests (with 5 sources).
  * F3a anchor-permutation null: k ∈ {20, 100}, 50 draws, seed `[20261006, 16, …]`.
  * F3b shared-control diagnostic: k ∈ {5, 20, 100}, 50 paired draws.
  * F4 and F6 as in N1.

### 1.1 Per-context classification

This is the N1 rule (§2.7), applied per target with k_ref = 743.

* **PASS:** E3 or E4 meets C1 ∧ C2 ∧ C3.
* **FAIL-A-type:** not PASS, and for E3 or E4: C2 holds and the k_ref perturbation-bootstrap lower bound of M3
  is > 0.
* **FAIL-B-type:** otherwise.

### 1.2 Pattern components (descriptive, per context)

These mirror the three claims under replication.

* **k_T50.** The smallest grid k ≥ 1 at which the median M0 gain of E2 over E0s reaches 50 % of its k_ref
  gain. If the k_ref gain is ≤ 0, report "no template/scale gain".
  * P1 "template/scale learnable with small k" holds iff k_T50 ≤ 20.
* **k_γ50.** The smallest grid k at which E4's median M3 is ≥ 50 % of its median M3 at k_ref. If M3(k_ref) ≤ 0,
  report "none".
  * P2 "γ slower" holds iff k_γ50 > k_T50.
* **P3 "low budget recovers a minority".** E4's median M3(k = 20) / median M3(k_ref) < 0.25.

### 1.3 Primary replication decision (R-A1)

**FAIL-A replicates externally** iff both HCT116 and HEK293T are FAIL-A-type.

* If either is PASS: "small-k γ works externally". FAIL-A did not replicate.
* If both are FAIL-B-type: "γ⊥ not learnable externally". FAIL-A did not replicate.

The four original contexts are re-run with 5 sources. They are **not** independent replications; they feed
N3-B.

---

## 2. N4 replication on six folds

N1/N4 protocol §3, with deviation D1 (per-repeat energies), and these substitutions:

* The agreement is the mean of the 10 pairwise source Pearsons.
* b6 is averaged over the 10 pairs.
* T3 uses `τ²_p = max(0, Σ_s‖S̃_s,p − m_p‖² / ((5 − 1)·G) − mean_s σ²_s,p)`, with 200 simulations, seed
  `[20261006, 15, t, i]`.
* The bootstrap seed is `[20261006, 14, t, …]`.
* **Stability threshold:** the frozen 1.96738 × (6,499 / 6,640) = **1.92563**. Signal energy scales with the
  number of genes.

**Pass rule.**

* PASS iff (i) and (ii) of N1/N4 §3.5 hold in **≥ 5 of 6** contexts. This is the N1/N4 proportion of 3/4,
  rounded up.
* WEAK iff (i) holds in ≥ 5 of 6.
* FAIL otherwise.
* Per-context results are reported regardless.

---

## 3. N3-B: source-count ablation, R_γ(k, m)

**Units.** For each target t (6), every non-empty subset S of its 5 sources: 31 subsets.

| m | subsets |
|---|---|
| 1 | 5 |
| 2 | 10 |
| 3 | 10 |
| 4 | 5 |
| 5 | 1 |

All subsets are enumerated. **None is selected.**

**Fixed per target, shared across all subsets (paired design):**

* the test set (the N3-A T_t);
* the anchor draws: the first 8 draws per repeat of N3-A (d = 0..7), i.e. 40 draws per k, with identical seeds;
* the budgets k ∈ {0, 5, 10, 20, 50, 100, 200, k_ref}.

**Estimators with source set S:**

* E0s (k = 0), with `s = fit_scale(D, S)`; for m = 1, `fit_scale` returns 1.0.
* E2.
* **E4 (the primary calibration estimator).** Its features use only S.

**Metrics:**

* **M3_ref (primary for R_γ):** γ⊥ relative to the **fixed** consensus direction of all 5 sources, identical for
  every m. The evaluation target therefore never changes with m. Predictors built from fewer sources may then
  score non-zero, or negative, at k = 0.
* **M3_own (secondary):** γ⊥ relative to S's own consensus.
* **M1:** template-removed full response (definition-free across m).
* M0 for E2.

**Primary object.**

* `R_γ(k, m)` = for each draw, the mean over subsets of size m of E4's M3_ref; then the median over the 40 draws,
  with the 2.5–97.5 % draw interval.
* `R_full(k, m)` is the same with M1.
* For k = 0 the estimator is E0s.

**Consistency check F-B1.** For m = 5, the draw-level values equal N3-A's E4 values on the same draws (to
1e-10).

### 3.1 Q1. At fixed k, does m materially help?

* `Δ(k) = R_γ(k, 5) − R_γ(k, 3)`, paired per draw (3 = N1's source count), and also `R_γ(k, 5) − R_γ(k, 1)`.
* **"Material" at k = 20** iff Δ(20) ≥ 0.03, its draw lower bound is > 0, and Δ(20) ≥ 0.25 · R_γ(20, 3)
  (the last required only if R_γ(20, 3) > 0).

### 3.2 Q2. At a fixed γ target, how much budget does m save?

* Target level `g*_t` = median R_γ(50, 3).
* `k*_m` = the smallest budget at which median R_γ(k, m) ≥ g*. It is log-linearly interpolated between adjacent
  grid points k ≥ 5. Report "> k_ref" if never reached.
* Reduction factor `ρ_t = k*_3 / k*_5`. `k*_m` is reported for all m.
* **Material** iff ρ_t ≥ 1.5.

### 3.3 Outcome-B criterion

Q1 material in ≥ 4/6 targets **and** ρ_t ≥ 1.5 in ≥ 4/6 targets.

### 3.4 Q3. Averaging/reliability vs complementary contexts

Source-only subset descriptors, computed without any target response:

| descriptor | definition |
|---|---|
| **Rel_S** | mean over s ∈ S of the median over perturbations of source split-half reliability (Pearson of F vs (E1+E2)/2, averaged over repeats) |
| **Div_S** | mean over pairs in S of the disattenuated dissimilarity `1 − r_ab / √(ρ_a ρ_b)`. `r_ab` = median over perturbations of the Pearson between the panel-centred canonical responses of a and b; `ρ` = the Spearman–Brown full-depth version of Rel. Defined for m ≥ 2 |
| **Partner_S** | max over s ∈ S of the basal similarity (Pearson of control means) to the target. Target **controls** are legal |

**Response.** `R_S` = E4's median M3_ref over draws at **k = 50** (primary); k = 20 is secondary.

**Test.**

* Blocks are (t, m) for m ∈ {2, 3, 4}: 25 subsets per target, 150 points.
* Within each block, rank-transform R, Rel, Div and Partner. Pool the within-block ranks.
* Compute the pooled partial Spearman of R with Div given Rel, of R with Rel given Div, and (descriptively) of R
  with Partner given Rel and Div.
* p-values from 2,000 within-block permutations of R, seed `[20261006, 17]`.

**Reading:**

* **"complementarity"** iff partial(R, Div | Rel) > 0 with p < 0.05;
* **"averaging/reliability"** iff partial(R, Rel | Div) > 0 with p < 0.05;
* both, either or neither may hold.

**Combination check.** At k = 50: R_γ(k, 5) − max over single sources of R_{s}(k), paired per draw with the draw
interval. This asks whether the full set beats the best single source.

---

## 4. N3-C: rank structure (fixed procedure, written before any X-Atlas rank result)

For each target (primary interpretation on HCT116 and HEK293T; the 4 originals are reported but are not
independent of the N1 exploratory finding), with E4 and 5 sources:

* **Draws:** k ∈ {20, 50, 200, k_ref}. 20 draws for k < k_ref (the first 4 per repeat of N3-A); the 5 repeats
  for k_ref.
* **Estimator side.** `P⊥` = E4's test-set prediction, centred over T and made orthogonal to u. Its SVD is
  `P⊥ = U S Vᵀ`.
  * `f(r) = median_d M3(P⊥_r) / median_d M3(P⊥)`, with `P⊥_r` = the top-r reconstruction, re-orthogonalised to u,
    for r ∈ {1, 2, 3, 5, 10}. **f(1) is "the fraction of calibration gain explained by the first axis".**
  * The cumulative gain of the first r axes is `median_d M3(P⊥_r)`.
  * `r_eff90` = the smallest r ∈ {1..50} with f(r) ≥ 0.9 (otherwise "> 50").
  * Participation ratio `(Σ s²)² / Σ s⁴` (median over draws).
  * f is undefined, and reported as "no gain", if median M3(P⊥) ≤ 0.01.
* **Axis stability.** For all pairs of draws at the same k:
  * `|cos|` between the first gene-space axes `v1`;
  * the top-3 subspace overlap `‖V3aᵀ V3b‖²_F / 3`.
  * Medians are reported. Null: the same statistics from 10 **permuted-anchor** E4 fits per k
    (F3a-style permutation, seed `[20261006, 18, t, k, d]`).
* **Data side** (estimator-independent). V_r = the top-r right singular vectors of the target's **pool**
  perturbations' fit-part γ⊥ (F cells vs CF controls, centred over the pool, orthogonalised to each
  perturbation's own consensus direction). Evaluation-only, and disjoint from T and from E-cells.
  * `o(r) = Σ_p ⟨V_rV_rᵀ e1⊥_p, V_rV_rᵀ e2⊥_p⟩ / Σ_p ⟨e1⊥_p, e2⊥_p⟩` on T: the unbiased fraction of reliable
    test γ⊥ energy inside the top-r pool subspace.
  * r ∈ {1, 3, 10, 30}, repeat 0.
  * This separates "γ⊥ is low-rank in the data" from "E4 makes it low-rank": at k anchors, rank(P⊥) ≤ k by
    construction.

**Hierarchy criteria, per context:**

| id | criterion |
|---|---|
| h1 | template/scale early: k_T50 ≤ 20 (§1.2) |
| h2 | one axis dominates early: f(1) at k = 20 ≥ 0.40 |
| h3 | the first-axis share falls with budget: f(1) at k_ref < f(1) at k = 20 |
| h4 | effective rank grows: r_eff90(k_ref) > r_eff90(k = 20) |

The hierarchy template → scale → low-dimensional program(s) → perturbation-specific γ "reproduces" in a context
iff h1–h4 all hold.

**Outcome-D criterion.** h2 holds in **both** X-Atlas contexts, **and** in both the median top-1 axis stability
at k = 20 exceeds the 97.5th percentile of its permuted-anchor null.

**Novelty discipline.** Low-rank perturbation-response structure and gene programs are established. Examples:

* Replogle 2022 programs;
* Systema's systematic variation;
* low-rank linear baselines (Ahlmann-Eltze 2025);
* the C6 spectrum in this repo.

N3-C results are labelled "replication of a structural observation", never "novel", until compared against that
literature in the results.

---

## 5. Leakage controls (in addition to N1 F1/F2/F4)

* Target perturbation responses are never used to build the axes, subsets, descriptors, hyperparameters, test
  sets or anchor draws.
* Target controls are used only through the disjoint CF part (fitting) and CE parts (truth), and, for
  Partner_S, through the control mean. That is basal information, which is legal.
* Subset enumeration is exhaustive and fixed.
* Rel_S, Div_S and Partner_S read no target response.
* **Test:** the noise-replacement test from N1, run with 5 sources and with a 2-source subset.

---

## 6. Interpretation rules (fixed now)

| outcome | criterion |
|---|---|
| **A** | R-A1 holds (FAIL-A replicates in both X-Atlas contexts), and the Outcome-B criterion fails |
| **B** | the Outcome-B criterion (§3.3) holds |
| **C** | any of: (c1) an X-Atlas context's classification differs from FAIL-A-type while ≥ 3 of the 4 originals are FAIL-A-type in N3-A; (c2) either X-Atlas context's median E4 M3(k_ref) is < 0.5 × the minimum, or > 2 × the maximum, of the four originals; (c3) either X-Atlas context's zero-shot E0s M1 (k = 0) is ≤ 0 while all four originals are > 0 |
| **D** | §4 Outcome-D criterion |

* Every outcome whose criterion holds is reported.
* If C holds, any A/B/D statement is qualified as possibly assay- or lab-specific.

**Decision on the calibration-publication direction:**

* **Retain** (small-budget calibration as a practical solution) only if an X-Atlas context is PASS, or B
  holds.
* **Revise** (to a sample-complexity / limits paper, plus structure if D) if A holds.
* **Kill** if both X-Atlas contexts are FAIL-B-type and B fails; or if C holds via (c1) or (c2) and no
  component of the N1 pattern replicates in either X-Atlas context (P1–P3 all fail).
* Otherwise, report "unresolved" with the reason.

## 7. Outputs

* `outputs/n3/` (git-ignored): data object, draws, summaries, decisions, and a SHA-256 manifest.
* Aggregate figure sources only under `data/figure_sources/n3/` (license memo §7).
* `reports/n3_results.md`, carrying the X-Atlas attribution block.
* **Stop after N3. No acquisition model.**
