# N1 / N4 preregistered protocol

Written 2026-10-06, **before** any N1/N4 data object was built and before any N1/N4
result existed. Its SHA-256 is recorded in
`data/provenance/research_v3/n1_n4_protocol_digest.txt` at the moment of writing.

After that, no metric, estimator, grid, threshold or pass rule here may change.

* Implementation details not fixed here (for example vectorisation) may be chosen freely, but they must
  implement exactly these definitions.
* Any deviation discovered while running is reported in the results as a deviation, together with its
  consequence. It is never silently absorbed.

Source of the design: `reports/publication_novelty_audit.md` §8.1, §11 (N1, N4).

Scope:

* **Run:** N1 and N4 only.
* **Not run:** acquisition functions or selection rules (N2), the X-Atlas 6-context design (N3), or the
  shift ladder (N5).
* **Not modified:** any frozen result, frozen module or frozen script. All new code is new files.

---

## 0. Frozen inputs

| input | path | role |
|---|---|---|
| balanced design | `data/splits/four_context_v1/` (4 contexts × 1,264 perturbations × 6,640 genes) | perturbation and gene axes |
| LOCO folds | `data/splits/loco_v1/` | 4 outer folds: each context held out once, the other 3 are sources |
| canonical response tensor | `data/processed/four_context_v1/delta_tensor.npy` | **source-context** responses only (full cells, full controls) |
| canonical control means | `data/processed/four_context_v1/control_means.npy` | target basal expression (N4 baseline b5); shared-control diagnostic (N1 F3b) |
| cell counts | `data/processed/four_context_v1/cell_counts.npy` | N4 baseline b4 |
| raw cells | `data/raw/scperteval/*_processed_complete.h5ad` (hash-verified) | new disjoint cell splits (§1) |
| frozen predictor B | `virtual_cell.modelling.pathway_residual.fit_scale` / `baseline` | E0s; the N4 transfer predictor |
| frozen agreement | `virtual_cell.modelling.transferability.source_agreement` (mean pairwise gene-level Pearson among the 3 raw source responses) | the N4 score under test |
| frozen stability rule | `outputs/transferability_v1/stability_rule.json`, min signal energy 1.96738 | N4 D masking |

All existing freeze manifests must verify before and after the run.

---

## 1. New data object: fully disjoint cell splits (built once, used by N1 and N4)

For **each** of the 4 contexts and each repeat r ∈ {0..4} (R = 5), with seed `[20261006, context_index, r]`:

**Perturbed cells** of perturbation p (n_p ≥ 30 cells), randomly permuted:

| part | size | use |
|---|---|---|
| **F** (fit) | the first ⌊n_p/2⌋ | fitting |
| **E1** | the next ⌊n_p/4⌋ | evaluation |
| **E2** | the next ⌊n_p/4⌋ | evaluation |
| (dropped) | any remainder | unused |

**Control cells** (m cells), randomly permuted, **before any target-derived quantity is computed**:

| part | size |
|---|---|
| **CF** | the first ⌊m/2⌋ |
| **CE1** | the next ⌊m/4⌋ |
| **CE2** | the next ⌊m/4⌋ |

The three parts are disjoint (asserted in code).

**Deltas** (log1p CP10K, mean over cells, exactly the canonical aggregation `mean_log`):

```
fit[r,c,p] = mean(F cells) - mean(CF cells)
e1[r,c,p]  = mean(E1 cells) - mean(CE1 cells)
e2[r,c,p]  = mean(E2 cells) - mean(CE2 cells)
```

The perturbed-cell part means and the control part means are stored separately, so that the shared-control
diagnostic (F3b) can be formed. Output: `outputs/n1_n4/splits/` (git-ignored), with SHA-256 manifest.

**Consequences, stated in advance:**

* The fit side uses half the target cells. Anchor measurements are therefore noisier than a full-depth pilot
  screen, which makes the N1 test conservative.
* Evaluation uses two independent quarter-depth copies, so energies are estimated unbiasedly (§2.5).

---

## 2. N1: calibration-budget learning curve

### 2.1 Split

* **Outer:** the 4 LOCO folds. The held-out context is the *target*; the other 3 are sources.
* **Test set T:** per target context, a fixed random 30 % of the 1,264 perturbations (379), drawn once with
  seed `[20261006, 1, context_index]`. T is identical for every k, estimator and repeat.
* **Anchor pool:** the remaining 885 perturbations.
* **Budgets:** k ∈ {0, 1, 2, 5, 10, 20, 50, 100, 200}, plus a **reference point k_ref = 885** (the whole
  pool). k_ref is reported as the large-budget reference used in the pass rule; it is not part of the
  requested grid.
* **Draws:** for each k ≥ 1, each repeat r has 40 random anchor sets K ⊂ pool, |K| = k, drawn uniformly
  without replacement with seed `[20261006, 2, context_index, k, r, d]`. That gives **200 draws per
  (context, k)**. k = 0 and k_ref have 1 draw per repeat (5 values).

### 2.2 Information legally available to a prediction

1. Source responses `S_s` (canonical tensor, 3 source contexts, all 1,264 perturbations).
2. Perturbation identity.
3. For the target, **only** `fit[r, target, K]`, i.e. the anchors' fit-part deltas (fit cells against CF
   controls).

Nothing else from the target is read by any estimator:

* no evaluation-part cell;
* no CE1/CE2 control cell;
* no non-anchor fit delta;
* no canonical target delta;
* no four-context β/γ.

Hyperparameters are chosen only from the anchors (§2.3).

### 2.3 Estimators (fixed ladder)

Notation:

* `A = mean_s S_s` (P × G).
* `t_A = mean over all 1,264 perturbations of A` (the source template).
* `s = fit_scale(D, sources)` (frozen; source-only).
* `S̃_s = S_s − mean_p S_s` (each source centred over the full panel).
* `Y_K = fit[r, target, K]`.

| id | estimator | uses anchors? | can predict γ⊥? |
|---|---|---|---|
| **E0** | `P = A` | no | no |
| **E0s** | `P = t_A + s (A − t_A)` (frozen predictor B) | no | no |
| **E1** | `P = E0s + α̂`, with `α̂ = mean_{a∈K} (Y_a − E0s_a)` (target template from anchors) | k ≥ 1 | no |
| **E2** | `P = c + s' (A − t_A)`, with `s' = Σ_a ⟨Y_a − Ȳ_K, Ã_a − Ã_K⟩ / Σ_a ‖Ã_a − Ã_K‖²` (anchor-centred least squares, `Ã = A − t_A`) and `c = Ȳ_K − s' · mean_K Ã`. For k = 1, `s' := s` (reduces to E1) | k ≥ 1 | no |
| **E3** | gene-wise source reweighting (below) | k ≥ 1 | **yes** |
| **E4** | E2 + kernel ridge of anchor residuals over perturbation similarity (below) | k ≥ 1 | **yes** |

**E3.** For each gene g separately:

* Model: `Y_{a,g} = w0_g + Σ_s w_{s,g} S̃_{s,a,g}`.
* Ridge penalty `λ_g ‖w_g − w_prior‖²` with `w_prior = (s/3, s/3, s/3)`. The intercept is unpenalised.
* `λ_g = λ_rel × mean diagonal of the anchor-centred Gram matrix of gene g`.
* `λ_rel ∈ {0.01, 0.1, 1, 10, 100}`, **one λ_rel for all genes**, chosen by exact leave-one-anchor-out squared
  error pooled over genes and anchors.
* For k < 5, `λ_rel = 100` (no selection).
* Prediction for the test perturbations: `P = w0 + Σ_s w_s S̃_s`.

**E4.**

* Base: E2.
* Anchor residuals: `R_K = Y_K − E2_K`.
* Perturbation feature: `z_p = concat_s(S̃_{s,p} / ‖S̃_{s,p}‖) / √3`. Kernel `K(p,q) = z_p · z_q`.
* Prediction: `R̂_T = K_{TK} (K_{KK} + λ I)^{-1} R_K`, with `λ ∈ {0.01, 0.1, 1, 10, 100}` chosen by exact
  kernel-ridge leave-one-anchor-out squared error pooled over genes.
* For k < 5, `λ = 100`.
* `P = E2 + R̂`.

No other estimator is evaluated for the pass rule. E3 and E4 are the two predeclared γ-capable estimators.

### 2.4 Evaluation quantities (evaluation-only; never inputs)

For one (target, r): `e1 = e1[r,target,T]`, `e2 = e2[r,target,T]`, `y = (e1 + e2)/2`.

**Centring for template removal.** Each matrix X ∈ {e1, e2, y, P, A} is centred by its **own** mean over
the 379 test perturbations: `X̃ = X − mean_T X`. Template-removed metrics are therefore exactly invariant to
any constant-over-perturbations template estimate. That includes E1's α̂ and E2's c, so **template estimation
cannot produce a gain in M1–M4**.

**Conserved / context-specific split, per test perturbation p.** Let `u_p = Ã_p / ‖Ã_p‖`, the direction of
the test-centred source consensus. For any centred vector x_p:

* `x∥_p = (x_p · u_p) u_p`: the component along the conserved source response ("β-part");
* `x⊥_p = x_p − x∥_p`: the target-specific component orthogonal to the source consensus ("γ⊥-part").

Any predictor of the form `f(p) · Ã_p` has `P⊥ ≡ 0`. That covers E0, E0s, any global or per-perturbation
rescaling, and any template estimate. So **only target-specific gene patterns can score on γ⊥.**

The split is linear in x given u (u is source-only), so it preserves unbiasedness.

Also `γ_ANOVA ∝ ỹ − Ã` (Identity 1 of `zero_shot_recoverability_v1.md`) = `(b_p − 1) Ã_p` + `y⊥_p`. That is,
the ANOVA interaction is an amplitude part plus γ⊥.

### 2.5 Metrics (per draw; computed separately per target context)

Unbiased energy explained, using `⟨e1 − P, e2 − P⟩` for the reliable squared error and `⟨e1, e2⟩` for the
reliable signal energy (noise in e1 and e2 is independent, because the cells and controls are disjoint).
Sums are pooled over the test perturbations. Then:

| id | name | definition |
|---|---|---|
| **M0** | full response, raw | `1 − Σ_p⟨e1−P, e2−P⟩ / Σ_p⟨e1, e2⟩` on **uncentred** matrices (the template counts) |
| **M1** | full response, template-removed | as M0 on centred matrices |
| M1r | — | median over p of Pearson `r(P̃_p, ỹ_p)` |
| **M2** | conserved β-part | `1 − Σ⟨e1∥−P∥, e2∥−P∥⟩ / Σ⟨e1∥, e2∥⟩` |
| **M3** | **γ⊥ (primary endpoint)** | `1 − Σ⟨e1⊥−P⊥, e2⊥−P⊥⟩ / Σ⟨e1⊥, e2⊥⟩` |
| M3r | — | median over p of `r(P⊥_p, y⊥_p)`, defined as 0 when `P⊥_p ≡ 0` |
| **M4** | γ_ANOVA (residual of conserved transfer) | `1 − Σ⟨ẽ1−P̃, ẽ2−P̃⟩ / Σ⟨ẽ1−Ã, ẽ2−Ã⟩`: the fraction of the reliable residual of plain source-mean transfer that P removes. It can be raised by scale alone; its γ-specific reading is **M4(E3/E4) − M4(E2)** at the same k |

Identity check (asserted): the M2 and M3 numerators and denominators sum to those of M1.

For predictors with `P⊥ ≡ 0`, M3 = 0 exactly (the numerator equals the denominator).

### 2.6 Uncertainty

* **Draw interval:** the 2.5–97.5 percentiles of the per-draw metric over the 200 draws. This captures anchor
  sampling and cell-split noise.
* **Perturbation bootstrap:** for (context, k, estimator), the per-perturbation numerator and denominator
  terms are averaged over draws. Then 2,000 bootstrap resamples of the 379 test perturbations
  (seed `[20261006, 3, context_index, k]`) give a 95 % interval of the pooled ratio. This captures which
  perturbations happen to be in T.
* Reported per context. **No pooled-across-context number is used for any decision.**

### 2.7 N1 pass rule (primary)

For γ-capable estimator E ∈ {E3, E4} and target context c, **context c passes for E** iff all three hold:

* **C1:** at k = 20, the lower bound of the draw interval for M3 is > 0, **and** the lower bound of the
  perturbation-bootstrap interval for M3 is > 0.
* **C2:** median M3 at k_ref is > 0.
* **C3:** median M3 at k = 20 is ≥ 0.25 × median M3 at k_ref.

**N1 PASSES** iff, for at least one of E3 or E4, ≥ 3 of the 4 contexts pass. Both estimators are reported in
full. The multiplicity (2 estimators) is acknowledged; no correction is applied, and the report must say so.

Outcome classes if N1 does not pass:

* **FAIL-A ("γ⊥ learnable only with a large budget"):** for E3 or E4, in ≥ 3/4 contexts, C2 holds and the
  perturbation-bootstrap lower bound of M3 at k_ref is > 0, but the small-k criteria fail.
* **FAIL-B ("γ⊥ not learnable by these estimators even at k_ref"):** otherwise.

The **smallest k** at which C1 holds is reported per (context, estimator) as the descriptive "γ⊥ onset
budget". It does not change the pass rule.

### 2.8 Secondary (descriptive, does not decide)

* Learning curves of M0–M4 and M1r/M3r for all six estimators at every k.
* Template/scale gain: M0(E1), M0(E2) − M0(E0s).
* γ_ANOVA gain beyond scale: M4(E3/E4) − M4(E2), with the draw interval.
* Selected λ distribution per k.

### 2.9 N1 failure checks (results are declared invalid if F1, F2, F4 or F5 fails)

| id | check |
|---|---|
| **F1** | **Leakage, non-anchor rows.** Replacing every non-anchor target fit row, every target e1/e2 row and the canonical target row with Gaussian noise leaves every prediction bit-identical (test) |
| **F2** | **Source-only k = 0.** E0 and E0s are bit-identical when all target data is replaced by noise (test) |
| **F3a** | **Anchor-permutation null** (reported, flagging): k ∈ {20, 100}, 50 draws per context. The anchors' fit rows are randomly permuted among the anchors before fitting E3/E4, which breaks perturbation identity. **Flag** as a possible artefact any (context, estimator) where the permuted median M3 is ≥ 50 % of the unpermuted median M3 while the latter is > 0 |
| **F3b** | **Shared-control diagnostic** (reported): k ∈ {5, 20, 100}, 50 draws per context. Anchors and truths are both re-referenced to one **common** control mean (CF ∪ CE1 ∪ CE2). Report M0 and M3 for E1, E3, E4 against the split-control values. The difference quantifies how much shared control noise would have inflated the gains |
| **F4** | **Disjointness:** F/E1/E2 and CF/CE1/CE2 index sets are pairwise disjoint (asserted at build time, recorded) |
| **F5** | **Synthetic recovery** (test): on synthetic data with a planted low-rank, perturbation-structured γ⊥, E4 at k = 100 attains M3 > 0.2. With γ⊥ = 0 and pure noise, E3/E4 M3 ≤ 0.02 in expectation |
| **F6** | **Decomposition identity:** M2 + M3 parts sum to M1 (asserted to 1e-8 relative) |

---

## 3. N4: does source agreement predict transfer beyond simple statistics?

### 3.1 Units and score

* Per held-out context c (4 LOCO folds), all 1,264 perturbations.
* **Score:** the frozen source agreement `a_p` (mean pairwise Pearson of the three canonical source responses;
  `transferability.source_agreement`), unchanged.

### 3.2 Held-out transfer quality (target; evaluation-only)

* Predictor: the frozen B (= E0s) for that fold.
* Target halves, averaged over the 5 repeats:
  * `h1 = fit[·, c]` (F cells vs CF controls);
  * `h2 = (e1 + e2)/2 [·, c]` (E1∪E2 cells vs CE1/CE2 controls).
  * These are disjoint in both perturbed and control cells.
* **Primary quality:** `q_p = −D_p`, with `D_p = ⟨h1 − B, h2 − B⟩ / ⟨h1, h2⟩` (raw, uncentred, as in
  `transferability_confidence_model_v1.md`). Masked to NaN where `⟨h1, h2⟩ ≤ 1.96738` (the frozen stability
  rule).
* **Secondary quality:** `r_p = Pearson(B_p, (h1 + h2)/2)` over genes, all perturbations.

### 3.3 Simple competing statistics (all source-only or target-basal; inference-available)

| id | statistic |
|---|---|
| b1 | source magnitude `‖A_p‖` |
| b2 | source reliability: mean over the 3 sources of Pearson(`fit[r,s,p]`, `(e1+e2)/2[r,s,p]`), averaged over r |
| b3 | source reliable energy: mean over sources of `⟨fit_s, (e1+e2)/2_s⟩`, averaged over r |
| b4 | minimum source cell count |
| b5 | target-gene basal expression in the target's canonical control mean (NaN for the 192 off-axis targets; those rows are excluded only from analyses using b5) |
| b6 | noise-implied agreement ceiling: mean over source pairs of `√(ρ_a ρ_b)`, with `ρ` the Spearman–Brown full-depth reliability from b2's half correlation (clipped to [0, 1]) |

### 3.4 Tests (all per context; never pooled for a decision)

* **T1, marginal comparison.**
  * Spearman(statistic, q) for the agreement and for each of b1–b6.
  * Paired bootstrap (2,000 resamples of perturbations; seed `[20261006, 4, c]`) 95 % interval of
    ρ(agreement) − ρ(b_j).
  * For b1, b4 and b5 the sign of ρ is taken as estimated (the competitor gets its best orientation): the
    comparison uses |ρ(b_j)|.
* **T2, incremental information.** Partial Spearman of the agreement with q controlling for b1, b2, b3, b4 and
  b5 jointly (`foundations.partial_spearman`), with a 2,000-resample bootstrap 95 % interval. Rows with NaN b5
  are excluded from T2.
* **T3, exchangeable random-effects null** (descriptive classification):
  * All responses are centred over the 1,264 perturbations within each context (template removed).
  * For each perturbation p, using the sources:
    * `m_p = Ã_p` (centred source mean);
    * per-gene noise variance `σ²_{s,p} = ‖fit − (e1+e2)/2‖² / (4G)` for source s, averaged over r;
    * between-context variance `τ²_p = max(0, Σ_s ‖S̃_{s,p} − m_p‖² / (2G) − mean_s σ²_{s,p})`;
    * target noise `σ²_{c,p}` from the target's own halves (a measurement property, used only by the null).
  * 200 simulations (seed `[20261006, 5, c, i]`). Draw independently:
    * `u ~ N(0, τ²_p I)` for each of the 4 contexts;
    * source responses `S'_s = m_p + u_s + N(0, σ²_{s,p} I)`;
    * target halves `h'_i = m_p + u_c + N(0, 2σ²_{c,p} I)`.
  * Recompute, in the simulated world: the agreement `a'`, `B' = s · mean_s S'_s` (same frozen s), and
    `q' = −D'` with the same stability mask rule.
  * Null distribution = Spearman(a', q') over the 200 simulations.
  * The **observed** statistic for T3 is Spearman(a, q_centred), with a and q recomputed on the same centred
    responses, for like-for-like comparison.
  * **Classification:**
    * "below null": observed < 2.5th percentile. Contexts are less exchangeable than the random-effects model;
      agreement over-promises.
    * "within null": agreement behaves exactly like a random-effects heterogeneity statistic.
    * "above null": observed > 97.5th percentile. Agreement carries information beyond exchangeable
      heterogeneity.
  * Plug-in `m_p` double-counts some noise and heterogeneity. This known limitation is stated in the results.

### 3.5 N4 pass rule

**N4 PASSES** iff, in ≥ 3 of the 4 contexts, both hold:

* **(i)** the T2 partial-Spearman 95 % interval lower bound is > 0; and
* **(ii)** the T1 difference intervals ρ(agreement) − |ρ(b_j)| have lower bounds > 0 for **each** of b1, b2,
  b3 and b6.

These are the signal-strength and reliability explanations. b4 and b5 are reported but are not required.

Reported tiers:

* **"N4 weak"**: (i) holds in ≥ 3/4 contexts but (ii) does not. Agreement adds information but does not
  dominate the simple statistics.
* **T3 class per context:** reported. It decides the framing (heterogeneity statistic vs beyond), not pass or
  fail.

---

## 4. What each outcome will be taken to mean (fixed now)

**N1 PASS.** A few target perturbations improve held-out γ⊥, not just the template and scale. P1 (calibration
budget) proceeds, and N2 (selection) becomes the next experiment, but only after a new predeclaration.

**N1 FAIL-A.** γ⊥ is learnable, but only at budgets near k_ref. That is the regime already reported by State
and Molina & Zhang. P1 shrinks to a limits result: "the template and scale are cheap, γ is expensive". N2 is
not justified.

**N1 FAIL-B.** These low-capacity estimators cannot learn γ⊥ even from 885 anchors at half depth. That says
nothing about higher-capacity models. It argues against P1 as framed.

**N4 PASS.** Agreement is a genuine transferability signal beyond signal strength and reliability.

**N4 weak, or N4 FAIL.** Agreement is mainly a signal-strength / reliability proxy. It must be presented as a
baseline, not as a contribution.

**T3 "within null".** Agreement behaves as random-effects heterogeneity theory predicts. It is not novel, but
it is well understood.

## 5. Outputs

* `outputs/n1_n4/` (git-ignored): splits, per-draw metric tables (`n1_draws.parquet`), bootstrap tables,
  F3a/F3b tables, N4 tables (`n4_t1.csv`, `n4_t2.csv`, `n4_t3_null.csv`), figures, `decision.json`, and a
  SHA-256 manifest.
* Small figure-source CSVs under `data/figure_sources/n1_n4/`.
* `reports/n1_n4_results.md`.
