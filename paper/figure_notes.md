# Figure notes

Populated only from frozen outputs, frozen reports and the figure captions. Numbers are quoted, not recomputed.
Figure files: `reports/paper_figures/`. Rebuild: `bash scripts/paper_figures/run_all.sh`.

---

## Figure 1 — Conserved and context-specific perturbation effects differ in reproducibility

* **Scientific question:** what biological information does a cross-context predictor try to transfer?
* **Supported claim:** the context-specific interaction γ is a substantial part of the response (21.0 % of energy
  after noise correction) but is much less reproducible than the conserved effect β (49.5 % vs 80.8 %).
* **Evidentiary status:** background / replication of the Molina & Zhang (2026) decomposition (independent
  re-derivation; predeclared decomposition and 21-variant sensitivity battery). Not a novelty claim.
* **Primary numerical results:**
  * energy shares (noise-corrected): template 20.3 %, β 30.1 %, γ 21.0 %, noise 28.6 %;
  * uncorrected β 37.2 %, γ 42.5 %;
  * reproducibility: β 80.8 % (21 variants 77.4–85.0 %), γ 49.5 % (45.6–57.7 %);
  * γ carries about three-quarters of all measurement noise.
* **Caveats:** four essential-gene CRISPRi screens from two labs; one public preprocessing release; the variant
  ranges are robustness ranges, not confidence intervals.
* **Frozen source:** `outputs/four_context_v1/summary.json`, `outputs/four_context_sensitivity/variant_table.csv`;
  reports `four_context_decomposition_v1.md`, `four_context_decomposition_sensitivity.md`.
* **Caption:** `reports/paper_figures/fig1_caption.md`

## Figure 2 — Context-specific perturbation effects require substantially more target data to estimate

* **Scientific question:** how many measured target-context perturbations does the frozen estimator need to recover
  each response component?
* **Supported claim:** under the frozen estimator and random-anchor design, template/scale adaptation is largely
  complete within tens to about 100 anchors. Context-specific γ⊥ is recovered slowly and is incomplete at
  743 anchors.
* **Evidentiary status:**
  * N1: confirmatory, four contexts, FAIL-A.
  * N3-A: preregistered replication on six contexts with the protocol unchanged; R-A1 holds.
  * Panel D draw intervals and all cross-context medians are descriptive.
* **Primary numerical results:**
  * γ⊥ recovery is 0.02–0.15 at k = 20 and 0.26–0.41 at k_ref = 743.
  * k20/k_ref is 0.09–0.36. The preregistered C3 bar (0.25) is met only in RPE1, so FAIL-A holds in 5/6 contexts
    (and in 3/4 in N1).
  * Template + scale reaches half of its k_ref gain at k_T50 = 10–100, and ≥ 86 % by k = 200 in every context.
  * G(k) for γ⊥ is 0.09–0.36 at k = 20, 0.23–0.55 at k = 50 and 0.40–0.68 at k = 100.
  * G(k) for template + scale is 0.35–0.93 at k = 50 and 0.70–0.96 at k = 100.
  * A single anchor is harmful (normalised template gain down to −29).
* **Caveats:**
  * One deliberately simple linear estimator on half-depth anchors; random anchors only.
  * Not a bound for other estimators.
  * Lab and assay are partly confounded for the X-Atlas contexts.
  * Recovery of γ at ~30 % of the panel was already reported by Molina & Zhang and by State.
* **Frozen source:** `outputs/n3/n3a/{n3a_summary,n3a_boot,n3a_draws}.csv`, `outputs/n3/n3_decision.json`;
  reports `n3_results.md` §1, `n1_n4_results.md` §2.
* **Caption:** `reports/paper_figures/fig2_caption.md`

## Figure 3 — Broader source panels improve zero-shot transfer across most target contexts

* **Scientific question:** does adding source contexts improve zero-shot transfer to a held-out context?
* **Supported claim:** increasing source breadth improves zero-shot transfer in most contexts (5/6), with HEK293T
  as an important exception. No claim about whether breadth can substitute for target-context measurements.
* **Evidentiary status:** preregistered N3-B quantity (zero-shot template-removed full response, frame-free), read
  descriptively; no preregistered trend test. Panel C re-displays panel B's values.
* **Primary numerical results:**
  * zero-shot M1 change from m = 2 to 5: K562 +0.15, HepG2 +0.14, Jurkat +0.13, HCT116 +0.11, RPE1 +0.10;
  * HEK293T −0.16 (+0.05 → −0.11);
  * m = 1 values are −0.19 to −5.17 (single unshrunk source).
* **Caveats:**
  * Six targets from three labs.
  * The repeat ranges reflect cell-split variation only.
  * HEK293T is atypical (lowest knockdown, weakest correlation with other contexts) for reasons that cannot be
    separated.
  * The preregistered Outcome-B criterion and the exploratory own-frame analyses are in ED 5, not here.
* **Frozen source:** `outputs/n3/n3b/n3b_R_gamma.csv` (k = 0 rows), `outputs/n3/n3b/n3b_draws.csv`
  (repeat ranges); report `n3_results.md` §3.
* **Caption:** `reports/paper_figures/fig3_caption.md`

## Figure 4 — Source choice strongly influences perturbation transfer

* **Scientific question:** if source quantity is not the whole story, how much does the choice of source matter?
* **Supported claim:** within one study, the choice of source strongly influences transfer to the K562-essential
  target. The difference is not explained by source reliability, depth or perturbation signal strength alone. Cell
  identity cannot be separated from study/lab effects.
* **Evidentiary status:** confirmatory (N5 Outcome A: primary contrast and preregistered controls). Panel D
  summarises design facts, including N6 Outcome D (uninformative).
* **Primary numerical results:**
  * C_S: K562 GWPS 0.822 [0.808, 0.837]; depth-matched GWPS 0.818; Jurkat 0.558; HepG2 0.490; HCT116 0.456;
    HEK293T 0.387; RPE1 0.379.
  * ΔC_S (GWPS − RPE1) = +0.444 [0.420, 0.469]; Adj-4 +0.439.
  * Adj-2 +0.221 at full depth, +0.181 depth-matched; Adj-3 +0.233 (n = 108); raw median r difference +0.138.
  * 85.8 % (full) and 70.5 % (depth-matched) of perturbations favour GWPS.
  * The preregistered ladder fails: RPE1 is last.
* **Caveats:**
  * One target and one same-cell source, from the same study and lab.
  * No technical replicate; GWPS differs in library and timepoint.
  * Largely expected from Replogle 2022 and Nadig 2025.
  * N6 (VIPerturb-seq) failed its reliability gate (0.0866 < 0.10).
* **Frozen source:** `outputs/n5/{n5_decision,n5_summary}.json`, `outputs/n5/per_perturbation_terms.csv`;
  reports `n5_results.md`, `n6_results.md`.
* **Caption:** `reports/paper_figures/fig4_caption.md`

---

## Extended Data Figure 1 — Decomposition robustness

* **Scientific question:** does Fig. 1 depend on a single preprocessing choice?
* **Supported claim:** across 21 predeclared variants, β never falls below 28.0 % and corrected γ never below 20.5 %.
  Uncorrected γ is about twice corrected γ in every variant. Reliability rises with cells per half.
* **Evidentiary status:** predeclared sensitivity battery (background).
* **Primary numerical results:**
  * β 28.0–30.8 %; γ 20.5–22.8 %;
  * reproducibility: β 77.4–85.0 %, γ 45.6–57.7 %;
  * uncorrected γ 38.5–45.1 %;
  * median split-half reliability 0.10 / 0.18 / 0.26 / 0.42 at 15 / 30 / 50 / 100 cells per half.
* **Caveats:** preprocessing variants on one data release only. The depth CIs are transcribed from the frozen
  report; their medians were asserted against `depth_experiment.csv`.
* **Frozen source:** `outputs/four_context_sensitivity/{variant_table,depth_experiment}.csv`; report
  `four_context_decomposition_sensitivity.md`.
* **Caption:** `reports/paper_figures/ext_fig1_caption.md`

## Extended Data Figure 2 — Validity controls for the target-budget experiment

* **Scientific question:** could Fig. 2 be produced by leakage, shared control noise or chance pairing?
* **Supported claim:** no. Permuted-anchor E4 recovery is about zero; the γ⊥ endpoint is identical with shared vs
  split controls; a single anchor is harmful; N1 and N3-A agree.
* **Evidentiary status:** preregistered failure checks (F3a, F3b); the panel D intervals are descriptive.
* **Primary numerical results:**
  * permuted E4 M3 medians lie between −0.0014 and +0.0007;
  * shared-control M0 inflation is up to +0.027, while M3 is identical;
  * E1 ΔM0 at k = 1 is −0.47 to −2.51;
  * k20/k_ref, N1: 0.13, 0.39, 0.16, 0.13 (K562, RPE1, HepG2, Jurkat);
  * k20/k_ref, N3-A: 0.10, 0.36, 0.12, 0.10, 0.09, 0.13 (adding HCT116, HEK293T).
* **Caveats:** F1, F2, F4, F5 and F6 are tests in the repository and are not plotted.
* **Frozen source:** `outputs/n3/n3a/{n3a_f3a,n3a_f3b_summary,n3a_draws}.csv`, `outputs/n1_n4/n1/n1_draws.csv`.
* **Caption:** `reports/paper_figures/ext_fig2_caption.md`

## Extended Data Figure 3 — Correction of the source-agreement estimate

* **Scientific question:** how strongly does source agreement predict transfer quality, once correctly estimated?
* **Supported claim:** the v1 estimate (0.55–0.79) was inflated by averaging halves over repeats; the corrected
  estimate is 0.29–0.59. Agreement carries only modest incremental information (N4 WEAK) and is a baseline, not a
  contribution.
* **Evidentiary status:** panel A exploratory; panels B–C confirmatory (N4, N3 six folds).
* **Primary numerical results:**
  * per-repeat ρ: K562 0.472, RPE1 0.589, HepG2 0.285, Jurkat 0.438;
  * partial ρ (N1/N4): 0.12–0.22;
  * partial ρ (N3): positive in 5/6 folds, HCT116 −0.07 [−0.16, 0.03];
  * agreement beats all four required baselines in 4/6 N3 folds and is worse than source reliability in HCT116.
* **Caveats:** the frozen v1 report is unchanged and carries a dated erratum.
* **Frozen source:** `data/figure_sources/n1_n4/{exploratory_halves_bias,n4_t1,n4_t2}.csv`,
  `outputs/n3/n4/{n4_t1,n4_t2}.csv`; report `n1_n4_results.md` §3.
* **Caption:** `reports/paper_figures/ext_fig3_caption.md`

## Extended Data Figure 4 — N6 reliability gate

* **Scientific question:** why is the independent-lab K562 test (N6) uninformative?
* **Supported claim:** VIPerturb-seq failed the preregistered pooled-reliability gate (0.0866 < 0.10), so N6 is
  Outcome D. The GWPS positive control at matched depth was stable. VIPerturb compatibility values are
  non-interpretable and are omitted.
* **Evidentiary status:** gate items are confirmatory. The pooled GWPS values (0.109 at VIPerturb depth, 0.341 at
  full depth) are exploratory post-hoc diagnostics.
* **Primary numerical results:**
  * d1 0.0866 (fail); d2 2,481 (pass); d3 |ΔC| 0.007 (pass); d4 637 perturbations (pass);
  * C_S GWPS 0.818 [0.796, 0.838] vs 0.811 [0.785, 0.834] at VIPerturb depth;
  * fraction of perturbations with reliability ≥ 0.10: VIPerturb 22 %, GWPS at VIPerturb depth 27 %, GWPS full
    82 %.
* **Caveats:** the same-cell / different-study question remains unanswered, not answered negatively. No other
  eligible public K562 dataset was identified (literature search; N6 results §5).
* **Frozen source:** `outputs/n6/{source_gate,n6_decision}.json`, `outputs/n6/per_perturbation_terms.csv`; report
  `n6_results.md`.
* **Caption:** `reports/paper_figures/ext_fig4_caption.md`

## Extended Data Figure 5 — Source count with target anchors: preregistered fixed frame and exploratory own frame

* **Scientific question:** with target anchors, does source breadth change how much γ⊥ the anchors recover? How do
  the preregistered and exploratory endpoints differ?
* **Supported claim:**
  * The preregistered fixed-frame Q1 criterion fired: Δ(20) is material in 5/6 targets.
  * The same pairing in each subset's own frame is about zero.
  * The fixed-frame increase with m tracks a zero-shot offset created by the fixed reference.
  * The disagreement is reported, not resolved. Neither frame alone establishes whether source breadth can
    substitute for target measurements.
* **Evidentiary status:**
  * Panel A fixed frame: confirmatory (criterion fired; Outcome B mechanically met).
  * Panel A own frame, panel B, the own-frame curve in C, and panel E: exploratory.
  * Panel D: preregistered quantities, read descriptively.
* **Primary numerical results:**
  * fixed-frame Δ(20) 0.027–0.073; own-frame Δ(20) −0.000 to +0.009, with every interval including 0;
  * own-frame gain change from m = 1 to 5: −0.005 to +0.045 at k = 20, +0.003 to +0.045 at k = 50;
  * fixed-frame zero-shot offsets are negative for m < 5;
  * the 20-anchor M1 gain is constant or smaller with more sources in 5/6 targets (K562 0.070 → 0.027).
* **Caveats:** the own-frame endpoint was defined after the preregistered analysis, and its target depends on
  each subset's consensus.
* **Frozen source:** `outputs/n3/n3b/{n3b_R_gamma,n3b_draws}.csv`, `outputs/n3/n3_decision.json`; report
  `n3_results.md` §3, §6–7.
* **Caption:** `reports/paper_figures/ext_fig5_caption.md`
