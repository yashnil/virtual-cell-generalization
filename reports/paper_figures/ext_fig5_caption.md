# Extended Data Figure 5 | Source count with target anchors: preregistered fixed-frame test and exploratory own-frame analyses

**Purpose.** Report, side by side, the preregistered N3-B source-count test and the exploratory analyses that led us
to judge it confounded, so that the disagreement between them is visible. Data and design as in Fig. 3 (N3-B; six
held-out targets; 31 source subsets; 40 random anchor draws per budget, identical across subsets; 319 test
perturbations per target; X-Atlas/Orion CC BY-NC-SA 4.0, research use).

**(A) Preregistered result (primary).** Δ(20) = R(k = 20, m = 5) − R(k = 20, m = 3), the preregistered Q1 quantity,
paired over the 40 draws (median and 2.5–97.5 % draw interval). Black: the preregistered fixed frame (γ⊥ relative to
the full 5-source consensus; values copied from `n3_decision.json`): 0.027–0.073, material by the preregistered rule
(Δ ≥ 0.03 with lower bound > 0) in 5 of 6 targets; together with Q2 this **fired the Outcome-B criterion**.
Vermilion: the same pairing with γ⊥ measured in each subset's own frame (exploratory): −0.000 to +0.009, interval
including 0 in all six targets. The two frames disagree; both are reported.

**(B) Own-frame calibration gain versus m (exploratory; formerly Fig. 3C).** γ⊥ recovered from k = 20 and k = 50
anchors relative to each subset's own consensus, so that zero-shot is 0 and the value is the calibration gain.
Open diamond m = 1, ticks m = 2–4, filled diamond m = 5; bars: 95 % draw intervals; right column: change from m = 1
to 5 (−0.005 to +0.045 at k = 20; +0.003 to +0.045 at k = 50; the larger increases are in K562 and RPE1).

**(C) Why the frames differ.** Per target, fixed-frame γ⊥ at k = 0 (grey open) and k = 20 (black; band 95 % draw
interval) and own-frame γ⊥ at k = 20 (vermilion). In the fixed frame a predictor built from fewer sources has
conserved-estimate error orthogonal to the 5-source consensus, which scores as negative γ⊥ before any anchor is used
(0 at m = 5 by construction; m = 1 values printed, off-scale). The fixed-frame k = 20 curve rises with m largely
because its k = 0 baseline does.

**(D) Frame-free full response with anchors.** Template-removed full response (M1), zero-shot (open; bars: range
over 5 repeats; the Fig. 3B values) and after 20 anchors (filled; 95 % draw interval). From 2 to 5 sources the
20-anchor gain is constant or smaller in five targets (e.g. K562 0.070 → 0.027) and larger only in HEK293T, whose
zero-shot start deteriorates.

**(E) Own-frame k × m surface (exploratory; formerly Fig. 3D).** Own-frame γ⊥ for every tested budget (k = 5–743)
and source count, one heatmap per target, identical colour limits (0–0.42); grey: below 0 (minimum −0.053, k = 5);
k = 0 omitted (0 by construction). No uncertainty is drawn on cells; intervals for k = 20 and 50 are in B.

**Interpretation constraints.** The own-frame endpoint was defined after the preregistered analysis (N3 results
§3.3–3.4) and its target depends on each subset's consensus, so it changes slightly with m. Neither frame alone
establishes that source breadth can or cannot substitute for target-context measurements; the main text (Fig. 3)
therefore relies only on the frame-free zero-shot result. The mechanical preregistered verdict and the judgement to
depart from it are recorded in `reports/n3_results.md` §6–7.

**Why not main text.** Exploratory analyses and an audit trail for a confounded preregistered criterion.

*Source data:* `ext5_delta20.csv`, `ext5_reference_frame.csv` (frozen `n3b_R_gamma.csv`), `fig3_source_count.csv`.
Script: `scripts/paper_figures/plot_ext_fig5.py`.
