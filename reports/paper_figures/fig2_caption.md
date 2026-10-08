# Figure 2 | Context-specific perturbation effects require substantially more target data to estimate

**Panel A reports recovery relative to each context's own high-budget gain, whereas Panel B reports absolute
reliable γ⊥ recovery; Panel D provides the directly comparable normalized summary.** Recovery is measured under
the frozen estimator and random-anchor design; the curves characterize this evaluation framework, not a universal
bound for other estimators or anchor-selection strategies.

**Question.** How many measured target-context perturbations does this estimator need to recover each component of
the response?

**Data and design.** Preregistered N3-A experiment: the confirmatory N1 protocol (four contexts) re-run unchanged on
six contexts — K562, RPE1, HepG2, Jurkat (scPertEval) and HCT116, HEK293T (X-Atlas/Orion, Xaira Therapeutics; Huang
et al. 2025; CC BY-NC-SA 4.0, used for non-commercial research only; aggregate statistics only). 1,062
perturbations × 6,499 genes. Each context is held out in turn with the other five as sources; per target a fixed
test set of 319 perturbations and an anchor pool of 743. For k ∈ {1, 2, 5, 10, 20, 50, 100, 200}, anchors are drawn
uniformly at random from the pool, 200 draws per k (40 anchor sets × 5 independent cell-split repeats); k_ref = 743
(the whole pool). Anchors are measured on half of each perturbation's cells and half of the controls; the truth
uses the two remaining disjoint quarters, so all energies are estimated without noise bias. k = 0 is the zero-shot
predictor, set off by an axis break and not connected to k ≥ 1. Markers show the tested budgets; lines only guide
the eye; vertical guides mark k = 20, 50, 100 and k_ref.

**(A) Template + scale (normalised).** Estimator E2 (target template plus least-squares scale from the anchors):
gain in full-response energy explained (M0, template included) over the zero-shot predictor, divided by the same gain
at k_ref, per context (median over draws). Thin lines: six held-out contexts; thick line: their median (descriptive;
no bootstrap exists for this ratio or for the cross-context median). Values below −1 (k ≤ 10; minimum −29, a single
anchor's own response mistaken for the template) are drawn at −1 (▼). Half of the k_ref gain is reached at
k_T50 = 10–100 depending on the context, ≥ 86 % by k = 200 in every context.

**(B) Context-specific γ⊥ (absolute).** Preregistered γ-capable estimator E4 (E2 plus kernel-ridge transfer of
anchor residuals across perturbations similar in the sources). Endpoint M3: unbiased fraction of the reliable
held-out γ⊥ energy explained, where γ⊥ is the target response orthogonal, perturbation by perturbation, to the
source consensus. Template, global-scale and per-perturbation-scale predictors score exactly 0 on M3 by construction,
so any gain is context × perturbation information. Not normalised (1 = all reliable γ⊥ energy). Median over draws;
lines as in A. 0.02–0.15 at k = 20; 0.26–0.41 at k_ref and still rising.

**(C) Six held-out contexts.** M3 point estimate and 95 % percentile interval from a perturbation bootstrap (2,000
resamples of the 319 test perturbations on draw-averaged terms; the point can differ slightly from the draw median
in B). Identical axes, ticks and budgets. Corner text: median M3(k = 20)/M3(k_ref) and the class under the
preregistered N1 rule (PASS needs ≥ 0.25 at k = 20 with both lower bounds > 0; FAIL-A = recovery positive at k_ref but not at small budgets).

**(D) Like-for-like comparison.** G(k) = [R(k) − R(0)] / [R(k_ref) − R(0)] at k = 20, 50, 100, for γ⊥ (E4 M3;
R(0) = 0; filled diamonds) and template + scale (E2 M0 gain; open squares). Point: ratio of draw medians (for γ⊥
at k = 20, the preregistered C3 quantity). Interval: 2.5–97.5 % over the 200 draws of the per-draw ratio with
k_ref taken from the same cell-split repeat (descriptive; it reflects anchor and cell sampling, not test-perturbation
sampling); intervals below −0.75 truncated (◀). Bottom row: median of the six contexts (descriptive, no interval).
Dashed line at k = 20: preregistered C3 bar (0.25). γ⊥ reaches 0.09–0.36 of its k_ref value at k = 20, 0.23–0.55 at
k = 50 and 0.40–0.68 at k = 100; template + scale 0.35–0.93 at k = 50 and 0.70–0.96 at k = 100.

**Reading.** Template/global adaptation is largely complete within tens to ~100 random anchors, whereas context-
specific perturbation information accumulates roughly log-linearly and is incomplete at 743 anchors (70 % of the
panel). The preregistered small-budget criterion failed (FAIL-A) in 5/6 contexts here and 3/4 in the confirmatory
N1 run (Extended Data Fig. 2D); RPE1, the most reliable screen, is the exception.

**Analysis status.** Confirmatory (N1) and preregistered replication (N3-A); D's intervals are descriptive.
**Limitations.** Random anchors only (selection untested); a deliberately simple linear estimator; essential-gene
screens from three labs, with lab and assay partly confounded for the X-Atlas contexts. Recovery of γ from ~30 %
of target perturbations was reported by Molina & Zhang (2026) and State; the contribution here is the
component-resolved budget curve across six contexts.

*Source data:* `fig2_curves.csv`, `fig2_gain_fraction.csv`. Script: `scripts/paper_figures/plot_fig2.py`.
