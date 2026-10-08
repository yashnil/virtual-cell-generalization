# Extended Data Figure 3 | Correction of the source-agreement estimate; agreement is a baseline, not a contribution

**(A) Exploratory re-estimation (not preregistered).** Spearman correlation between source agreement (mean pairwise
Pearson of the source responses) and held-out transfer quality (−D, reliability-normalised error of the frozen
zero-shot predictor) on the same disjoint N1 splits. Left: target halves averaged over the five repeats before
forming energies (the construction of the frozen v1 report; grey ticks: the v1 report's values, 0.55–0.79).
Right: per-repeat unbiased energies (the N4 estimator, deviation D1 declared before any N4 output). Averaging
halves over repeats re-uses cells across halves, inflates the stability filter (≈ 99.6 % vs 65–84 % of
perturbations passing) and the correlation. Corrected range 0.29–0.59; the frozen v1 report is unchanged and
carries a dated erratum.

**(B) Incremental information (preregistered T2).** Partial Spearman of agreement with transfer quality controlling
for source magnitude, source reliability, source reliable energy, minimum source cell count and target-gene basal
expression; 95 % bootstrap intervals (2,000). N1/N4: four folds, three sources; N3: six folds, five sources.
Positive in 5/6 N3 folds (0.12–0.27); not distinguishable from 0 in HCT116.

**(C) Marginal comparison (preregistered T1, N3 six folds).** ρ(agreement) − |ρ(baseline)| with paired bootstrap
95 % intervals for the four baselines required by the pass rule. Filled black circle: agreement
better; open black triangle: agreement worse; grey: not distinguishable. Agreement beats all four only in four of six folds and is worse than plain
source reliability in HCT116.

**Reading.** Preregistered verdict N4: WEAK in both runs. Agreement is a mostly signal-strength-driven trust
baseline; it should not be presented as stronger than source reliability.

**Purpose / why not main text.** A transparent correction of a previously reported estimate. Source agreement is a
baseline (N4 WEAK), not a contribution, so it does not belong in the main argument.

*Source data:* `data/figure_sources/paper/ext3_*.csv`. Script: `scripts/paper_figures/plot_ext_fig3.py`.
