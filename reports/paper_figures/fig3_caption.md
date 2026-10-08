# Figure 3 | Broader source panels improve zero-shot transfer across most target contexts

**Question.** Does adding source contexts improve zero-shot transfer to a held-out context?

**Data and design.** Preregistered N3-B source-count ablation on the six contexts of Fig. 2 (K562, RPE1, HepG2, Jurkat
from scPertEval; HCT116, HEK293T from X-Atlas/Orion, Xaira Therapeutics, CC BY-NC-SA 4.0, non-commercial research
use, aggregate statistics only). 1,062 perturbations × 6,499 genes; for each held-out target, the 319 N3-A test
perturbations. Every non-empty subset of the target's five candidate sources is evaluated (31 subsets, m = 1–5);
values are averaged over all subsets of size m, and no subset is selected. **This figure uses only the zero-shot
(k = 0) predictor and a frame-free endpoint.** All analyses with target anchors (k > 0) and all own-frame
(exploratory) analyses are in Extended Data Fig. 5.

**(A) Design.** Source breadth m and target measurements k were varied independently in N3-B; this figure shows k = 0.

**(B) Zero-shot full response versus source count.** Template-removed full-response energy explained (M1; unbiased,
pooled over the 319 test perturbations) by the zero-shot predictor (source consensus with a source-fitted scale),
one line per held-out target; markers at each tested m with bars showing the range over the 5 cell-split repeats
(mostly narrower than the markers); thick line: median of the six targets (descriptive; no interval). m = 1 values
are off-scale (−0.19 to −5.17) because a single source is used without the source-fitted shrinkage (fixed at 1 for
m = 1); they are excluded from the median line.

**(C) Change from two to five sources, per target.** The m = 2 (grey) and m = 5 (black) values of panel B, with the
same repeat ranges, and their difference printed at right. No new quantity is computed. Zero-shot M1 rises by +0.10
to +0.15 in five targets (K562 +0.15, HepG2 +0.14, Jurkat +0.13, HCT116 +0.11, RPE1 +0.10) and falls in HEK293T
(−0.16; from +0.05 to −0.11).

**Reading.** Increasing source breadth improves zero-shot transfer in five of six held-out contexts; HEK293T is an
important exception, in which more sources make zero-shot prediction worse. This figure makes no claim about whether
source breadth can or cannot substitute for target-context measurements: the analyses bearing on that question are
exploratory and are shown, together with the preregistered fixed-frame test and its disagreement with them, in
Extended Data Fig. 5.

**Analysis status.** Zero-shot M1 per source count is a preregistered N3-B quantity (`R_full(0, m)`), read here
descriptively; no test was preregistered for its trend.
**Limitations.** Six targets from three labs; repeat ranges capture cell-split variation only, not test-perturbation
sampling; HEK293T is atypical in several respects (lowest knockdown, weakest correlation with every other context)
that cannot be separated from its source-count behaviour.

*Source data:* `data/figure_sources/paper/fig3_source_count.csv` (k = 0 rows; from frozen
`outputs/n3/n3b/n3b_R_gamma.csv` and `n3b_draws.csv`). Script: `scripts/paper_figures/plot_fig3.py`.
