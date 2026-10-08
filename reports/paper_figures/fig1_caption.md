# Figure 1 | Conserved and context-specific perturbation effects differ in reproducibility

**Question.** What biological information is a cross-context predictor trying to transfer? *Background and
replication; not the paper's novel contribution.*

**Data.** Four public CRISPRi Perturb-seq contexts (K562, RPE1: Replogle et al. 2022, CC BY 4.0; HepG2, Jurkat:
Nadig et al. 2025), scPertEval preprocessing; 1,264 perturbations shared by all four × 6,640 shared genes.

**(A) Prediction setting and information boundary.** Source contexts contribute measured perturbation responses
(pseudobulk mean log-expression change versus non-targeting controls). In the held-out target context only control
cells are available. *Zero shot*: source responses + target controls. *Target calibration*: the same information plus
the measured responses of k target perturbations ("anchors"); the evaluated perturbations are never anchors. Key:
grey = measured and available to the predictor; black = measured target anchors; hatched = unknown, held out for
evaluation.

**(B) Response decomposition.** δ(c, p) = μ + α(c) + β(p) + γ(c, p): global response μ and context-wide component
α(c) (together the *template*), conserved perturbation effect β(p) (the same in every context) and context-specific
perturbation interaction γ(c, p). The strip is illustrative, not data. The framework follows Molina & Zhang (2026);
it is independently re-derived here as background, not presented as a new method, and is not a reproduction of
their data.

**(C) Response energy.** Share of total response energy (sum of squares) per component after split-half noise
correction (filled markers; mean of 50 split-half resamples; s.d. over resamples 0.010–0.044 percentage points,
smaller than the markers). Open markers and arrows: shares before correction (β 37.2 %, γ 42.5 %). Shaded bars:
range over the 21 predeclared preprocessing variants (Extended Data Fig. 1) — a robustness range, not a confidence
interval.

**(D) Reproducibility.** Fraction of each component's raw sum of squares reproducible across disjoint split halves
of cells (cross-half signal / raw SS), on the same rows as C. Bars: range over the 21 variants (β 77.4–85.0 %,
γ 45.6–57.7 %); template ≥ 98.7 % in every variant (range not drawn).

**Reading (conservative).** The context-specific interaction is not negligible — 21.0 % of response energy, as large
as the template — but only 49.5 % of it is reproducible, against 80.8 % for β; γ carries about three-quarters of
all measurement noise. This motivates Fig. 2: how much target data is needed to estimate the reproducible part of γ.

**Analysis status.** Predeclared decomposition and sensitivity battery (descriptive; no hypothesis test).
**Limitations.** Four essential-gene screens from two labs; one public preprocessing release.

*Source data:* `data/figure_sources/paper/fig1_components.csv`. Script: `scripts/paper_figures/plot_fig1.py`.
