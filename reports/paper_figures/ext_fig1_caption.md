# Extended Data Figure 1 | The β/γ decomposition is robust to preprocessing

Predeclared sensitivity battery for the four-context decomposition (Fig. 1): 21 variant × feature-space
combinations (shared control mean vs independent control split; mean of log1p vs log of mean CP10K; five seeds;
all 6,640 genes, 4,000 or 2,000 control-defined HVGs). Filled marker: canonical analysis; open markers: the other
variants (vertical jitter only for visibility).

**(A)** Noise-corrected energy share per component; printed numbers are the min–max over variants. β never falls
below 28.0 % and γ never below 20.5 %. Template and noise shares move with the gene set (HVG subsets carry more
template and less noise), as expected.
**(B)** Split-half reproducibility of β and γ: β 77.4–85.0 %, γ 45.6–57.7 %; the ordering holds in every variant.
**(C)** γ share before and after noise correction, one line per variant (canonical bold): uncorrected γ
(38.5–45.1 %) is roughly twice corrected γ in every variant.
**(D)** Controlled depth experiment: the same 643 (context, perturbation) pairs re-estimated at 15, 30, 50 and
100 cells per split half (25 repeats each); median split-half reliability with 95 % bootstrap interval over pairs
(interval values from the frozen report). Reliability more than quadruples from 15 to 100 cells, justifying
reliability-aware endpoints throughout.

**Constraints.** These variants test preprocessing choices on one public data release; they do not test other
screens, labs or perturbation panels.

**Purpose / why not main text.** Shows that Fig. 1 does not depend on one preprocessing choice; it is a robustness
audit of background (replication) analyses, with no new claim.

*Source data:* `data/figure_sources/paper/ext1_variants.csv`, `ext1_depth.csv`. Script:
`scripts/paper_figures/plot_ext_fig1.py`.
