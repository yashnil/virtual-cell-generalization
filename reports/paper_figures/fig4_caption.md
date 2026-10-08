# Figure 4 | Source choice strongly influences perturbation transfer

**Question.** If source quantity is not the whole story, how much does the choice of source matter?

**Data.** Preregistered N5 experiment (Outcome A). Target: K562 essential-gene screen (Replogle et al. 2022).
Sources: K562 GWPS (an independent genome-wide K562 screen from the same study; Replogle et al. 2022, figshare+
10.25452/figshare.plus.20029387, CC BY 4.0), RPE1 (same study and library), Jurkat and HepG2 (Nadig et al. 2025;
companion study from the target's lab), HCT116 and HEK293T (X-Atlas/Orion, other lab; CC BY-NC-SA 4.0, research use).
1,054 perturbations × 6,408 genes. Panels A–C fit no model.

**(A) Source → K562 transfer similarity.** C_S (called "compatibility" in the N5 protocol; a descriptive similarity,
not a biological attribution): pooled noise-corrected latent cosine between template-removed source and
target response fields, disattenuated on both sides by split-half reliable energy (1 = identical noise-free
responses). Points with 95 % paired perturbation-bootstrap intervals (2,000 resamples; the same resample for every
source; most intervals are narrower than the markers). Open marker: GWPS subsampled per perturbation to RPE1's cell
count. Right column: median per-perturbation source split-half reliability — RPE1 is the most reliable (0.53) and
lowest-C_S (0.38) source; depth-matched GWPS is less reliable than RPE1 (0.15) and has the same C_S as full GWPS.

**(B) The GWPS − RPE1 advantage under each preregistered control.** Upper block, pooled latent-cosine difference
ΔC_S = C_S(GWPS) − C_S(RPE1): primary +0.444 [0.420, 0.469]; GWPS subsampled to RPE1 depth (Adj-4) +0.439. Lower
block, per-perturbation template-removed Pearson r: raw median difference; Adj-2, source coefficient difference from
a regression with perturbation fixed effects plus source split-half reliability (linear and quadratic) and log source
magnitude, at full and matched depth (perturbation-cluster bootstrap); Adj-3, median difference over the 108
perturbations whose source reliabilities differ by ≤ 0.05. The blocks are different estimands; compare within a
block. Every interval excludes 0.

**(C) Per-perturbation distribution.** Empirical CDF of r_GWPS − r_RPE1 over all 1,054 perturbations, full depth
(solid) and depth-matched (dashed); markers at the medians. 85.8 % and 70.5 % of perturbations favour GWPS. No
subgroups were defined.

**(D) Identification limit.** Sources arranged by whether they share the target's cell line and study. The only
usable same-cell source shares the target's study and lab, so cell identity and study/lab/protocol effects are
confounded. The same-cell / different-study cell was tested in preregistered N6 with VIPerturb-seq K562 (Satija lab;
Bradu et al. 2026; Zenodo 10.5281/zenodo.18460279, CC BY 4.0), which failed its reliability gate (0.087 < 0.10;
Outcome D): uninformative, not evidence against same-cell transfer (Extended Data Fig. 4).

**Reading.** Within one study, the choice of source strongly influences perturbation transfer, and the difference is
not explained by generic measurement reliability, depth or perturbation signal strength alone. Whether it reflects
cell identity or shared study/lab/protocol cannot be determined from these data. Among non-K562 sources, sharing the
target's study and library confers no advantage (RPE1 is last; the preregistered ladder fails). The magnitude is
largely expected from Replogle et al. 2022 and Nadig et al. 2025; N5 quantifies it with noise correction. Two K562
screens from one lab share only about two-thirds of their latent response variance (0.82² ≈ 0.68).

**Analysis status.** Confirmatory (N5 primary and preregistered controls); D summarises design facts.
**Limitations.** One target and one same-cell source; no technical replicate; GWPS differs in library scale and
timepoint; the biological versus technical origin of the advantage is unresolved, and no public independent K562
dataset meeting the N6 eligibility rules was found.

*Source data:* `fig4_compatibility.csv`, `fig4_robustness.csv`, `fig4_per_perturbation.csv`. Script:
`scripts/paper_figures/plot_fig4.py`.
