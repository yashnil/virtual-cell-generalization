# Extended Data Figure 4 | The independent-lab K562 test (N6) failed its preregistered reliability gate

N6 asked whether an independent-lab K562 screen (VIPerturb-seq, Satija lab; 637 panel perturbations × 6,083 genes)
keeps the within-study compatibility of Fig. 4. Its preregistered gate decides whether any compatibility value can
be interpreted.

**(A) Gate item d1.** Pooled split-half reliability (mean over 5 repeats of the pooled Pearson of centred split
halves). VIPerturb-seq: 0.0866 (frozen gate statistic) < 0.10 → Outcome D. Open squares: the same statistic for
K562 GWPS subsampled to VIPerturb's per-perturbation cell counts (0.109) and at full depth (0.341), from an
exploratory, read-only post-hoc diagnostic (N6 results §4; not part of the preregistered analysis). Depth accounts
for most of the gap; at equal depth VIPerturb retains ~20 % less reproducible energy.

**(B) Per-perturbation reliability.** ECDF over the 637 panel perturbations; legend gives the fraction of
perturbations with reliability ≥ 0.10 (VIPerturb 22 %, GWPS at VIPerturb depth 27 %, GWPS full depth 82 %).

**(C) Positive control (gate item d3).** C_S of K562 GWPS at full depth and at VIPerturb depth, with 95 % paired
bootstrap intervals; grey band: the preregistered tolerance (±0.10). The shift (0.007) shows that the endpoint
itself is stable at VIPerturb's cell counts, so low depth alone would not have prevented an informative test.

**(D) Gate checklist.** One failing item yields Outcome D.

**What is deliberately not shown.** All VIPerturb compatibility values and contrasts (N6 results §2) are omitted.
Under the frozen protocol they are non-interpretable and cannot be read as evidence for or against same-cell
transfer across studies. The same-cell / different-study question remains unanswered, not answered negatively.

**Purpose / why not main text.** Documents why N6 is Outcome D (uninformative) rather than a biological result.
VIPerturb-seq: Bradu et al., bioRxiv 2026; Zenodo 10.5281/zenodo.18460279; CC BY 4.0 (converted from Seurat;
pseudobulk statistics derived; `reports/viperturb_conversion_audit.md`).

*Source data:* `data/figure_sources/paper/ext4_*.csv` (from `outputs/n6/{source_gate,n6_decision}.json`,
`per_perturbation_terms.csv`; the two exploratory pooled values are transcribed from the frozen report and flagged
in the provenance sidecar). Script: `scripts/paper_figures/plot_ext_fig4.py`.
