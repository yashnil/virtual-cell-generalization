# Extended Data Figure 2 | The γ⊥ recovery curve is not produced by leakage, shared control noise or chance

Preregistered failure checks for N3-A (six contexts) and the confirmatory N1 run.

**(A) Anchor-permutation null (F3a).** E4 γ⊥ recovery with the correct anchor ↔ perturbation pairing (diamonds;
200 draws) and with the anchors' measured responses permuted among the anchors before fitting (×; 50 draws), at
k = 20 (open) and k = 100 (filled). Median and 2.5–97.5 % draw interval. Permuted medians lie between −0.0014 and
+0.0007 in all 12 cells: the γ⊥ gain requires the correct perturbation identity.

**(B) Shared-control diagnostic (F3b).** Median over 50 paired draws of (score with one control mean shared by
anchors and truth) − (score with disjoint control parts), E4, k = 5 (open), 20, 100 (light), per context.
Template-level M0 would be inflated by up to +0.027 if controls were shared; the γ⊥ endpoint M3 is identical to
numerical precision in every cell, because a control error is constant over perturbations and is removed by
centring.

**(C) Too few anchors are harmful.** Full-response gain over zero-shot (ΔM0) when the target template is estimated
from k anchors (E1), per draw paired with the same cell-split repeat's zero-shot value; thin lines: context medians
over 200 draws; thick: median across contexts; band: range across contexts (descriptive). With one anchor its own
response is mistaken for the template (median ΔM0 −0.47 to −2.51).

**(D) Confirmatory N1 and the N3-A replication agree.** γ⊥ recovery at k = 20 as a fraction of the k_ref value
(ratio of draw medians, the preregistered C3 quantity), with the 2.5–97.5 % interval of the per-draw ratio
(k_ref paired by repeat; descriptive). N1: 3 sources, 1,264 perturbations, k_ref = 885; N3-A: 5 sources,
1,062 perturbations, k_ref = 743. Dashed line: C3 bar 0.25. Only RPE1 clears it in either run.

**Also passed (not plotted; tests in the repository):** F1 noise-replacement leakage test (predictions
bit-identical when every non-anchor target value is replaced by noise), F2 source-only zero-shot test, F4 cell and
control disjointness (asserted at build), F5 synthetic recovery (planted γ⊥ recovered at M3 = 0.92 by E4 at
k = 100; no γ → ≤ 0), F6 decomposition identity.

**Purpose / why not main text.** Shows that the Fig. 2 γ⊥ curve is not produced by leakage, shared-control noise or
chance pairing, and that N1 and N3-A agree. Validity checks are methodological rigour, not findings (cf. Miller et
al. 2026 on experimental controls).

*Source data:* `data/figure_sources/paper/ext2_*.csv`. Script: `scripts/paper_figures/plot_ext_fig2.py`.
