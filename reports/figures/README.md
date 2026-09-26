# Figure index and result story

Nine figures, each drawn from a small figure-source table extracted from frozen
artifacts. Full provenance (source files, SHA-256, freezes, scripts, git HEAD):
[`FIGURE_MANIFEST.md`](FIGURE_MANIFEST.md).

Suitability key: **R** README hero · **P** paper main figure · **S** paper
supplement · **A** Arc engineering diagnostic.

| # | file | question | result | source report | kind | R | P | S | A |
|---|---|---|---|---|---|---|---|---|---|
| 1 | `fig1_decomposition` | How is perturbation-response energy split across template, conserved effect β, context-specific interaction γ and noise? | β holds 30.1% and γ 21.0% after noise correction, but γ is only 49.5% reproducible against 80.8% for β. | `four_context_decomposition_v1.md` | positive | ✓ | ✓ | | |
| 2 | `fig2_decomposition_robustness` | Does the split depend on one preprocessing choice? | No: across 21 variants β stays in 27.98–30.79% and corrected γ in 20.55–22.80%. | `four_context_decomposition_sensitivity.md` | positive | | | ✓ | |
| 3 | `fig3_gamma_recovery` | Does pathway aggregation reveal zero-shot interaction structure that gene level misses? | Yes in K562, RPE1 and Jurkat (p = 0.010 against 100 structure-preserving nulls; Reactome agrees); HepG2 is a clean negative. | `pathway_gamma_falsification_v1.md` | positive, context-dependent | | ✓ | | |
| 4 | `fig4_recoverability_vs_utility` | Does recovering γ improve the actual response prediction? | No: even at r = 0.78 (K562) the correction improves no context, and at the theoretically correct weight every context gets worse. | `pathway_residual_model_v2_clean_gamma.md` | negative (why pathway modelling was terminated) | ✓ | ✓ | | |
| 5 | `fig5_external_generalization` | Do network-prior predictions for never-perturbed genes survive a real context shift? | No: r 0.509 internally falls to 0.057 (arch1) and 0.018 (Feng), while direct transfer of measured perturbations holds 0.30–0.32 and is positive in 19/19 Feng lines. | `feng_multicontext_external_validation_v1.md` | negative (priors), positive (direct transfer) | ✓ | ✓ | | |
| 6 | `fig6_arc_target_support` | How much of the Arc validation panel has direct public perturbation evidence? | 86 of 300 targets (7 Tier 2, 79 Tier 1); 214 (71.3%) have none, and 73 of the 86 rest on Kaden alone. | `arc_count_space_baseline_v1.md` | descriptive | ✓ | | ✓ | ✓ |
| 7 | `fig7_source_reliability` | Is Kaden, the source behind 80 of 86 supported Arc targets and half of `m_hat`, reliable enough, and is its disagreement with arch1 just noise? | Kaden responses are weak (median reliability 0.17 vs 0.91 arch1), its main effect moderately reliable (0.76), and the arch1–Kaden main effects correlate at 0.095 against a noise ceiling of 0.86: not noise. CASE E. | `kaden_source_reliability_diagnostic_v1.md` | diagnostic | ✓ | | ✓ | ✓ |
| 8 | `fig8_count_generator_benchmark` | Which count generator best turns a predicted response into cells, scored against real held-out cells? | G1 control transport wins 5 of 6 metrics (pds_cosine 0.49 → 0.88); G2 adds capacity and loses; no generator gets the expression-error metric below 1.0. | `arc_count_space_baseline_v1.md` | positive, with a stated weakness | | | ✓ | ✓ |
| 9 | `fig9_transport_fidelity` | Does G1 deliver the mean effect it is given while keeping single-cell structure? | Mostly: realised/intended slope 0.935 (r 0.78), library sizes preserved, genes detected −3.5%, sparsity 0.587 vs 0.579; unsmoothed transport loses 22% of detected genes. | `arc_count_space_baseline_v1.md` | positive (engineering) | | | ✓ | ✓ |

## The story, in figure order

1. A real conserved effect and a real but noisy context-specific interaction
   (Fig 1), robust to preprocessing (Fig 2).
2. The interaction is partly recoverable zero-shot at pathway resolution, in some
   contexts only (Fig 3), but not accurately enough to improve prediction, so
   pathway modelling was terminated (Fig 4).
3. Predicting never-perturbed genes from priors looked good internally and
   collapsed externally; measured perturbations kept transferring (Fig 5).
4. That leaves most of the Arc panel with no direct evidence (Fig 6), and makes
   the reliability of the few sources that do exist decisive (Fig 7).
5. The count generator is not the bottleneck (Figs 8–9).

## Strongest candidates for a future paper

1. **Fig 1**: the decomposition and its reproducibility asymmetry.
2. **Fig 4**: recoverability ≠ utility; the cleanest negative result.
3. **Fig 5**: internal vs external generalization of priors vs direct transfer.
4. **Fig 3**: pathway γ recovery against structure-preserving nulls.
5. **Fig 7**: source reliability decides what a sparse panel can support
   (supplement or main, depending on the paper's framing).

No paper has been drafted.
