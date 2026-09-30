# Project synthesis

This is a consolidated account as of 2026-09-30. Each claim is limited to what a frozen
experiment in this repository shows. Numbers are quoted from the linked reports.

## The question

Can one predict the transcriptional response to a CRISPRi perturbation in a cellular
context where it was never measured?

That is the research question (zero-shot cross-context prediction), and it is also the
Arc Virtual Cell Challenge 2026 task (six undisclosed cell lines, controls only).

## Seven findings

### 1. A conserved perturbation effect β exists and is reproducible

On a balanced four-context CRISPRi design (1,264 perturbations × 6,640 genes; K562,
RPE1, HepG2, Jurkat), after split-half noise correction:

* β holds **30.1 %** of response energy and is **80.8 %** reproducible.
* The split is stable across 21 preprocessing variants (β 28.0–30.8 %).
* Conserved source-only transfer reaches a median per-perturbation Pearson of **0.306**.

Reports: `four_context_decomposition_v1.md`, `four_context_decomposition_sensitivity.md`,
`zero_shot_recoverability_v1.md`.

*Caveat:* this is an independent re-derivation on public scPertEval data, not a
reproduction of Molina & Zhang, whose processed data are not public.

### 2. The interaction γ exists but is not useful for point prediction

* γ holds **21.0 %** of energy but is only **49.5 %** reproducible.
* It is partly recoverable zero-shot, at pathway resolution and only where a similar
  partner context exists (K562, Jurkat, RPE1; HepG2 is a clean negative). The effect
  holds against 100 structure-preserving nulls.
* A clean, β-free γ model reached r = 0.78 for γ in K562, yet improved **response**
  prediction in **no** context. At the theoretically correct weight every context got
  worse. Pathway modelling was terminated under a predeclared rule.

Reports: `pathway_gamma_falsification_v1.md`, `pathway_residual_model_v2_clean_gamma.md`.

### 3. Annotation-derived prediction of never-perturbed genes fails externally

Prior-based predictors (STRING / DepMap / pathway features) score r = 0.509 internally,
but **0.057** on arch1 and **0.018** across 19 Feng iPSC lines. On the same lines, direct
transfer of *measured* perturbations holds **0.30–0.32** and is positive in **19 / 19**
lines.

Reports: `unseen_perturbation_generalization_v1.md`,
`external_unseen_perturbation_validation_v1.md`,
`feng_multicontext_external_validation_v1.md`.

*Caveat:* two external datasets. The failure is of these priors and model classes, not
a proof that no prior can work.

### 4. Broad direct perturbation evidence dramatically improves Arc performance

| entry | model | direct evidence on the Arc panel | Overall | PDS | rank |
|---|---|---|---|---|---|
| V1 | our tiered sparse model | 86 / 300 targets | **−0.062** | 0.022 | 883 |
| C1 | equal-weight direct-atlas fusion (reimplemented AtlasShift backbone) | 287 / 300 targets | **0.139** | 0.602 | 370 / 1207 |

Both are official hidden-validation scores. C1's are user-reported and the entry id is
pending.

Reports: `arc_submission_v1_result.md`, `competition_v2/c1_official_result.md`.

*Caveat:* the backbone concept is AtlasShift's (MIT, attributed). The V1 → C1 jump
confounds coverage with the fusion method and generator, so it shows that *this* direct
atlas works, not how much coverage alone contributes. Part of the answer is in finding 7.

### 5. Source agreement predicts transferability, but acting on it does not help

* **Research track:** raw agreement among source-context responses predicts held-out
  transfer quality in all four contexts (Spearman **0.55–0.79**). It beat every fitted
  alternative (`transferability_confidence_model_v1.md`).
* **Replication:** on three new atlases (C1 folds H1, K562, CD4) the correlation is
  ρ = **0.34–0.50** (`competition_v2/license_clean_c1_v1.md`).
* **Acting on it fails:**
  * shrinking responses by agreement (C1b) lost PDS;
  * redistributing weight among donors by reliability, agreement or sign consensus (C3)
    did not improve gene directions;
  * the best candidate gave +0.003 Overall with no sign-accuracy gain
    (`competition_v2/c3_mean_response_fusion.md`).

**Conclusion:** agreement is a *trust score*, not an amplitude or a weight.

### 6. A new direct source helps only if its biological context transfers

KOLF2.1J iPSC (GREEN, 2.66 M cells, reliable signal) was added as one more equal-weight
donor:

* H1 fold (pluripotent): PDS **+0.152** (scaled), Overall +0.036.
* K562 fold: PDS −0.044, Overall −0.030.
* CD4 fold: effect PDS −0.053.

On its own, KOLF transfers to non-pluripotent contexts at cosine ≈ **0.007**, and it
agrees only with H1 (`competition_v2/c4_new_direct_evidence_audit.md`).

### 7. Coverage alone is insufficient

KOLF raises Arc coverage from 287 to 298 targets and fills 11 of the 13 unsupported ones,
and still fails the held-out transfer rule. Donor-context quality matters more than
target count (Fig. 11).

*Caveat:* this is one new source and three held-out atlases. A context-gated use
(pluripotent targets only) was not tested.

## Competition-engineering findings (not scientific claims)

* **The generator is not the lever** (C2).
  * The frozen generator over-calls DE: 3,837 null DE genes on H1.
  * Under `vcc2026` scoring, that over-calling *raises* FID, which rewards calls at
    ~55 % precision, and the MSE noise credit favours its averaged template.
  * Realistic real-donor generators matched cell heterogeneity and raised PDS, yet lost
    Overall (`competition_v2/c2_expression_de_calibration.md`).
* **The mean is not "too small" where it is fused** (C3). It is mostly orthogonal to the
  truth (cosine 0.03–0.085). Direction is the bottleneck.
* **With three atlases, global donor weights are not identifiable** (C3).

## Current state

* **Champion:** C1 (`competition_v2/current_champion.md`).
* **X-Atlas permission:** PENDING.
* **C5 (C1 + X-Atlas):** predeclared and gated in code (`competition_v2/c5_predeclaration.md`).
