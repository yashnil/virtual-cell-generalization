# Paper outline (working; not a manuscript)

**Working title:** *What transfers across cellular contexts? Conserved perturbation
effects, the limits of priors, and direct-evidence transfer under a real challenge.*

Each section lists its claim, the supporting experiment, the figure or table, and the
caveat. Figures refer to `reports/figures/` unless prefixed `competition_v2/`.

## 1. Introduction

* **Claim:** zero-shot prediction of perturbation responses in new cellular contexts is
  the core problem behind "virtual cell" models. Public benchmarks and challenge
  leaderboards conflate three different sources of predictability:
  * conserved effects;
  * context-specific interactions;
  * priors for unmeasured perturbations.
* **Support:** literature (`literature_notes.md`); the Arc 2026 task definition.
* **Figure:** schematic (to be drawn).
* **Caveat:** Molina & Zhang (2026) already decompose responses this way. We
  re-derive the decomposition; we do not claim it as new.

## 2. Problem formulation

* **Claim:** the response decomposes as `δ = μ + α + β + γ` (template, context main
  effect, conserved effect, interaction). Split-half noise correction separates signal
  from measurement noise. The Arc score is normalised to a context mean-response
  baseline (0) and a replicate (1).
* **Support:** `virtual_cell.decomposition.anova` (41 mathematical tests);
  `arc2026_submission_requirements.md`.
* **Figure/table:** notation table; metric definitions.
* **Caveat:** the local `vcc2026` metric reimplementation is not the official scorer.

## 3. Perturbation-response decomposition

* **Claim:** β is substantial and reproducible (30.1 % of energy, 80.8 % reproducible).
  γ is real but half noise (21.0 %, 49.5 % reproducible). The split is robust to
  preprocessing.
* **Support:** four-context balanced design (1,264 × 6,640); 50 split-half resamples;
  21 sensitivity variants.
* **Figure:** Fig 1 (main), Fig 2 (supplement).
* **Caveat:** four contexts, one public preprocessing. This is not a Molina & Zhang
  reproduction.

## 4. Limits of zero-shot context-specific prediction

* **Claim:** γ is partly recoverable zero-shot at pathway resolution in some contexts,
  against structure-preserving nulls. Recovering it does **not** improve response
  prediction, so recoverability ≠ utility.
* **Support:** pathway falsification (100 nulls; Reactome replication); nested-LOCO
  residual models v1/v2 with a predeclared stopping rule.
* **Figure:** Fig 3 (supplement or main), Fig 4 (main).
* **Caveat:** the result holds for low-capacity models on four contexts. A
  higher-capacity or better-sourced γ model is not excluded in principle.

## 5. Limits of unseen-perturbation priors

* **Claim:** network and annotation priors for never-perturbed genes look strong
  internally and collapse externally. Direct transfer of measured perturbations stays
  positive across 19 / 19 external lines.
* **Support:** two-axis held-out design; external arch1 and Feng (19 iPSC lines),
  matched on signal strength.
* **Figure:** Fig 5 (main).
* **Caveat:** two external datasets; specific prior families (STRING / DepMap /
  pathways).

## 6. Direct atlas transfer

* **Claim:** broad direct perturbation evidence is what moves a real held-out score.
  From the sparse tiered model (V1, 86 / 300 targets) to equal-weight direct-atlas
  fusion (C1, 287 / 300), official Overall goes from −0.062 to 0.139 and PDS from 0.022
  to 0.602.
* **Support:**
  * official hidden-validation scores;
  * license-clean reimplementation, bit-exact against upstream;
  * public leave-one-atlas-out folds (H1, K562, CD4).
* **Figure:** Fig 10 plus `competition_v2/figures/c2_A_official_scorecard` (main);
  `competition_v2/figures/c1_B_public_scores` (supplement).
* **Caveat:** the backbone is AtlasShift's (attributed). V1 → C1 changes the method as
  well as the coverage. The C1 hidden score is user-reported, with its entry id pending.

## 7. Transferability depends on source context

* **Claim 1:** source agreement predicts which responses transfer. It holds in 4 / 4
  research contexts and replicates on 3 new atlases (ρ 0.34–0.50). But using it to
  shrink or reweight does not improve the predictor.
* **Claim 2:** a new direct source improves transfer only where its context transfers.
  KOLF iPSC raises coverage (287 → 298) yet helps only the pluripotent fold and hurts
  the others.
* **Support:**
  * `transferability_confidence_model_v1.md`;
  * C1b (agreement shrinkage), C3 (reliability, consensus and scale reweighting; rule O
    failed);
  * C4 (KOLF held-out transfer, rule L failed).
* **Figure:** `competition_v2/figures/c1_C_agreement_vs_transfer`, C3 figures
  (supplement), **Fig 11** (main).
* **Caveat:** three license-clean atlases make global weights unidentifiable. There is
  one new-source test, and context gating was not tested.

## 8. Arc challenge validation

* **Claim:** every successor to C1 was judged against a predeclared rule on public folds
  and rejected, for identified reasons:
  * the generator's score-positive over-calling (C2);
  * a direction bottleneck in the fused mean (C3);
  * context-specific new evidence (C4).
* **Support:** C2 / C3 / C4 reports and `competition_v2/current_champion.md`.
* **Figure:** C2 amplitude and null-DE panels, C3 transfer matrix (supplement).
* **Caveat:** validation phase only. The final contexts (D / E / F) use a different panel.

## 9. Limitations

* **The hidden truth is observed only through a leaderboard:** two submissions, one
  with its entry id pending.
* **Public folds are few** (2 VCC-scored plus 1 effect-level), and the K562 and CD4
  truths are 90–95 % noise energy, so only H1 measures direction cleanly.
* **The local scorer is a reimplementation.** Its anchors are emitted through our own
  generator, and the official anchors are not public.
* **Licensing shapes the evidence:** X-Atlas is blocked, Jurkat's license is unknown,
  and the CD4 license caveat stands.
* **No deep model was trained,** by design. Conclusions about "priors" are limited to
  the families tested.

## 10. Discussion

* **What transfers:** conserved β measured directly, from a context that transfers.
* **What does not:** zero-shot γ as a correction, annotation priors for unmeasured
  perturbations, and more coverage from a non-transferring context.
* **Practical recommendation:** judge prediction by trust (source agreement) and by
  donor-context transfer, not by target coverage.

---

## Contribution statements

### Scientific contributions (ours)

1. **Decomposition with reliability accounting.**
   * An independent four-context re-derivation of the β / γ decomposition.
   * Per-component reproducibility (β 80.8 % vs γ 49.5 %) and a 21-variant robustness
     battery.
   * The finding that γ carries three quarters of the noise.
2. **Recoverability vs utility.**
   * γ is genuinely recoverable at pathway level against structure-preserving nulls.
   * Recovering it does not improve response prediction, even at r = 0.78.
3. **Source agreement as a transferability predictor.**
   * Raw agreement beats fitted alternatives in 4 / 4 contexts.
   * It replicates on 3 independent atlases.
   * It works as a trust score but *not* as shrinkage or weighting (C1b, C3).
4. **External falsification of prior-based unseen-perturbation inference.**
   * Internal r = 0.51 collapses to ≈ 0 on two external datasets (19 lines).
   * Measured transfer stays positive in 19 / 19.
5. **Direct measurements beat generic priors.** Direct multi-context perturbation
   evidence, not biological priors, is what moves held-out performance: V1 → C1 on the
   hidden leaderboard, and direct vs prior transfer externally.
6. **Donor-context quality over coverage.** A large, reliable, GREEN new atlas (KOLF)
   increases coverage but degrades transfer outside its own context class.

### Competition-engineering contributions (ours, not scientific novelty)

* **Reimplementation:** a license-clean, bit-exact reimplementation of the AtlasShift
  backbone. *The backbone itself is not our contribution.*
* **Scoring:** a local `vcc2026` scorer, a public count-space benchmark (exact inversion
  of log-normalised scPertEval counts) and a public-fold harness.
* **Generator anatomy:** why the over-calling generator is score-optimal under
  `vcc2026` (FID and MSE-credit mechanics), and the null-DE qualification test.
* **Data pipeline:** streaming source statistics for 190 GB CSC atlases, and remote
  range-read auditing of h5ad target lists.
* **Discipline:** predeclared pass rules, freeze manifests, and code-enforced license
  gates for every phase.

### Not claimed

* **The decomposition idea** is Molina & Zhang's.
* **The direct-atlas backbone** is AtlasShift's.
* **"We predict context-specific responses":** we do not. The defensible framing is
  *conserved prediction plus a trust score*.
