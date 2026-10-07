# Publication novelty audit

Date: 2026-10-06 · Branch `research/c6-uncertainty-aware-transfer` (HEAD `47bca2b` + uncommitted C6 files)
Scope: audit and research planning only. **No modelling code was changed and no experiment was run for this
document**, apart from one read-only count of perturbation identifiers (§2.2, marked).

Evidence tags:
* **[repo]**: a frozen report in this repository. The file is named.
* **[R]**: a primary source read in this session (abstract or text).
* **[A]**: an abstract or search summary only, not verified in full.
* **[M]**: taken from `reports/literature_notes.md` (verified on 2026-09-06, not re-read today).

---

## 1. Executive summary

**Where the project is.** The repository holds a careful, reliability-corrected, leakage-tested study of
zero-shot CRISPRi response transfer across four public cell lines, plus an Arc 2026 competition track. Its
*solid* results are:

1. An independent re-derivation of the Molina & Zhang decomposition.
   * β is 30 % of energy and 81 % reproducible.
   * γ is 21 % of energy and 50 % reproducible.
   * This is a **replication**, not a discovery.
2. Zero-shot conserved transfer captures direction.
   * Median r is 0.31, about half of the noise ceiling.
   * It explains *negative* energy unless shrunk by a source-fitted scalar (about 0.44).
3. Raw source–source agreement ranks per-perturbation transfer quality.
   * Risk-coverage curves are monotone in 4/4 LOCO folds.
   * The ranking replicated on 3 Arc-like atlases (ρ 0.34–0.50).
   * No fitted model beat it.
4. γ is partially recoverable zero-shot at pathway level, but only with a similar partner context. Using it
   as a correction never improved prediction.
5. Annotation priors for never-measured perturbations collapse externally (arch1; 19 Feng iPSC lines).

Six competition successors (C1b, C2, C3, C4, C6, and the γ models) failed predeclared rules. C1 scored 0.139
officially (rank 370/1207, user-reported).

**Novelty verdict.** As the repository stands, **no headline claim is publication-level novel**:

* The decomposition and "β transfers, γ does not zero-shot" belong to Molina & Zhang (2026).
* "Simple baselines win" and "metrics are confounded" are crowded (Ahlmann-Eltze, Systema, PerturBench,
  VCC 2025, and others).
* "Basal similarity is insufficient" is pre-empted by Molina & Zhang, Li et al. and Mechanisms Matter.
* Uncertainty for perturbation models exists (PRESCRIBE, GPerturb, ConfPert).
* Arc's PIE (posted 2026-10-05) now claims state of the art on unseen contexts on the same Replogle–Nadig data.

The **closest-to-novel asset** is source agreement as a *prospective* per-perturbation transferability score
with risk-coverage evaluation. Even that reads to a skeptical reviewer as random-effects heterogeneity (an
I²-like statistic), unless it is shown to do more than detect weak perturbations.

**The real open question.** Every serious paper (Molina & Zhang, State, PIE's setting, PerturbMap) finds that
γ becomes learnable once *some* target-context perturbations are measured. Nobody has characterised the
*budget curve* for a new CRISPRi cell line:

* how many target perturbations are needed;
* which response component (template, scale, γ) each additional measurement buys;
* whether choosing them well beats choosing them at random.

The closest work:

* State and Molina & Zhang use a fixed 30 %.
* Pan et al. report about 100 perturbations, plateauing.
* PerturbMap uses about 128 random anchors on melanoma Perturb-CITE-seq conditions.

**Recommendation.**

* **Primary:** *the calibration budget of a new cellular context.* Measure a decomposition-explained
  learning curve from 0 to about 200 target perturbations. Test whether an informed choice of anchor
  perturbations beats random. Use the zero-shot trust score as the k = 0 endpoint. Replicate on ≥ 6
  contexts, which requires the research-track (non-commercial) use of X-Atlas, already on disk.
* **Backup:** *how much of "context-specificity" is lab/assay shift rather than biology.* Build a shift-type
  ladder: replicate → same line, other screen → other donor, same lineage → other lineage.

**Highest-information next experiment (N1).** The k-anchor learning curve on the frozen four-context tensor.
It is CPU-only and takes minutes. It can kill the primary direction in one run if γ-space gains appear only at
k ≳ 30 % of the panel, which is the regime already reported by State and Molina & Zhang.

**Biggest risk to the whole project.** Every scientific conclusion rests on **four essential-gene screens from
two labs**. Context, lab, library and assay are confounded. The Arc 2026 contexts are further apart (basal r
0.75–0.84 vs 0.89–0.93) and use a different perturbation regime (> 10,000 perturbations per line, ≥ 80 %
knockdown). Results may not travel. On the Arc-like atlases, even "good" transfer is cosine 0.03–0.09.

---

## 2. Current project state (Phase 1)

### 2.1 A. Scientific question

| item | exact current definition | source |
|---|---|---|
| prediction target | pseudobulk response δ[c,p] = mean log1p(CP10K) of perturbed cells − control mean, per gene (research track). Arc: raw counts, 400 cells × 18,533 genes per (context, target) | [repo] `four_context_decomposition_v1.md` §4; `arc2026_submission_requirements.md` |
| available at prediction time | (i) the held-out context's **control cells only**; (ii) perturbation identity; (iii) source-context responses; (iv) external priors free of target outcomes (STRING, DepMap, MSigDB, GENCODE) | [repo] `zero_shot_recoverability_v1.md` §1 (leakage contract, tested) |
| held out | every perturbation response of one context (LOCO). Also forbidden: the four-context β and γ, because they contain the target | same |
| zero-shot definition | no measured perturbation of any kind in the target context; only its controls | same |
| context setup | research: K562, RPE1 (Replogle 2022), HepG2, Jurkat (Nadig 2025). Competition: VCC-2025 H1, K562 GWPS, CD4 GWCD4i DE, KOLF2.1J (C4) | reports above |
| perturbation setup | single-gene CRISPRi; research panel = 1,264 essential-library perturbations shared by all four lines | `data/splits/four_context_v1/` |
| metrics | research: per-perturbation Pearson, energy explained, reliability-normalised r/√ρ, reliable residual fraction D = ⟨h1−A,h2−A⟩/⟨h1,h2⟩. Competition: vcc2026 six-member mean (PDS, MSE, NMAE, FID, REACH, JAC), scaled 0 = context mean-response baseline, 1 = replicate | `transferability_foundations_v1.md` §D; `literature_notes.md` §1 |
| Arc relationship | Arc 2026 is the same task (six undisclosed lines; controls plus a 300-target list). Validation A/B/C is live; final D/E/F is released 2026-10-22 with a different panel; due 2026-11-05 | `literature_notes.md` §1 [M] |

### 2.2 B. Dataset and evaluation setup

| dataset | context | role | scale | status |
|---|---|---|---|---|
| scPertEval `replogle22k562`, `replogle22rpe1`, `nadig25hepg2`, `nadig25jurkat` | 4 lines | research core (balanced 4 × 1,264 × 6,640) | 133k–309k cells each; median 57–139 cells/pert | local, hashed |
| `arch1` (VCC 2025 H1 subset) | H1 hESC | external unseen-perturbation test | 150 perts, ~1,045 cells/pert, reliability 0.84–0.95 | local |
| `kaden25rpe1` | RPE1 (second screen) | study-shift test | 1,836 perts; reliability 0.12–0.17 (**dead positive control**) | local |
| `wessels23` | — | audited, unused | 157 perts | local |
| Feng et al. 2026 (Cell Genomics) | **19 iPSC lines (donors)** | external unseen-perturbation test | log fold changes only; counts not downloaded | local LFC |
| VCC 2025 H1 (Training/Validation/Test) | H1 | competition source | ~491k cells | local |
| K562 GWPS (`K562_gwps_raw_singlecell_01`) | K562 | competition source | genome-wide | local |
| GWCD4i DE stats (Marson 2025) | primary CD4 T | competition source | DE summary statistics only (no cells, no split-half) | local |
| KOLF2.1J pan-genome | iPSC | C4 candidate | 2.66 M cells, 189 GB | local |
| **X-Atlas/Orion** | **HCT116, HEK293T** | competition **blocked** (prize use needs Xaira permission) | genome-wide FiCS; ~18,300 targets each | **local**, CC BY-NC-SA 4.0 |

**Read-only identifier count (done for this audit, no expression values read).** Perturbations with ≥ 30
cells:

* The 4-context panel ∩ HCT116 = **1,084** perturbations.
* The 4-context panel ∩ HCT116 ∩ HEK293T = **1,062** perturbations.

So a balanced **6-context × ~1,062-perturbation** research design is feasible from local data. The precomputed
KOLF statistics cover only the 437 competition-relevant targets; a full KOLF intersection needs a recount from
its obs.

**Licensing note (to resolve, not to assume).** CC BY-NC-SA permits non-commercial academic research.
X-Atlas is blocked only for the prize-bearing entry. A research-track use still needs:

* a written license memo in `data_license_register.md`;
* separation from all competition candidates, which is already enforced in code;
* attention to ShareAlike for any redistributed derived matrices.

**Splits and protections**

* LOCO: four folds, frozen in `data/splits/loco_v1/` with hashes. Competition folds are leave-one-atlas-out
  (H1, K562, CD4).
* The leakage test replaces the whole target row with noise and requires bit-identical predictions
  (`tests/test_loco.py`, `tests/test_foundations.py`).
* The shared perturbation and gene sets are identifier-only intersections. There is no HVG or
  outcome-dependent filtering.
* Every phase has a SHA-256-frozen predeclaration (`data/provenance/**/c*_predeclaration_digest.txt`,
  `protocol_freeze.txt`). Freeze manifests are re-verified at the start of each phase.
* Licensing is enforced in code (`competition_v2/licensing.py`).

**Noise and reliability machinery**

* 50 split-half resamples of the decomposition.
* An independent-control-split robustness check.
* Disattenuation r/√ρ, validated on planted signals.
* The unbiased reliable-energy estimator ⟨h1,h2⟩ and the residual fraction D.
* A stability rule (minimum signal energy 1.967), derived by simulation before scoring.

### 2.3 C. Current empirical findings

| # | Experiment (report) | Hypothesis | Method | Result | Pass/Fail | Interpretation |
|---|---|---|---|---|---|---|
| 1 | Four-context decomposition v1 | β and γ both materially non-zero | Two-way ANOVA on δ, 50 split halves | template 20.3 %, β **30.1 %**, γ **21.1 %**, noise 28.6 %; β 80.8 % and γ 49.5 % reproducible; γ carries 75 % of all noise | PASS | Replicates Molina & Zhang's split (their 27.8/29.4/23.5/19.3) on an independent pipeline. **Replication** |
| 2 | Decomposition sensitivity | Split robust to preprocessing | 21 variants (control split, gene space, aggregation order, seeds), controlled depth | β 28.0–30.8 %, γ ≥ 20.6 %; reliability rises with depth (0.097 → 0.419) | PASS | Robust within one data family |
| 3 | Zero-shot recoverability | Source-only transfer recovers held-out responses | LOCO; source mean, basal-weighted, nearest-basal | median r 0.306 (0.28–0.36); **energy explained negative in 3/4 folds**; averaging beats picking; nearest-basal is worst | Mixed | Direction transfers; magnitude and template do not. Basal similarity does not help pick a source |
| 4 | γ zero-shot recovery (gene level) | γ predictable from sources | r(γ, γ̂) by fold | K562 0.19, Jurkat 0.22, RPE1 0.00, HepG2 0.03 | Partial | γ shared only by the K562–Jurkat pair (excess over the −1/3 null +0.196; other pairs ≈ 0). Basal similarity vs γ-sharing Spearman 0.71, **n = 6 pairs, one point drives it** |
| 5 | Foundations A: template from basal | Basal deviation encodes α | correlation, subspace energy, LOCO estimators | r ≈ 0, sign-flipping; 0.15–8.2 % of α in the basal span | FAIL | Same conclusion as Molina & Zhang (controls carry little of the template) |
| 6 | Foundations A: scale | Source-fitted shrinkage fixes energy | inner leave-one-source-out scalar | shrink 0.43–0.47 in every fold; energy −0.47 → −0.21 (K562), positive in HepG2 | PASS | Variance-optimal shrinkage follows from unpredictable γ. A simple, robust fix |
| 7 | Pathway γ falsification | Pathway γ beats matched random | Hallmark/Reactome vs 100 geometry-preserving nulls | normalised ~2× gene level in K562/Jurkat/RPE1 (p = 0.01 floor); HepG2 null; Jurkat wholly K562-dependent (−0.195 without K562) | PASS 3/4 | γ structure is biological but partner-dependent |
| 8 | Pathway residual model v1/v2 | Predicted γ improves response prediction | nested LOCO, clean β-free target | γ̂ r = 0.78 (K562) but Δr ≤ 0 in 4/4; λ_theory hurts every fold; Hallmark = random | **FAIL** (terminated by rule) | "Recoverable" is not "useful". The error injected exceeds the signal at r ≈ 0.6–0.8 |
| 9 | Transferability foundations B | Source agreement predicts transfer | Spearman vs reliability-normalised success, confound partials | 0.62–0.84 per fold; joint partial 0.14–0.55 | PASS (attenuated) | Partly a weak/noisy-perturbation detector; residual signal is real but modest, weakest in HepG2 |
| 10 | Transferability confidence v1 | Learned confidence beats raw agreement | M0 raw vs M1 isotonic vs M2 ridge, nested LOCO; risk-coverage; calibration quintiles | M0 ≥ M2 on ranking 4/4, risk-coverage 16/20; monotone curves; D falls 18–41 % at 10 % coverage | PASS for M0; learned FAIL | Simple trust score. **Trust ≠ error magnitude:** least-confident selection captures *half* of random's error; magnitude ranking captures 2× random |
| 11 | Unseen-perturbation priors (internal → arch1) | STRING/DepMap kNN predicts unseen β | two-axis held-out | internal r 0.51 → arch1 0.06 (energy worse than zero) | FAIL externally | Internal benchmark measures within-design interpolation |
| 12 | Kaden external | Separate study shift from context shift | same RPE1 line, other screen | positive control dead (r −0.03; reliability 0.12) | Uninformative | Study-vs-context question **still open** |
| 13 | Feng 19 iPSC lines | Prior failure is real | per-line unseen-perturbation vs measured transfer | priors 0.007 median; measured 0.128/line, 0.320 pooled, 19/19 positive; between-line variation tracks reliability (r 0.956) | PASS (confirms #11) | Two external datasets. Context distance barely varies within Feng |
| 14 | Kaden source-reliability diagnostic | arch1/Kaden disagreement is noise | 50 block split halves | main effects reliable (0.76) yet correlate 0.095; Kaden ≠ Replogle RPE1 more than lines differ | CASE E (inconclusive) | Suggests a large lab/screen effect. **Suggestive only** |
| 15 | Arc V1 (official) | Tiered sparse model scores > 0 | 86/300 targets with direct evidence | Overall −0.062, rank 883 | FAIL | FID accounts for 94 % of the deficit |
| 16 | C1 license-clean atlas (official, user-reported) | Broad direct evidence moves the hidden score | equal-weight fusion of H1/K562/CD4 (AtlasShift backbone, attributed) | Overall **0.139**, PDS 0.602, rank 370/1207 | PASS | Engineering. Confounds coverage, method and generator |
| 17 | C1b agreement shrinkage | Shrinking low-agreement targets helps | nested folds | lower PDS 0/3 | FAIL | Agreement is a trust score, not an amplitude |
| 18 | C1 agreement replication | #10 generalises to new atlases | ρ(agreement, transfer cosine) | 0.50 / 0.49 / 0.34; high-tercile cosine only 0.05–0.10 | PASS | **Replicated across independent atlases**, but the absolute transfer is tiny |
| 19 | C2 generator/DE calibration | Generator realism raises score | G0–G3, amplitude | realistic generators raise PDS, lose Overall; scorer rewards DE over-calling | FAIL | Scorer property; no scientific content |
| 20 | C3 source fusion | Reweighting donors fixes direction | reliability, consensus, scale weights | best +0.003 Overall; sign accuracy flat; transfer cosine 0.03–0.085 | FAIL | Direction is the bottleneck; 3 atlases cannot identify weights |
| 21 | C4 KOLF new source | More coverage helps | KOLF as an extra donor | H1 fold +0.036; K562 −0.030; CD4 −0.053; KOLF→non-pluripotent cosine 0.007 | FAIL | **A lineage-matched source helps; a mismatched one hurts.** One source |
| 22 | C6 uncertainty-aware fusion | The more precise source wins conflicts | delta-method σ², IVW/RE/EB | P(lower-σ correct) 0.47–0.51; no candidate beats C1 | FAIL | Disagreement is not measurement noise (H1–CD4: τ² ≈ 4× σ²). Precision ≠ transportability |
| 23 | C6 spectrum | Responses are low-rank | noise-ceiling spectrum | H1: 20 components, 67 % of energy; K562: 1 | Mixed | Deep sources are low-rank; noisy ones are not |

**Established** (replicated within the repository or across datasets):

* The β/γ split and its robustness. This replicates Molina & Zhang.
* γ is about half noise.
* Conserved transfer is directionally positive but magnitude-miscalibrated.
* Scale shrinkage of about 0.44 is consistent across folds.
* Source agreement ranks transfer quality: 4/4 research folds plus 3/3 atlas folds.
* Raw agreement ≥ learned confidence.
* Annotation-prior unseen-perturbation prediction collapses externally (2 datasets).

**Suggestive** (single dataset, n ≤ 6 pairs, or one source):

* γ is shared only between similar partners.
* Pathway γ is biological, not geometric.
* Lineage-matched sources help and mismatched ones hurt (KOLF).
* Lab or screen effects may rival lineage effects (Kaden).
* Trust and error magnitude behave as opposite objectives.
* Source disagreement is context rather than noise (one pair).

**Failed hypotheses:**

* γ correction improves prediction (twice).
* Basal state recovers the template.
* Basal-weighted source choice beats averaging.
* Learned confidence beats raw agreement.
* Agreement shrinkage helps.
* Reliability, consensus or uncertainty weighting helps (C3, C6).
* A new non-matched source helps (C4).
* Generator realism helps the score (C2).
* Prior-based unseen-perturbation prediction transfers externally.

**Not actually tested:**

* Any **few-shot / k-anchor** curve in a new context.
* Active vs random choice of target measurements.
* More than 4 contexts in the research design (X-Atlas is local).
* Non-essential perturbations (the panel is essential-library only).
* A clean study-vs-context separation (Kaden failed as a test). K562 essential vs K562 GWPS is available and
  unused for this.
* Donor-genetic γ (Feng counts not downloaded).
* Foundation-model embeddings.
* Molina & Zhang's own models, State, or PIE under our splits.
* Prospective transfer to Arc D/E/F.
* Context-gated sources.
* Cross-perturbation denoising.
* Single-cell distributional transfer.

### 2.4 D. Strongest supported claim and the tempting unsupported one

**Strongest supported claim.** *In CRISPRi Perturb-seq across distinct human cell lines, a perturbation's
zero-shot cross-context transfer quality can be ranked in advance from source data alone. Raw agreement among
source-context responses gives monotone risk-coverage curves in every held-out context: 4/4 public cell lines
and 3/3 independent atlases. No fitted confidence model improved on it. The conserved predictor it qualifies is
directionally right (about half the noise ceiling) but needs about 0.44× shrinkage to be a sensible point
prediction.*

**Tempting, but not supported:**

* *"We can predict when transfer will succeed."*
  * After controlling for source magnitude and reliability, the partial association is 0.14–0.55.
  * Much of agreement's power is detecting perturbations with little reproducible signal anywhere.
  * In an exchangeable random-effects world, agreement predicts transfer *by construction*. Nothing yet shows
    it predicts *non-exchangeable* failures (the partner and lineage effects).
* *"Context-specific effects are fundamentally unpredictable."*
  * Only four contexts were tested.
  * Molina & Zhang, and our own pathway result, show γ *is* learnable with partner or target data.
* *"Basal state is insufficient."*
  * The evidence is 6 context pairs, and it is already Molina & Zhang's claim.
* *"Coverage is not transferability."*
  * The evidence is one new source and three folds.
* *"Our approach generalises to Arc's unseen lines."*
  * There is one hidden score, user-reported, from an attributed backbone.

---

## 3. Novelty assessment

| level | assessment |
|---|---|
| competition novelty | low. C1 is AtlasShift's backbone; rank 370/1207 |
| engineering novelty | moderate. The frozen predeclaration, licensing-in-code and reliability tooling are exemplary but not publishable alone |
| methodological novelty | **low to moderate.** The reliable-energy residual fraction D, the disattenuated LOCO protocol, the pathway null battery and risk-coverage for unseen contexts are useful. Each is a standard tool (split-half reliability, disattenuation, selective prediction) applied carefully |
| biological novelty | **low.** Conserved housekeeping (POLRMT, TFAM) vs context-specific responses (BCR, GATA1-like) is already described in Nadig 2025, Jiang 2025 and Pan 2026 |
| publication-level novelty | **not yet.** The current package reads as "careful replication of Molina & Zhang plus negative results plus a heterogeneity-based trust score" |

---

## 4. Literature comparison (Phase 2)

| Paper | Core contribution | Already answers | Leaves unresolved | Overlap with us | Opportunity |
|---|---|---|---|---|---|
| **Molina & Zhang 2026**, bioRxiv 10.64898/2026.07.24.740459 [R abstract; M full text] | μ+α+β+γ decomposition with split-half noise; component-aligned Ridge/MLP on DepMap; 4 lines | β transfers; γ **not predicted zero-shot by any model** (linear, MLP, STATE, MORPH, even with target controls); γ learnable when 30 % of target perturbations are measured | Budget curve below 30 %; which perturbations to measure; per-perturbation triage; > 4 contexts | **Very high.** Our decomposition and γ-negative replicate theirs | Their 30 % point is one point on a curve nobody has drawn |
| **PIE** (Verma, …, Roohani; Arc), bioRxiv 10.64898/2026.10.02.756297, 2026-10-05 [R abstract] | Perceiver-style model over knowledge, basal expression and observed responses; predicts DEGs/LFC per (context, perturbation) | SOTA on Replogle–Nadig unseen context, unseen perturbation and both; DEG AUPRC 1.2–3.2× best baseline; cross-dataset | Per-context difficulty, uncertainty and few-shot are not in the abstract; whether γ (vs β) is what improved is unknown | High: same data, same setting | Must be run under our reliability-normalised, component-attributed evaluation. Does PIE's gain live in β or γ? A cheap, sharp side analysis |
| **State** (Adduri et al., Cell 2026) [M; R partial] | Set transformer + SE embedding; ~70 contexts | Zero-shot context: matches the perturbation-mean baseline on genetic data; few-shot at a fixed 30 % | No budget curve; no selection policy | Moderate | Same: budget curve |
| **TxPert** (Wenkel et al., arXiv 2505.14919 → Nat Biotech 2026) [R] | Graph-prior model; leave-one-line-out over the same 4 lines | Beats a mean/additive "general baseline" in 4/4 | Nearest-line baseline only partly reported (unverified snippet: 2/4); no uncertainty | Moderate | Our reliability-normalised ladder is a fair-comparison tool |
| **Nadig et al. 2025** (Nat Genet; TRADE) [R] | HepG2/Jurkat screens; transcriptome-wide effect estimation | Cross-line effect correlations (K562–Jurkat 0.74, HepG2–RPE1 0.75, K562–RPE1 0.40); 56 % consistent vs 44 % pair-specific | No prospective predictor; lab and context confounded | High on biology | Our partner result re-states theirs. Do not claim it |
| **Li et al.** (bioRxiv 10.1101/2024.12.23.630036 v3) [R] | Benchmark of response models including cell-type transfer | "Transferability not determined by basal-state distance alone" | Stimulation/drug data only; descriptive | Moderate | Pre-empts "basal is insufficient" as a headline |
| **Mechanisms Matter** (Qi & Chapfuwa, ICML GenBio 2026; bioRxiv 10.64898/2026.05.08.723625) [R abstract] | Transportability framing + semi-synthetic simulator | Transfer is governed by shared mechanisms, not distributional similarity | No prospective per-perturbation predictor; no real-data budget | Moderate | Pre-empts the transport framing. Cite it; do not claim it |
| **Virtual Cells Need Context, Not Just Scale** (Dibaeinia et al., bioRxiv 10.64898/2026.02.04.703804) [M] | Position paper; context diversity > cell count (CD4 donor × time) | Context coverage matters | Single dataset; no design policy for a new context | Low–moderate | Our context-level design question extends it |
| **PerturbMap** (Cui, Liu, Sun, arXiv 2607.28090) [R abstract; agent read of HTML] | Cross-context transfer via paired anchor perturbations in the recipient context | Anchors (~128) enable transfer; capping anchors hurts; route reliability varies 6.2× | **Random anchors only**; no zero-shot endpoint; melanoma Perturb-CITE-seq conditions (+ Jiang 6 lines per agent, unverified); +4.1 % MSE | **Closest to the primary direction** | Shows the question is live and unsolved: small-k regime, CRISPRi cell lines, selection policy, decomposition attribution |
| **Pan, Saunders, Replogle, Weissman, Zhuang 2026** (bioRxiv 10.64898/2026.05.16.725005) [M] | Mixture-of-experts manifold; K562/RPE1/hESC | Conserved vs context-specific programs; few-shot r ≈ 0.4–0.5 at 100 perturbations, plateau; zero-shot < baseline | No selection policy; no attribution of what anchors buy | Moderate | Gives an external number to compare our curve against |
| **ConfPert** (Alwani & Wang, ICML 2026 workshop) [R abstract] | Model-agnostic conformal coverage on distributional discrepancies; one cross-cell-line split | Finite-sample coverage | No risk-coverage or abstention; exchangeability under context shift is not addressed | Moderate for selective prediction | Selective prediction alone is no longer clean novelty |
| **PRESCRIBE** (NeurIPS 2025) [R abstract]; **GPerturb** (Nat Commun 2025) [A] | Evidential / GP uncertainty | Filtering low-confidence predictions helps (unseen genes) | Unseen contexts | Moderate | "Uncertainty for perturbation prediction" is crowded |
| **IterPert** (Huang et al., RECOMB 2024) [A]; **Panagopoulos et al.** arXiv 2503.14571 [R abstract] | Active / one-shot selection of perturbations | Active selection saves about 3× *within* a context (unseen perturbations) | Selection for calibrating a **new context** | Moderate | The cross-context version is the gap. Their within-context priors are baselines |
| **MapPFN** (ICLR 2026, arXiv 2601.21092) [R abstract] | In-context adaptation to new contexts | Adapts from interventional evidence at inference | No budget curve; synthetic pretraining | Low–moderate | A competitor for the k > 0 regime |
| **Wang, Kuipers, …, Beerenwinkel 2026** (bioRxiv 10.64898/2026.08.11.744177) [M; could not re-find today] | Per-perturbation reliability ρ and specificity φ across 29 datasets | 65 % unreliable, 11 % shared, 24 % specific (retrospective labels) | Prospective prediction of specificity in an unseen context | Moderate on the trust score | Our trust score must beat "reliability-only" triage. Already partly shown (partial 0.14–0.55) |
| **Ahlmann-Eltze 2025; Systema; PerturBench; PertEval-scFM; Csendes; Kedzierska; Mao 2026** [A/M] | Simple baselines ≥ deep / foundation models; metric confounds | Crowded consensus, with a counter-claim (Shift Bioscience, Nat Biotech 2026, calibrated metrics) | — | High if we headline it | Avoid as a claim; use as the evaluation hygiene section |
| **Jiang et al. 2025** (Nat Cell Biol; Satija) [R] | 6-line pathway-regulator screens | IFNγ regulators conserved; TGFβ/insulin regulators line-specific | Prospective transfer | Low (different modality) | A possible 3rd dataset for replication (if counts are available) |
| **Arc VCC 2026** (Cell commentary; Arc data post) [R] | Six lines (stem, immortalised healthy, cancer); > 10,000 perturbations per line; ≥ 80 % knockdown selection | — | Whether perturbed data will be released (not stated) | Direct | If released: the ideal prospective, independent replication |
| **X-Cell / X-Atlas/Orion** (Xaira 2025–26) [A] | Large diffusion model; HCT116/HEK293T genome-wide CRISPRi | Zero-shot to new primary cells (claimed) | Independent evaluation | Data overlap | X-Atlas is our route to 6 contexts |

Verification gaps:

* I could not read the full texts of Molina & Zhang (re-read today), PIE or Mechanisms Matter. **Read PIE in
  full before any framing is fixed:** it is from Arc and posted yesterday.
* COMPASS (Liang & Singh 2026), cited in our notes as using 6 CRISPRi lines including HCT116/HEK293T, was
  not re-found by today's search. If it exists, it overlaps any 6-context design.
* "scPerturBench Nature Methods" was not re-found today. Our notes cite a DOI; re-verify before citing.

---

## 5. Crowded directions to avoid as headline claims

1. **"Simple or linear baselines match or beat deep and foundation models."** Crowded and contested.
2. **"Metrics are confounded / systematic variation."** Systema, Diversity by Design, cell-eval, Heidari.
3. **"Decompose responses into β and γ; β transfers, γ does not zero-shot."** Molina & Zhang own it. We
   replicate it.
4. **"Basal state / similarity is insufficient for cross-context transfer."** Molina & Zhang, Li et al.,
   Mechanisms Matter.
5. **"Transportability framing."** Mechanisms Matter; the context-not-scale position paper.
6. **Another knowledge-prior model for unseen contexts.** PIE, TxPert, State, X-Cell, MORPH, C3TL, AdaPert,
   CellFlow. We should not build one.
7. **"Uncertainty for perturbation prediction" in general.** PRESCRIBE, GPerturb, ConfPert.
8. **"Foundation-model embeddings do not help."** PertEval-scFM, Kedzierska, Zinchenko.
9. **Descriptive taxonomy of conserved vs context-specific perturbations.** Nadig, Jiang, Pan, and Wang
   reliability/specificity labels.
10. **Annotation priors for unseen perturbations fail out of distribution.** Our external result is clean
    but confirms a crowded message. It works as a supporting figure, not a headline.

---

## 6. Open scientific questions (Phase 3 brainstorm, 13 candidates)

| ID | Direction | One-line question |
|---|---|---|
| D1 | Prospective transferability | Can source-only features rank per-(perturbation, context) transfer better than "is this perturbation reliable/strong", under LOCO? |
| D2 | Selective prediction | Can a zero-shot model abstain so that the accepted subset meets a stated accuracy, with valid coverage under context shift? |
| D3 | Information limits | Are there context pairs with near-identical basal states and divergent responses, i.e. an empirical non-identifiability of γ from controls? |
| **D4** | **Calibration budget** | **How many target-context perturbations does a new cell line need, what does each buy (template → scale → γ), and does informed selection beat random?** |
| D5 | Structure of γ | Do context-specific perturbations cluster by pathway, lineage or target expression in a way that predicts *which* γ transfers? |
| D6 | Source selection / transport compatibility | Which source context is most relevant for a target, measured by something other than precision? |
| D7 | Foundation-model embeddings | Does embedding similarity predict response transfer better than pseudobulk basal similarity? |
| D8 | Failure taxonomy | Decompose zero-shot errors into abundance, low effect, noise, lineage, rewiring and nonlinearity |
| **D9** | **Shift-type ladder** | **How much cross-context discrepancy is replicate noise vs lab/screen vs donor genetics vs lineage?** |
| D10 | Prospective Arc validation | Predeclare per-target trust scores for D/E/F; test them on whatever Arc releases |
| D11 | Context-level design | Given a target lineage, which *source context* should be screened next to maximise transfer? |
| D12 | Component-attributed audit of SOTA | Does PIE/State/TxPert improvement come from β, template or γ, under reliability-normalised metrics? |
| D13 | Cross-perturbation denoising | Can low-rank structure in deep sources improve gene direction (C6 §I)? |

Aggressive filtering:

* **D3 rejected as a headline.** It is pre-empted by Molina & Zhang, Li et al. and Mechanisms Matter. With
  4–6 contexts an "impossibility" result is not credible: there are no near-twin basal pairs except within
  Feng, and those have near-identical responses. Keep it as a figure in D4's paper.
* **D5 rejected.** It is crowded and descriptive. Use it as interpretation inside D4/D9.
* **D6 rejected for now.** C3 and C6 showed that 3 atlases cannot identify donor weights. It needs ≥ 8
  contexts. It is subsumed by D11.
* **D7 rejected.** Crowded, and low information gain.
* **D8 rejected as standalone.** "Where Simple Baselines Fail" [M/S] covers it. Fold it into D9.
* **D10 opportunistic only.** Arc's leaderboard gives one aggregate number, not per-target outcomes. It is
  valuable only if the data are released.
* **D13 rejected.** Competition-engineering; the C6 gate failed on K562.
* **D2 demoted to a component.** Repo results exist, and ConfPert and PRESCRIBE are close. It would be a
  methods note, not a paper.
* **D1 kept, as a component of D4,** unless N4 shows it beats the exchangeable null.

---

## 7. Ranked future directions (Phase 4)

Scores are 1–5. For "data req.", "difficulty", "leakage risk" and "lit. overlap", **5 = favourable** (little
new data, easy, low risk, little overlap).

| Direction | importance | novelty | feasibility (current data) | data req. | difficulty | leakage risk | decisive | interpretability | experimental relevance | pub. fit | lit. overlap | **total /55** |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **D4 calibration budget** | 5 | 4 | 5 | 5 | 4 | 3 | 5 | 5 | 5 | 4 | 3 | **48** |
| **D9 shift-type ladder** | 5 | 4 | 3 | 3 | 3 | 4 | 3 | 4 | 4 | 4 | 4 | **41** |
| D12 component audit of SOTA | 4 | 3 | 3 | 4 | 2 | 4 | 4 | 5 | 3 | 3 | 4 | 39 |
| D1 prospective transferability | 4 | 2 | 5 | 5 | 5 | 4 | 3 | 4 | 4 | 3 | 2 | 41 |
| D11 context-level design | 4 | 4 | 2 | 2 | 3 | 3 | 2 | 4 | 5 | 3 | 4 | 36 |
| D2 selective prediction | 3 | 2 | 5 | 5 | 5 | 4 | 4 | 4 | 3 | 2 | 2 | 39 |
| D10 Arc prospective | 3 | 3 | 1 | 1 | 4 | 5 | 2 | 3 | 3 | 3 | 4 | 32 |
| D3 information limits | 3 | 1 | 3 | 4 | 4 | 4 | 2 | 3 | 2 | 2 | 1 | 29 |
| D5 γ structure | 3 | 1 | 4 | 4 | 4 | 4 | 2 | 4 | 3 | 2 | 1 | 32 |
| D6 source selection | 3 | 2 | 1 | 2 | 3 | 3 | 1 | 3 | 3 | 2 | 3 | 26 |
| D7 FM embeddings | 2 | 1 | 3 | 3 | 2 | 4 | 3 | 3 | 2 | 2 | 1 | 26 |

**Ranking:** D4 > D9 ≈ D1 (only as part of D4) > D12 > D2 (as part of D4) > D11 > the rest.

**Top three programs.** D4, D9 and D12.

* D1 and D2 survive only as **components** of D4. They are its k = 0 endpoint and its abstention rule.
* D12 is a cheap, sharp companion, but it depends on running others' code (PIE, State) faithfully.

---

## 8. Detailed top-3 research plans (Phase 5)

### 8.1 Program P1 (D4): the calibration budget of a new cellular context

**Research question.** How many measured perturbations in a new cell line are needed before cross-context
prediction of the remaining perturbations becomes reliable? Which response component does each measurement
buy, and can the anchor set be chosen better than at random using source data only?

**Novel claim, if it succeeds.** The existing claims are "zero-shot fails on γ" (Molina & Zhang) and "30 % or
~128 random anchors help" (State, Molina & Zhang, PerturbMap). We could say:

* *"In CRISPRi cell lines, the first k₁ ≈ 5–10 target measurements recover the context template and scale.
  γ needs k₂ measurements.*
* *A source-only selection rule reaches the random-anchor accuracy with x % fewer measurements.*
* *The zero-shot trust score predicts which perturbations will never need measuring."*

That is a quantitative, decomposition-explained budget law with a selection policy. No current paper states
one.

**Why it matters.**

* Screen designers deciding how big a pilot screen in a new line must be.
* Virtual-cell developers deciding whether zero-shot is the right target or few-shot is the honest product.
* Arc-style benchmarks, which could add a few-shot track with defined budgets.

**Primary hypothesis H1.**

* The k-anchor learning curve in template-removed (perturbed-reference) space has a knee at k ≤ 20
  perturbations.
* At k = 20, gain over the k = 0 scale-calibrated conserved predictor is ≥ 25 % of the gap to the
  k = all-but-test ceiling, in ≥ 4 of 6 held-out contexts.
* **Secondary hypothesis H2:** a source-only selection rule needs ≤ 0.7× the anchors of random to reach random's
  k = 20 performance, in ≥ 4/6 contexts.

**Null hypotheses.**

* **H1₀:** gains below k ≈ 30 % of the panel are confined to the template/scale. Template-removed γ recovery at
  k = 20 is within the bootstrap CI of k = 0.
* **H2₀:** no source-only selection rule beats random by more than its bootstrap CI over 200 random draws.

**Data.**

* *Sufficient now:* the frozen 4-context tensor (4 × 1,264 × 6,640).
* *Needed for a paper:*
  * The 6-context balanced design (+ X-Atlas HCT116, HEK293T; about 1,062 shared perturbations).
    Research-track license memo required.
  * Ideally a 7th context (KOLF; recount needed).
  * At least one independent dataset with distant contexts: Jiang 2025 six-line screens if counts allow, or
    the Arc 2026 data if released.
* Genome-wide K562 GWPS can widen the perturbation panel within K562 but is not a new context.

**Experimental design.**

* **Split:** LOCO over contexts (outer).
  * Within the held-out context, a random or selected anchor set A of size k is drawn from the panel.
  * Evaluation is on a **fixed test set T** (the same 30 % of perturbations for every k and every rule, drawn
    once by seed).
  * Anchors are drawn from the remaining 70 %.
  * 200 random draws per k.
* **Allowed:** source responses; target controls; target responses for anchors only.
* **Noise-coupling control (critical).** δ = pert mean − control mean, so control-mean error is shared by every
  perturbation of a context. An α estimated from anchors would absorb that shared error and "predict" it in T.
  Fix: split target control cells into halves.
  * Anchors and model fitting use control half 1.
  * Test truths use control half 2.
  * Test truths use perturbed-cell half B, and the reliability ceilings use the A/B halves.
* **Estimator ladder** (each a closed form or ridge, hyperparameters by inner leave-one-*source*-out that
  simulates anchors in source contexts):
  * **E0:** zero-shot source mean.
  * **E0s:** E0 × source-fitted scale (the current best).
  * **E1:** E0s + target template α̂ = mean over anchors of (δ_target − E0s).
  * **E2:** E1 + per-target scale refit on anchors.
  * **E3:** per-gene ridge regression of target on [source responses, 1], fitted on anchors (4 coefficients
    per gene, shrunk toward E2).
  * **E4:** reduced-rank regression of the anchor residuals on source-response principal components (rank
    ≤ 10), i.e. a low-rank γ map.
  * **E5 (reference only):** the target-only kNN/ridge within-context model on anchors (the Ahlmann-Eltze style
    ridge on PCA + STRING), which ignores sources. This shows when sources stop mattering.
* **Frozen features:** the gene set, perturbation panel, PCs, pathway sets and STRING version are fixed
  before any target is opened.
* **Selection rules** (source-only and anchor-free; frozen before the target is opened):
  * R0 random;
  * R1 k-center / diversity in source-response PC space;
  * R2 low source agreement (expected-γ-rich);
  * R3 high expected magnitude;
  * R4 pathway coverage (Hallmark);
  * R5 DPP / D-optimal design on source-response PCs.
  * The repo's finding that "low agreement = low signal" predicts R2 will *lose*. It is predeclared as the
    contrarian test.
* **Metrics** (on T, per context):
  * reliability-normalised perturbed-reference Pearson (template removed);
  * reliable residual fraction D (pooled);
  * energy explained;
  * component attribution: projection on held-out-only β and γ, used for **evaluation only**;
  * the Arc-like PDS on the 6-context design as a secondary metric.
* **Summary statistics:** k₅₀, the anchors needed to close 50 % of the E0s → ceiling gap, and the area under
  the learning curve.
* **Uncertainty:** 200 anchor draws × bootstrap over T perturbations; per-context reporting (no pooled-only
  numbers); Holm correction across the 5 selection rules.
* **Leakage controls:**
  * The noise-replace test extends to non-anchor target rows: they must not change predictions.
  * Selection rules may read only source arrays (a test replaces the target with noise and requires an identical
    anchor set).
  * The decomposition β/γ are evaluation-only.
* **Replication criterion:** H1 and H2 must hold in the 4-context design **and** the 6-context design **and**
  one independent dataset, all with the same frozen protocol.

**Critical baselines.**

* Zero-shot source mean (E0) and its scaled version (E0s).
* The "context mean-response" baseline (Arc's 0).
* Template-only calibration (E1). This is the baseline the paper must beat to claim anything about γ.
* Nearest-source.
* Target-only ridge (E5).
* Random anchors (R0).

**Killer experiment.** For each held-out context, plot the reliability-normalised template-removed accuracy vs
k (log axis) with a stacked attribution: template, scale, β-refinement, γ. Overlay random vs the best selection
rule with 200-draw bands. Show that the same curve shape replicates in an independent dataset. A reviewer
convinced by this figure believes that the budget law and the selection gain are real, and not artefacts of
the template.

**Failure condition** (predeclared; abandon or pivot):

* **(a)** Template-removed gain at k = 20 is not distinguishable from k = 0 in ≥ 3/4 (≥ 4/6) contexts. That
  means γ needs ≳ 30 % of the panel, which is already known. Pivot to P2.
* **(b)** No selection rule beats random after Holm correction in the 6-context design. Report the null:
  "random is optimal". This is publishable only as part of a negative/limits paper.
* **(c)** The 6-context design is unavailable (license) **and** no independent dataset is found. Then there is
  no paper-level replication; keep it as a workshop paper.

**Expected figures.**

1. **Fig 1.** Decomposition and noise budget across 6 contexts; zero-shot ceiling (replication, brief).
2. **Fig 2.** Zero-shot endpoint: trust-score risk-coverage per context; what never needs measuring.
3. **Fig 3.** The budget curve with component attribution (killer figure), per context.
4. **Fig 4.** Selection policy vs random; which source-only features make anchors informative; R2's predicted
   failure.
5. **Fig 5.** Replication on an independent dataset, and a practical table: "for a new line, measure k
   perturbations chosen by rule R to reach accuracy X".

**Positioning.** An ML-for-biology / experimental-design analysis paper with a methods component. Venues:

* **Stretch:** Nature Methods; Cell Systems.
* **Strong realistic:** Genome Biology; RECOMB or ISMB (Bioinformatics proceedings); a NeurIPS/ICML main track
  only if the selection policy is methodologically new.
* **Fallback:** Bioinformatics; PLoS Computational Biology; NAR Genomics & Bioinformatics; ICML GenBio / MLCB
  workshop.

### 8.2 Program P2 (D9): how much of "context-specificity" is the lab?

**Research question.** Is the cross-context discrepancy of CRISPRi responses dominated by biology (lineage,
genetic background) or by study effects (lab, library, chemistry, depth), once replicate noise is removed?

**Novel claim.** *"After reliability correction, a same-line/other-screen pair disagrees by X % of the
reproducible energy, compared with Y % for other-lineage pairs and Z % for other-donor pairs."*

* If X ≈ Y, a large share of benchmarked "context generalisation" is study generalisation. Every cross-context
  benchmark built from multi-lab data (Replogle–Nadig, X-Atlas, Arc vs public) is mis-attributing error.
* If X ≪ Y, γ is biological, and the case for context-conditioned models is stronger.
* Nobody has published this ladder for CRISPRi.

**Why it matters.**

* Benchmark designers.
* Anyone training on pooled public atlases.
* The Arc challenge, where public sources come from other labs and the target is Arc's own assay.

**Primary hypothesis H1.** The reliability-normalised discrepancy of same-line, different-screen pairs is
≥ 50 % of that of different-lineage, same-lab pairs.

**Null hypothesis.** Same-line discrepancy ≤ 20 % of the different-lineage discrepancy, i.e. study effects are
minor.

**Data.**

| rung | pair |
|---|---|
| replicate | split halves within every screen (have) |
| same line, same lab, other screen | K562 essential (Replogle) vs K562 GWPS (local) |
| same line, other lab | K562 Replogle vs VIPerturb-seq K562 (GREEN, not downloaded; median 47 cells/pert, so reliability is the risk); RPE1 Replogle vs Kaden (reliability dead, so only usable with block-reliability correction) |
| same lineage, other donor | Feng 19 iPSC lines (need counts for split halves); H1 vs KOLF (other lab too) |
| other lineage, same lab | K562 vs RPE1 (Replogle); HepG2 vs Jurkat (Nadig); HCT116 vs HEK293T (X-Atlas) |
| other lineage, other lab | all cross pairs |

**Design.**

* The quantity is reliability-normalised discrepancy: 1 − corr/√(ρ₁ρ₂) per shared perturbation, plus the
  reliable residual energy.
* A crossed random-effects model over pairs, with factors same-line, same-lab, same-lineage and same-donor,
  fitted to pair-level discrepancies.
* Match on perturbation set and depth: subsample cells to equal depth per pair.
* **Leakage:** none (descriptive). The **confounding** risk is high; mitigate by reporting each rung with ≥ 2
  pairs where possible.

**Baselines.** The replicate floor; depth-matched noise simulation; a permuted-perturbation null.

**Killer figure.** A "ladder" plot: reliability-normalised discrepancy (with CI) by rung, every pair shown as a
point.

**Failure condition.** Abandon if fewer than two usable (reliability ≥ 0.3) pairs exist at the
"same line, other lab" rung. Kaden and VIPerturb both risk this. The design then cannot separate lab from
line.

**Figures.**

1. The ladder.
2. Per-perturbation: which programs are lab-sensitive (e.g. stress/ribosomal) vs lineage-sensitive.
3. Consequence for benchmarks: rescoring published cross-context results after subtracting the study floor.
4. Feng donor-γ: the genetic-background contribution.
5. Recommendation for benchmark construction.

**Positioning.** An evaluation/benchmark analysis paper.

* **Stretch:** Nature Methods (Analysis); Genome Biology.
* **Strong realistic:** Cell Systems; Bioinformatics; RECOMB.
* **Fallback:** workshop; a short communication.

### 8.3 Program P3 (D12): where does SOTA improvement come from?

**Research question.** When PIE, State or TxPert beat baselines on unseen contexts in Replogle–Nadig, is the
gain in β (conserved), the template, or γ, after reliability normalisation?

**Novel claim.** *"Under component attribution and reliability normalisation, x of y published unseen-context
gains are explained by β/template recovery; none recovers γ beyond the trust-score-selected subset"*, or the
opposite. Molina & Zhang attributed MLP, STATE and MORPH; PIE and TxPert are new.

**Hypothesis.** PIE's held-out γ attribution is ≤ 0.1 (disattenuated) in ≥ 3/4 contexts.

**Null.** PIE's γ attribution ≥ 0.2 in ≥ 2/4 contexts. This would be a positive result for the field, and it
would also undercut P1's motivation.

**Data.** The 4-context tensor; PIE code (github.com/ArcInstitute/pie); State checkpoints. License checks are
needed (State weights are non-commercial).

**Design.** Run each model under our LOCO splits. Reuse their published checkpoints only if their training data
exclude the held-out context; otherwise retrain. Use our evaluation-only β/γ projections and reliability
ceilings.

**Baselines.** E0, E0s, source mean, Molina & Zhang MLP.

**Killer figure.** Per-model, per-context bars of disattenuated r with β and with γ, with the ceilings.

**Failure condition.** Models cannot be retrained faithfully under our splits within about 2 GPU-days, or
the authors' splits leak the held-out context. Then drop P3.

**Positioning.** A short analysis / Matters Arising-style piece or workshop paper. Better as **Fig. 2 of P1**
than as a standalone paper.

---

## 9. Recommended primary research direction

**P1: the calibration budget of a new cellular context.** Its components:

* zero-shot trust score (k = 0);
* decomposition-attributed few-shot learning curve;
* source-only anchor selection vs random;
* replication on ≥ 6 contexts plus one independent dataset.

**Why this, and not the proposed arc as written.** The arc proposed in the brief is decomposition → limits →
prospective transferability → selective prediction → minimal calibration. Challenged against the literature:

* **"Decomposition" and "limits"** are Molina & Zhang's. They are background, not contribution.
* **"Prospective transferability" and "selective prediction"** are things we already have. They are useful but
  thin as novelty: heterogeneity-as-confidence is natural, ConfPert and PRESCRIBE are adjacent, and our partial
  correlations are modest.
* **"Minimal calibration"** is the only link where:
  * our machinery (reliability-normalised, component-attributed, noise-coupling-aware evaluation) gives a
    real advantage;
  * the closest work (PerturbMap, State 30 %, Molina & Zhang 30 %, Pan 100) leaves the *small-k regime*, the
    *attribution* and the *selection policy* open;
  * the answer changes what experimentalists do.

So the arc is right in content but **wrong in emphasis**: invert it. Lead with the budget law, and use the
earlier links as the k = 0 boundary condition.

**Why it has the highest probability of a rigorous paper:**

* The first decisive test (N1) costs minutes on frozen data.
* The 6-context replication is on disk.
* Both outcomes are informative. A "γ needs ≳ 30 %" result plus "selection does not help" is still a clean
  limits result, with the trust score as the practical take-away.
* It reuses the predeclaration discipline unchanged.

## 10. Backup direction

**P2: the shift-type ladder (lab vs donor vs lineage).**

* **Why:** it attacks the biggest validity threat to *all* cross-context work, including ours. It has a real
  gap, and the repo already hit the phenomenon (Kaden CASE E).
* **Why it is the backup, not the primary:** its decisive rung, "same line, other lab, reliable", may not
  exist in public data. If VIPerturb K562 is too shallow and Kaden is dead, the ladder cannot separate lab from
  line. N5 tests that cheaply.

---

## 11. Next five experiments (Phase 7)

All five keep the frozen-predeclaration discipline. Each starts with a hashed predeclaration
(`reports/research_v3/nX_predeclaration.md`) and a freeze check. None changes modelling code in
`src/virtual_cell/competition_v2/` or the frozen research modules; new code goes in a new module with leakage
tests.

### N1: the k-anchor learning curve on the frozen 4-context tensor (highest information)

| field | specification |
|---|---|
| hypothesis | H1 of P1: template-removed gain at k = 20 ≥ 25 % of the E0s → ceiling gap, in ≥ 3/4 contexts |
| inputs | `data/processed/four_context_v1` tensor; per-context cell-level halves (perturbed and **control**) for noise decoupling; LOCO splits `loco_v1` |
| target | held-out-context responses on a fixed 30 % test set T (seeded) |
| split | outer LOCO; k ∈ {0,1,2,5,10,20,50,100,200,all}; 200 random anchor draws per k |
| baselines | E0, E0s, E1 (template-only), E5 (target-only ridge) |
| model | E2–E4 closed form / ridge; hyperparameters from inner leave-one-source-out |
| metric | disattenuated template-removed Pearson on T; pooled D; energy; β/γ attribution (evaluation only) |
| success | as in the hypothesis, with CI from draws × bootstrap; reported per context |
| compute | existing tensor (4 × 1,264 × 6,640 float); halves re-streamed once (~20 min, < 20 GB RAM, as in prior phases) |
| runtime | **CPU**, < 1 h total |

Outcomes:

* **Pass:** γ-space few-shot gain exists at small k → run N2 and N3.
* **Fail, with gains only in E1 (template):** a "template is cheap, γ is expensive" result. Record k₅₀ for the
  template. P1 shrinks to a limits note; run N5 for P2.
* **Noise-decoupling changes the answer** (gains vanish once controls are split): that is itself a methods
  warning worth reporting, since few-shot results may be inflated by shared control noise.

### N2: source-only anchor selection vs random (conditional on an N1 pass)

| field | specification |
|---|---|
| hypothesis | H2: the best predeclared rule reaches random's k = 20 accuracy with ≤ 0.7× the anchors in ≥ 3/4 contexts. Predeclared contrarian prediction: R2 (low agreement) is **worse** than random |
| inputs | as N1; source responses only for selection |
| target | as N1 |
| split | as N1; rules R0–R5 frozen and hashed before running |
| baseline | R0 random (200 draws) |
| model | the best N1 estimator, fixed |
| metric | k needed to match R0@20; area under the learning curve; Holm-corrected |
| success | as in the hypothesis |
| compute | tiny |
| runtime | CPU, < 30 min |

Outcomes:

* **Pass:** a selection-policy claim exists.
* **Fail:** "random is near-optimal" for new-context calibration. Report it; P1 becomes a budget-law paper
  without a policy.

### N3: 6-context replication with X-Atlas (research track)

| field | specification |
|---|---|
| prerequisite | a written license memo for non-commercial research use of X-Atlas (CC BY-NC-SA), recorded in `data_license_register.md`; code-level separation from competition candidates (already enforced) |
| hypothesis | N1 (and N2) conclusions hold with 6 contexts, ≥ 4/6 folds |
| inputs | 6 × ~1,062 × shared genes balanced tensor (K562, RPE1, HepG2, Jurkat, HCT116, HEK293T) built with the frozen v1 recipe; re-run the decomposition first (replication of #1 at 6 contexts) |
| target / split | LOCO over 6 contexts |
| baselines / metrics | identical to N1/N2 |
| success | same thresholds; same sign in ≥ 4/6 |
| compute | streaming X-Atlas lance files (~2 × multi-million cells); ~1–3 h, < 32 GB |
| runtime | **CPU** |

Outcomes:

* **Replicates:** the result is paper-level.
* **Fails in the X-Atlas folds only:** lab/assay shift (FiCS chemistry) is entangled. This feeds P2 directly.
* **Caveat:** X-Atlas differs in lab and chemistry, so the X-Atlas folds are not pure "new lineage".

### N4: is source agreement more than exchangeable heterogeneity? (decides the k = 0 component)

| field | specification |
|---|---|
| hypothesis | Observed agreement → transfer curves deviate from a fitted **exchangeable random-effects null** (β + iid γ + noise with fitted variance components). The deviation is concentrated in partner-dependent folds (K562/Jurkat) |
| inputs | frozen 4-context tensor; frozen confidence outputs |
| target | per-perturbation D in each LOCO fold |
| split | LOCO |
| baseline | simulated data from the fitted exchangeable model (1,000 replicates), with the observed depths and reliabilities |
| model | none; compare observed vs simulated Spearman(agreement, −D) and risk-coverage |
| metric | observed − null Spearman, per fold, with simulation CI |
| success | observed exceeds the null in ≥ 2/4 folds |
| compute | tiny |
| runtime | CPU, minutes |

Outcomes:

* **Agreement ≈ the null:** it is a heterogeneity statistic. Present it as a sanity baseline only; do not claim
  transferability prediction as novel.
* **Agreement > the null:** it captures non-exchangeable transport. Then it is a genuine, publishable
  component.

### N5: feasibility of the shift-type ladder (gate for P2)

| field | specification |
|---|---|
| hypothesis | At least two *reliable* (median split-half ρ ≥ 0.3 on shared perturbations) "same line, different screen" pairs exist in local or GREEN data |
| inputs | K562 essential (Replogle) vs K562 GWPS (local); VIPerturb-seq K562 (GREEN, ~needs download); Kaden RPE1 (local); Feng counts (MIT, needs download) |
| target | reliability-normalised discrepancy per pair on shared perturbations |
| split | none (descriptive) |
| baseline | the split-half replicate floor; the K562 vs RPE1 same-lab lineage pair |
| model | none |
| metric | 1 − r/√(ρ₁ρ₂); reliable residual energy fraction |
| success | ≥ 2 usable same-line pairs; then report X/Y (lab/lineage discrepancy ratio) as a pilot |
| compute | GWPS and essential already local; VIPerturb download size to be checked first |
| runtime | CPU, about 1 h |

Outcomes:

* **Usable pairs and X/Y ≥ 0.5:** a strong P2 signal, and a validity threat for P1 that must be reported.
* **Usable pairs and X/Y ≤ 0.2:** γ is biological; this strengthens P1.
* **No usable pairs:** P2 is infeasible with public data; drop it.

**Order.** Run N1 and N4 first: both are tiny and use only frozen data, and N4 could remove a claimed
component before anyone relies on it. Then run N2 if N1 passes, N5 in parallel, and N3 once the license memo
is written.

---

## 12. Potential title and abstract outline (recommended direction)

**Title options:**

* *"How many experiments does a new cell line need? Budget laws for calibrating perturbation-response
  predictions across cellular contexts"*
* *"From zero-shot to few-shot: what each measured perturbation buys in a new cellular context"*

**Abstract outline:**

1. **Problem.** Virtual-cell models are asked to predict CRISPRi responses in unseen cell lines, but the
   context-specific interaction is not predictable zero-shot (Molina & Zhang 2026). Practitioners therefore
   need a measurement budget, not just a model.
2. **Approach.** Reliability-corrected decomposition across N contexts; a zero-shot trust score; a few-shot
   learning curve with component attribution; source-only anchor selection, evaluated under LOCO with
   decoupled control noise.
3. **Result 1.** k₁ ≈ … measurements recover the template and scale; γ requires k₂ ≈ ….
4. **Result 2.** Selection rule R reaches random-anchor accuracy with … fewer measurements; uncertainty-style
   selection fails, as predicted, because low-agreement perturbations carry little signal.
5. **Result 3.** The zero-shot trust score identifies the … % of perturbations whose prediction does not
   improve with any anchor budget.
6. **Replication.** 6 contexts plus an independent dataset.
7. **Implication.** A design rule for pilot screens and a proposed few-shot track for virtual-cell benchmarks.

---

## 13. Evidence required before an honest submission

1. **≥ 6 contexts in the main design.** Include at least one context whose basal distance to all sources is
   Arc-like (r < 0.85), and at least two labs crossed with lineage. N5 or P2 must bound the lab effect.
2. **An independent dataset** (not Replogle/Nadig/X-Atlas) reproducing the curve shape and the selection
   result under the frozen protocol.
3. **Noise-decoupled evaluation** shown to matter or not matter (N1), with the control-split protocol in the
   main methods.
4. **Strong baselines run faithfully under our splits:**
   * Molina & Zhang MLP;
   * at least one of PIE, State or TxPert in its few-shot mode (or a documented inability, with reasons);
   * PerturbMap-style anchor transfer;
   * target-only ridge.
5. **N4 resolved,** so the k = 0 trust score is either claimed (beats the exchangeable null) or explicitly
   presented as a baseline.
6. **A non-essential perturbation panel tested.** Arc-like panels are not essential-gene libraries. Use K562
   GWPS / X-Atlas genome-wide perturbations in the 2–3 contexts where available.
7. **Full read of PIE, Molina & Zhang (latest version), PerturbMap and ConfPert**, with explicit differences
   stated in the related work.
8. **All predeclarations hashed before data access, and every negative arm reported.**

---

## 14. Reviewer #2 critique

**1. "This is Molina & Zhang plus PerturbMap with a smaller model."**

* *The criticism:* the decomposition and the γ-needs-target-data finding are published. Anchor-based
  cross-context transfer is published.
* *Neutralise with:*
  * Show something neither has: the small-k regime (1–50) with component attribution.
  * Show a selection policy that beats random under Holm correction.
  * Show the predicted failure of uncertainty-style selection, explained mechanistically.
  * Make a head-to-head with a PerturbMap-style estimator at matched k on the same data.

**2. "Four (or six) contexts from two or three labs cannot support a general law; context is confounded with
lab and assay."**

* *Neutralise with:*
  * The 6-context design plus an independent dataset.
  * P2's ladder, or N5, quantifying the lab floor.
  * Per-context (never only pooled) results.
  * A leave-one-*lab*-out analysis showing the curve shape is stable when the held-out context's lab is absent
    from the sources.

**3. "Gains at small k are just the template/scale, which is trivial and confounded by shared control-mean
noise."**

* *Neutralise with:*
  * The control-split noise decoupling.
  * Results reported in template-removed space.
  * An explicit E1 (template-only) baseline that every claimed γ gain must exceed.
  * A simulation with known template and γ showing the estimator does not create γ gains from template or noise.

**4. "Source agreement is just heterogeneity / perturbation strength; the trust score is not new and its
partial effect is weak."**

* *Neutralise with:*
  * N4's exchangeable-null comparison.
  * Partial correlations controlling for reliability and magnitude (already 0.14–0.55).
  * Comparison against Wang-style reliability/specificity triage and PRESCRIBE/ConfPert-style uncertainty
    under the same LOCO splits.
  * If agreement does not beat the null, demote it openly to a baseline.

**5. "Essential-gene screens with pseudobulk metrics do not represent real use; Arc-like panels and
single-cell distributions behave differently, and the absolute accuracies (cosine 0.03–0.09 on atlases) are
too low to matter."**

* *Neutralise with:*
  * A genome-wide (non-essential) perturbation subset in ≥ 2 contexts.
  * One Arc-like evaluation (vcc2026 metrics on the competition atlases, or the Arc 2026 data if released).
  * Absolute accuracies at each k alongside normalised ones.
  * A practical framing in "fraction of perturbations meeting a usable-accuracy threshold".
  * Single-cell distributional effects stated as out of scope, as a limitation.

---

## Sources consulted for this audit (primary, this session)

* Molina & Zhang 2026: https://www.biorxiv.org/content/10.64898/2026.07.24.740459v1
* PIE (Arc, 2026-10-05): https://www.biorxiv.org/content/10.64898/2026.10.02.756297v1 ; https://github.com/ArcInstitute/pie
* PerturbMap: https://arxiv.org/abs/2607.28090
* ConfPert: https://icml.cc/virtual/2026/72123
* State: https://www.biorxiv.org/content/10.1101/2025.06.26.661135v2.full
* TxPert: https://arxiv.org/abs/2505.14919
* Nadig et al. 2025: https://pmc.ncbi.nlm.nih.gov/articles/PMC11244993/
* Li et al. benchmark v3: https://www.biorxiv.org/content/10.1101/2024.12.23.630036v3.full
* Mechanisms Matter: https://icml.cc/virtual/2026/70734
* Jiang et al. 2025: https://pmc.ncbi.nlm.nih.gov/articles/PMC12083445/
* IterPert: https://www.biorxiv.org/content/10.1101/2023.12.12.571389v1.full
* Panagopoulos et al.: https://arxiv.org/html/2503.14571
* MapPFN: https://arxiv.org/abs/2601.21092
* Arc VCC 2026: https://arcinstitute.org/news/virtual-cell-challenge-2026 ; https://arcinstitute.org/news/behind-the-data-virtual-cell-challenge-2026
* VCC 2026 Cell commentary: https://www.sciencedirect.com/science/article/pii/S0092867426009311
* Other entries are as cited in `reports/literature_notes.md` (verified 2026-09-06).
