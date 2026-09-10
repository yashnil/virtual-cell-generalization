# Literature notes: zero-shot perturbation prediction across unseen cellular contexts

Last updated: 2026-09-06 (revised the same day after corrections; see the
change note at the end). Compiled from web research (official Arc documentation,
bioRxiv/journal pages, official GitHub repositories). Nothing here is quoted at
length; every entry is a summary.

Evidence tags used throughout:

- **[P]** verified against a primary source (paper page, official repo, official docs).
- **[S]** secondary source only (search snippets, third-party summaries).
- **[?]** could not verify, or the source was inaccessible.
- **Established** = finding reported with data in a peer-reviewed or primary source.
- **Hypothesis** = position, single-dataset result, or our own inference.

Where a bioRxiv page was rate-limited (HTTP 429) during checking, the existence,
title, and date of the preprint were confirmed through search indexes and PMC
mirrors; numbers inside such entries come from the earlier agent reads and are
tagged accordingly.

---

## 1. The 2026 Arc Virtual Cell Challenge (the external benchmark)

Sources: Arc announcement https://arcinstitute.org/news/virtual-cell-challenge-2026 [P];
CLI wiki https://vcc-cli-wiki.virtualcellchallenge.org/ [P];
metrics brief https://github.com/ArcInstitute/cell-eval2/blob/main/docs/vcc2026_metrics/vcc2026-metrics-brief.md [P];
PyPI `vcc-cli` 0.2.0 (2026-09-01) and `cell-eval2` 0.16.0 (2026-08-20) [P];
Cell commentary DOI 10.1016/j.cell.2026.08.004 (abstract via Europe PMC) [P].

**Task** (Established, [P])

- Six undisclosed human cell lines "from different tissues of origin" (Arc
  announcement). A Cell commentary abstract uses the phrase "unseen cancer
  cell lines", but until official metadata or rules confirm it we describe
  them only as six cell lines from different tissues of origin. CRISPRi
  Perturb-seq with 10x Flex, sequenced on Ultima.
- Participants get, per context, non-targeting control cells (about 18,400 per
  context) and the list of target genes (about 300 per context, one construct
  per gene). The perturbed measurements are hidden.
- Validation contexts `A`, `B`, `C` (released 2026-08-20, live leaderboard).
  Final contexts `D`, `E`, `F` released 2026-10-22, with a different
  perturbation panel; validation and final scores are not comparable.
- Final submissions due 2026-11-05 23:59 UTC. Prizes: 100k / 50k / 25k USD
  plus compute credits.
- Any modelling strategy and any training data (public or private) are
  allowed. Team size, publication obligations, and the daily submission limit
  live on a client-rendered rules page that could not be fetched [?].

**Submission format** (Established, [P], CLI wiki v0.2.0)

- One `.vcc` file built by `vcc prep` from an `.h5ad`, covering all three
  contexts of the phase. Per-context files are rejected.
- `obs` columns: `target_gene` (gene symbol) and `context` (`A/B/C` or
  `D/E/F`). Optional cell-type column.
- Exactly 400 cells per perturbation per context (the wiki stresses "not a
  minimum"). Exactly 18,533 genes from `gene_names.csv`, order fixed by the
  tool. Raw, finite, non-negative, whole-number counts; normalised values are
  rejected. No control cells may be included.
- Complete submission = 300 x 400 x 3 = 360,000 cells.
- Discrepancy: the cell-eval2 metrics brief (2026-08-19) says cells per
  perturbation are "unconstrained"; the CLI (2026-09-01) enforces exactly 400.
  Treat the CLI as authoritative and re-check before submitting.

**Metrics** (Established, [P], metrics brief + wiki)

Six metrics, computed per context, averaged with equal weight:

1. Perturbation discrimination (`pds_cosine`): cosine distance between
   predicted and true pseudobulk deltas (log1p, target sum 5e4), rank-based;
   target genes removed; chance = 0.5.
2. Expression accuracy (`expr_mse_unbiased_capped_norm`): jackknife
   noise-corrected MSE relative to measured effect size, clamped to [0, 1].
3. DE direction fidelity (`de_wilcoxon_direction_fidelity_yield_raw`).
4. DE direction reach (`de_wilcoxon_direction_reach_raw`), purity floor 0.9.
5. DE significance overlap (`de_wilcoxon_sig_jaccard`).
6. DE log-fold-change accuracy (`de_wilcoxon_lfc_nmae`), floored at -6.

DE definition: two-sided Wilcoxon rank-sum on CPM, genes above 5 CPM in
control, Benjamini-Hochberg at 0.05 per perturbation, target gene excluded.

Normalisation: score = (u - b) / (r - b), where b is the **context
mean-perturbation-response baseline** (scores 0) and r is the **split-half
replicate** of the reference (scores 1). Scores are not clipped. Reference
values for A-C: PDS b = 0.500, r = 0.93-0.98; Jaccard b = 0.02-0.04,
r = 0.38-0.42.

Implication for us (Hypothesis): the zero anchor is the official
mean-response baseline, which is a single context-level response shared by
all perturbations. A perturbation-specific conserved-effect model
(beta_p transferred from other contexts, gamma = 0) is a different predictor
and can score above 0 on discrimination and DE metrics whenever beta_p
carries perturbation-specific signal. How much of the achievable score sits
in beta versus gamma is an empirical question for Phase 3, not something the
normalisation settles.

**2025 challenge, for contrast** (Established, [P])

- H1 hESC, 300 CRISPRi targets, about 300k cells, 150/50/100
  train/val/test split; three metrics (DES, PDS, MAE). Full dataset including
  held-out test perturbations is now public in
  `gs://arc-institute-virtual-cell-atlas/virtual-cell-challenge/2025/`.
- Wrap-up (https://arcinstitute.org/news/virtual-cell-challenge-2025-wrap-up):
  almost all models were worse than baseline on MAE; winners were hybrid
  deep + statistical; the Generalist Prize went to Altos Labs' flow-matching
  approach (preprint "PRiMeFlow", arXiv 2604.13986, link to Altos is [S]).

---

## 2. State (Arc Institute)

Adduri, Gautam, et al. "Predicting cellular responses to perturbation across
diverse contexts with State." bioRxiv 10.1101/2025.06.26.661135 (June 2025);
published in Cell, August 2026. Code https://github.com/ArcInstitute/state [P].

- **Architecture** (Established): State Transition (ST) is a bidirectional
  transformer with self-attention over *sets* of cells (set size 256), mapping
  a control cell set plus a perturbation embedding to a predicted perturbed
  set, trained with an MMD distribution loss. State Embedding (SE) is a cell
  embedding model (ESM2 gene embeddings, dual-axis loss) pretrained on about
  167M observational cells; SE-600M weights are public under a non-commercial
  licence.
- **Scale** (Established): ST trained on more than 100M perturbed cells across
  about 70 contexts (Tahoe-100M, Parse-PBMC, Replogle-Nadig).
- **Two distinct cross-context tasks** (Established, [P] v2 full text):
  - *Underrepresented context generalization task*: 30% of the held-out
    context's perturbations are included in training. This is where the
    headline gains are reported (+54% absolute perturbation discrimination on
    Tahoe-100M, +29% on Parse-PBMC; on Replogle-Nadig State matched the
    perturbation-mean baseline and beat other learned baselines).
  - *Zero-shot context generalization task*: no perturbations from the
    held-out context are used; ST with State Embedding (ST+SE) ranked
    perturbations by effect size better than the perturbation-mean baseline
    and than ST on raw expression, with gains of about 17% on larger datasets
    and several-fold on small genetic datasets where baselines are near zero.
  State is therefore evaluated in both regimes; the two must not be conflated
  when quoting numbers.
- **Stated limitation**: not evaluated on entirely held-out datasets where no
  context was seen in training.
- Reported baselines: linear models, CPA, scGPT, control/mean-type; metrics
  from cell-eval (PDS, DE overlap, log-FC correlation, effect size).

Implication: when citing State cross-context, name the task. Our
leave-one-context-out benchmark corresponds to State's zero-shot context
task, so only those numbers are comparable, and State should be re-run under
our split (ideally ST+SE) to be a fair baseline. Molina and Zhang (section 5)
report that State generalises poorly on both the conserved and interaction
components in held-out cell lines under a zero-shot split.

---

## 3. Benchmarks and the simple-baseline debate

**scPerturBench** — Wei et al. (Tongji), Nature Methods, online 2025-12-11,
issue Feb 2026. https://www.nature.com/articles/s41592-025-02980-0,
https://github.com/bm2-lab/scPerturBench [P]. This is the "27 methods on 29
datasets" benchmark named in plans.MD.

- Two scenarios: **cellular-context generalisation** (known perturbation,
  unseen context; i.i.d. and o.o.d. variants) and perturbation generalisation.
- Established: reported efficacy does not hold across diverse unseen contexts
  and unseen perturbations, foundation models included. Their bioLord-emCell
  variant, which conditions on cell-line embeddings, improves the o.o.d.
  context case, i.e. explicit context conditioning helps.

**Ahlmann-Eltze, Huber, Anders** — Nature Methods 22, Aug 2025.
https://www.nature.com/articles/s41592-025-02772-6 [P].

- Baselines: no-change, additive (for double perturbations), ridge on PCA
  embeddings. scGPT, scFoundation, GEARS, CPA, scBERT, Geneformer, UCE did not
  consistently beat them on Norman, Adamson, Replogle.
- Directly relevant (Established): in K562 to RPE1 and RPE1 to K562 transfer, a
  linear model pretrained on Replogle perturbation data beat everything;
  pretraining on observational atlases gave marginal benefit over random gene
  embeddings; accuracy tracked how similar a gene's behaviour was between the
  two lines.

**Wong, Hill, Moccia (Pfizer)** — Bioinformatics 41(6), May 2025 [P].
A "CRISPR-informed mean" (training mean with target gene set to zero for
CRISPRi) beat GEARS and scGPT; pretrained scGPT was no better than
randomly initialised scGPT.

**Rebuttals / reconciliation**

- Miller et al., "Deep Learning-Based Genetic Perturbation Models Do
  Outperform Uninformative Baselines on Well-Calibrated Metrics," bioRxiv
  2025.10.20.683304 [P]. Introduces an interpolated-duplicate positive control
  and a dynamic-range calibration measure across 14 datasets x 13 metrics;
  MSE and Pearson are poorly calibrated; under calibrated (weighted,
  rank-based) metrics deep models beat mean/control/linear. Hypothesis-level
  until independently replicated; directly contradicts the paper above.
- "Where Simple Baselines Fail," OpenReview May 2026 (authors [?]): about 40%
  of per-perturbation evaluations across 9 CRISPR datasets are near-solved by
  baselines, about 30% resistant; deep models beat baselines only on the
  resistant subset, recovering direction and identity of responding genes but
  not magnitudes [S].

**Cross-context specific**

- Qi and Chapfuwa, "Mechanisms Matter: Transportability of Cellular
  Perturbation Effects," bioRxiv 10.64898/2026.05.08.723625 (July 2026) [P].
  Causal-transportability framing; K562-anchored pairs from Replogle/Nadig plus
  Zhu 2025 T cells and a semi-synthetic simulator. Cross-context performance of
  GEARS, CPA, scVI, State "drops substantially, often to simple baseline
  levels"; State most stable; when the causal interaction structure differs
  between contexts, models collapse to baseline. Numbers are qualitative in the
  text we could read.
- Mao et al., "Benchmarking virtual cell models for in-the-wild perturbation
  response," arXiv 2604.27646 (April 2026) [P]: performance drops markedly
  under unseen context / unseen perturbation / cross-dataset; metric choice
  reorders rankings.

---

## 4. Context, not scale

Dibaeinia et al. (CZ Biohub Chicago, UChicago, Northwestern), "Virtual Cells
Need Context, Not Just Scale," bioRxiv 10.64898/2026.02.04.703804 (2026-02-09),
PMC12919078 [P].

- Position: the binding constraint is coverage over biological contexts, not
  model expressivity; framed as causal transportability (Pearl/Bareinboim)
  rather than covariate shift.
- Empirical (single dataset, Zhu et al. 2025 CD4+ T-cell Perturb-seq, held-out
  donor x timepoint contexts): perturbations seen in 8 training contexts reach
  mean DEG-F1 about 0.19 versus below 0.1 at 3 or fewer contexts; cell count
  correlates only r = 0.11 with DEG recovery; at matched cell counts, context
  diversity is significant (p = 0.002). Aggregate metrics correlate poorly with
  DEG-F1.
- Status: **Hypothesis** (position paper, one dataset). Useful framing and
  metric recommendations (report cross-context separately, add DEG recovery,
  difficulty tiers, make context an explicit conditioning variable).

---

## 5. Decomposing responses into conserved and context-specific parts

**Primary reference: Molina and Zhang (2026).** Alexis Molina and Xinyi Zhang
(AITHYRA, Research Institute for Biomedical Artificial Intelligence, Austrian
Academy of Sciences, Vienna), "Perturbation response decomposition enables
biologically aligned generalization to unseen perturbations and cellular
contexts," bioRxiv DOI 10.64898/2026.07.24.740459, posted 2026-07-27
(submitted 2026-07-24), CC BY 4.0. Code:
https://github.com/xinyizhanglab/perturbation-decomposition. [P: full text
and PDF read directly.] An earlier version of these notes wrongly stated that
this paper could not be found; that statement was an error in our search and
has been removed.

**Decomposition** (Established, [P]). Pseudobulk response
delta(c, p, g) = mean expression in perturbed cells minus non-targeting
controls, decomposed into four orthogonal, split-half noise-corrected
components:

    delta(c, p, g) = mu(g) + alpha_c(g) + beta_p(g) + gamma(c, p, g)

mu = global mean response; alpha_c = cell-line effect; mu + alpha_c together
are called the *cell-line template*; beta_p = perturbation effect conserved
across cell lines; gamma = cell-line-by-perturbation interaction. The authors
describe it as a diagnostic description of response variation, not a
mechanistic model.

**Data** (Established, [P]). Four CRISPRi Perturb-seq screens: K562 and RPE1
(Replogle 2022), HepG2 and Jurkat (Nadig). Per the Methods (full-text fetch):
K562 188,590 cells / 1,383 perturbations; RPE1 173,737 / 1,499; HepG2
96,616 / 1,340; Jurkat 184,470 / 1,537; perturbations with fewer than 30
cells excluded. Supplementary: genome-wide K562, Norman 2019 CRISPRa doubles,
Tahoe-100M drugs (5 lines, 268 compounds x 3 doses).

**Variance attribution** (Established, [P], Fig. 1c). Reproducible variance
of the full response across the four lines: template (mu + alpha) 27.8%,
conserved perturbation beta 29.4%, interaction gamma 23.5%, measurement noise
19.3%. After removing the template: beta 35.0%, gamma 27.7%, noise 37.4%.
Template share varies from 23% to 65% of response energy by cell line
(largest in RPE1 and HepG2). Essential-gene perturbations have about
two-fold larger template contribution.

**What control expression carries** (Established, [P]). Cell-line identity
axes from non-targeting cells explain only 0.9% to 7.4% of template variance.
The top control-expression PCs capture most of the template but the top 100
control PCs capture under 40% of the template-removed residual. Their
conclusion: generalising from unperturbed controls alone is insufficient;
biological priors about the perturbed gene are needed to identify
perturbation-specific axes.

**Gene priors and response alignment** (Established, [P]). Raw DepMap
coessentiality similarity does not correlate with template-removed response
similarity (near zero, Fig. 2b), yet Ridge or MLP models that *map* DepMap
PCA profiles (20 to 50 components suffice) onto pseudobulk responses match
MORPH within cell lines and beat scGPT. MORPH's own perturbation encoder
reduced recoverable signal relative to the raw DepMap readout.

**Cross-cell-line, zero-shot** (Established, [P], section 2.3 and 2.4).
Train Ridge/MLP on three source lines with DepMap features, evaluate on the
held-out fourth. Both predict positive template-removed signal; MLP beats
Ridge; more source lines help; this holds even when both the target line and
the perturbed gene are held out. Projecting the target's perturbation-aligned
bases onto the source's captures only 20 to 30% of target variance at rank 20.
Component attribution: MLP aligns with the conserved component beta at
r = 0.25 to 0.39 (Ridge 0.10 to 0.13). **No model recovered the interaction
gamma in held-out cell lines**, including STATE and MORPH, and adding source
or target control expression (or fine-tuning MORPH on target controls) did not
change this.

**The key negative finding** (Established for the information sources they
tested, [P]). Quoted from the abstract and section 2.4: the
perturbation-by-cell-line interaction component "cannot be predicted without
information from the target perturbation context." When 30% of the target
line's perturbations were added to training, STATE and MORPH did align
positively with gamma, and comparable interaction performance was obtained
with disjoint perturbation sets in source and target, so gamma is learnable
from measured target-context data but not from basal state, DepMap priors,
or source-context responses. "This boundary persisted across linear,
nonlinear and conditional generative models."

**Within-line unseen perturbations** (Established, [P]). Template component
predicted at r = 0.62 to 0.75; template-removed component r = 0.20 to 0.26;
Ridge assigns 55% (K562) to 94% (RPE1) of predicted variance to the template
direction. Aggregate metrics therefore mostly reward template recovery.

**Metric** (Established, [P]). Perturbed-reference Pearson (template-removed
Pearson): correlation after subtracting the training-set template from both
prediction and truth. Splits: within-line 5-fold over perturbations;
cross-line hold out the target line; double-unseen (line and gene);
leave-one-mechanism-of-action-out for drugs.

**Stated limitations and openings** (from the paper, [P]). Pseudobulk only,
essential-gene-enriched screens, four lines. Suggested sources of context
information not present in baseline expression: regulatory state, lineage,
perturbation efficiency, temporal dynamics, and richer single-cell
distributions. Their closing recommendation is to evaluate models by which
components they predict and what information makes recovery possible.

**Consequence for this project.** The project's central decomposition is
exactly theirs, and its central negative result is aimed at our original
hypothesis: basal state plus DepMap-style gene priors did not recover gamma
zero-shot. Direct prediction of gamma from basal state and gene priors is
therefore a high-risk hypothesis, not an assumed achievable outcome. See the
revised research questions at the end of these notes.

**Related work: COMPASS** (separate paper, kept as related work). Liang and
Singh (Duke), "COMPASS: Component-Wise Inference of Shared and Gene-Specific
Perturbation Response," bioRxiv 10.64898/2026.08.03.742643 (2026-08-06) [P:
existence and abstract verified; numbers from agent read of full text, S].

- Different decomposition (multiplicative-shared plus residual):
  z_cp = beta_cp * u_c + r_cp, with r_cp = g_p + eps_cp. u_c is the
  cell-line-wide shared response direction, beta_cp a scalar coefficient that
  is strongly conserved across lines, g_p a conserved gene-specific component,
  eps_cp the cell-line-specific deviation (their analogue of gamma). Not the
  same as the Molina/Zhang two-way ANOVA.
- Data: 2,270 CRISPRi perturbations shared across six lines (K562, RPE1,
  HepG2, Jurkat, HCT116, HEK293T), 2,000 HVGs.
- Transfer (leave-one-line-out on four lines) [S]: shared coefficient about
  0.66 Pearson; gene-specific component about 0.22; beta predictable from
  STRING embeddings at CV R^2 about 0.35.
- Discrimination (cosine PDS above chance) [S]: COMPASS about 0.23, TabICL
  0.20, GEARS/scGPT about 0.0, training mean 0.0 by construction.
- Reads consistently with Molina/Zhang: shared/conserved structure transfers,
  line-specific structure largely does not, and GEARS/scGPT capture mainly
  the shared part.

**Pan, Saunders, Replogle, Weissman, Zhuang** — "Global cell-state and
gene-program representations reveal conserved and context-specific
perturbation responses of cells," bioRxiv 10.64898/2026.05.16.725005
(2026-05-18) [P]. Mixture-of-experts over a 12M-cell manifold; K562, RPE1,
hESC genome-scale CRISPRi. Conserved: lysosome/autophagy, p53 activation.
Context-specific: mitochondrial perturbations trigger the integrated stress
response in K562 but not RPE1; the p53 axis is absent in p53-null K562.
Few-shot Pearson about 0.4-0.5 at 100 perturbations, plateauing; **zero-shot
underperforms the baseline**. No variance partition.

**Nadig et al.** — Nature Genetics 57 (2025-04-21),
https://www.nature.com/articles/s41588-025-02169-3 [P]. New HepG2 and Jurkat
CRISPRi screens (2,393 perturbations each; median 45 and 83 cells per target);
241 perturbations with cell-type-specific large effects (e.g. GATA1 in K562,
HMGCR in Jurkat). Established biological existence of the interaction term.

---

## 6. Existing context-aware or cross-context methods

| Method | Context conditioning | Unseen context from controls only? | Source |
|---|---|---|---|
| GEARS (Nat Biotech 2024) | none; per-dataset | no (README says not designed for cross-cell-type) | [P] github.com/snap-stanford/GEARS |
| CPA (MSB 2023) | learned discrete covariate embedding | no (unseen line has no embedding) | [P] PMC10258562 |
| scGPT (Nat Methods 2024) | none for perturbation task | no | [P] abstract; [S] details |
| biolord / bioLord-emCell | attribute disentanglement; emCell adds cell-line embeddings | partially, via external embeddings | [S] |
| scPRAM (Bioinformatics 2024) | VAE + OT; unseen cell type of same lineage for a seen perturbation | scGen-style, yes for seen perturbations | [S] |
| CellOT (Nat Methods 2023) | per-perturbation OT map | transports new control cells; no perturbation or context embedding | [S] |
| CellFlow (bioRxiv 2025) | condition encoder over perturbation and context covariates; flow matching | docs claim it; not the headline result | [P docs] [?] |
| PrePR-CT (bioRxiv 2024) | per-cell-type co-expression graphs from control cells; shared gene embeddings; GAT | yes, chemical perturbations, few per dataset | [P] |
| State (2025/2026) | set transformer over control cells + SE embeddings | yes in its zero-shot context task (ST+SE); larger gains in the 30% underrepresented-context task | [P] |
| X-Cell (Xaira, bioRxiv 2026-03) | diffusion LM with cross-attention to text/protein/network priors; 4.9B params | claims zero-shot to unseen primary cells | [P abstract] |
| MultiFlow (bioRxiv 2026-08) | control-derived state representation; RNA+ATAC flow matching | claims unseen contexts | [P abstract] |
| MORPH (bioRxiv 2025-07) | discrepancy-VAE + attention | claims new contexts; details [?] | [P abstract] |

**PrePR-CT in detail** (Alsulami et al., KAUST/Karolinska, bioRxiv
10.1101/2024.07.24.604816) [P]: builds a co-expression graph per cell type from
unperturbed cells only, uses shared learnable gene embeddings so overlapping
genes transfer similarity between graphs, GAT encoder, EMD loss, held-out cell
type chosen as most dissimilar. Reports R^2 on mean and standard deviation of
post-perturbation expression (e.g. Kang B cells 0.995 vs CPA 0.959). Caveats:
chemical/cytokine perturbations only, one to eleven perturbations per dataset,
no control-mean baseline reported, and R^2 on a mostly-unchanged transcriptome
is a generous metric. **This establishes that "context-specific graph prior
from control cells" is not novel**, as plans.MD already warns.

**Encodings benchmark** (Zinchenko et al., bioRxiv 10.64898/2025.12.15.694331)
[P abstract]: for DepMap-style response prediction, raw basal expression beat
foundation-model embeddings as a cell-line encoding. Supports using pseudobulk
basal state as a strong context baseline.

---

## 7. Evaluation pitfalls

- **Systema** — Vinas Torne et al., Nature Biotechnology 44(6), June 2026 [P].
  "Systematic variation" (the shared perturbed-vs-control shift) correlates
  0.91-0.95 with measured method performance, so Pearson-delta against a
  control reference largely scores the confound. Fix: perturbed-centroid
  reference, centroid accuracy, landscape reconstruction. Under the corrected
  reference, scores drop sharply for GEARS/scGPT/CPA.
- **Heidari et al.**, "Evaluating Single-Cell Perturbation Response Models Is
  Far from Straightforward," bioRxiv 10.64898/2026.02.14.705879 [P]:
  correlation and distributional metrics are driven by scale, sparsity, and
  dimensionality; Wasserstein fails under variance scaling; energy distance
  misses disrupted gene-gene dependencies.
- **Agarwal and Bisht**, "The Metric Picks the Winner," arXiv 2606.12639 [P]:
  ranking inversion between deep and linear models depending on metric, on a
  scaffold-split chemical screen.
- **PertEval-scFM** (Wu et al., ICML 2025) [P]: zero-shot foundation-model
  embeddings give limited gain over simple baselines under distribution shift.
- **Squair et al.**, Nat Commun 2021 [P]: cell-level Wilcoxon/t-tests treat
  cells as replicates and inflate false positives; pseudobulk with
  edgeR/DESeq2 matches bulk ground truth better. Note Arc's scorer nonetheless
  uses cell-level Wilcoxon, so we need both: Arc's definition for the challenge
  and a pseudobulk DE for the research track.

Established across these: control-referenced Pearson-delta and MSE are
confounded; DE-recovery and discrimination metrics are more informative; the
metric determines the conclusion.

---

## 8. Reliability and uncertainty

- **Wang, Kuipers, Hugi, Platt, Beerenwinkel**, "Reliable single-cell
  perturbations explain and improve model performance," bioRxiv
  10.64898/2026.08.11.744177 (Aug 2026) [P]. Split-half reliability rho and
  specificity phi per perturbation; of 7,170 perturbations across 29 datasets,
  65% unreliable, 11% shared, 24% specific. Training on reliable perturbations
  only matches full-data performance with 55% of them; model rankings change
  in 7 of 8 benchmark configurations when restricted to specific
  perturbations. Recommends ceiling-normalised metrics (Arc's split-half
  normalisation is exactly this). Data-side reliability, not per-prediction
  abstention.
- **PRESCRIBE** (Cheng et al., NeurIPS 2025, arXiv 2510.07964) [P]:
  evidential regression separating epistemic and aleatoric uncertainty for
  unseen genes; filtering low-confidence predictions improves accuracy.
- **Iversen, Renard, Baum** (bioRxiv 10.64898/2026.04.03.715851) [P]:
  risk-coverage template on drug response; ensembles flag distribution shift;
  64% MSE reduction at 10% coverage.
- **CILANTRO-SL** (bioRxiv 10.64898/2026.02.25.708096) [P]: conformal
  prediction for synthetic lethality with finite-sample guarantees.
- Gap: no paper found that does selective prediction / risk-coverage for
  transcriptome-wide responses in **unseen contexts**. The pieces exist
  separately.

---

## 9. Data resources

- **Replogle et al. 2022**, Cell 185 [P]: K562 genome-wide (9,866 genes,
  dCas9-BFP-KRAB), K562 essential (2,057 genes), RPE1 essential (ZIM3-KRAB);
  more than 2.5M cells; median above 100 cells per perturbation; dual-sgRNA.
  Hosted at gwps.wi.mit.edu and Figshare+ (10.25452/figshare.plus.20029387,
  page returned 403 to automated fetch [S]); also via pertpy and scPerturb.
- **Nadig et al. 2025** [P]: HepG2 and Jurkat essential-gene CRISPRi, GEO
  GSE264667.
- **X-Atlas/Orion** (HCT116, HEK293T) used by COMPASS [S]; check hosting.
- **scPerturb** (Peidli et al., Nat Methods 2024) [P]: 44 harmonised datasets
  [S count], unified `obs['perturbation']`, Zenodo 7041849; whether every
  matrix is raw UMI counts was not confirmed [?].
- **DepMap** [P forum]: latest release 26Q1 (2026-04-01); next is 26Q3.
  `CRISPRGeneEffect.csv` (Chronos); DepMap-generated data are CC BY 4.0.
  Coessentiality per Wainberg et al. 2021 (GLS with whitening) [S].
- **Arc 2025 dataset** [P]: fully public on GCS including held-out test.
- **Arc 2026 controls** [P]: via `vcc datasets download controls`.

---

## 10. Technical background (brief, [P] unless noted)

- **CRISPRi**: dCas9 fused to KRAB targets the promoter/TSS; KRAB recruits
  KAP1, leading to H3K9me3 and blocked initiation. No double-strand breaks, no
  frameshift mosaicism, reversible, tolerable for essential genes. Knockdown is
  partial (roughly 40-80% depending on domain and guide; ZIM3-KRAB strongest;
  Arc 2025 reports 83% of cells with more than 80% knockdown). Perturb-seq
  prefers CRISPRi because the partial loss of function is uniform and avoids
  DNA-damage transcriptional confounds.
- **Perturb-seq**: guide identity via expressed barcode or direct sgRNA
  capture; low MOI (about 0.3-0.6) for single guides; guide assignment by
  mixture models on guide UMIs [S]; typical 45-1,000 cells per perturbation.
- **Counts**: UMI counts are integers; library size varies about 10x within a
  cell type; overdispersed relative to Poisson and well fit by negative
  binomial; droplet data are not zero-inflated beyond NB sampling (Svensson
  2020). Standard normalisation: normalize_total then log1p. Our synthetic
  data are Poisson and therefore under-dispersed relative to real data.

---

## Established findings

1. Cross-context and unseen-perturbation generalisation degrades to near
   simple-baseline levels across independent benchmarks (scPerturBench,
   Mechanisms Matter, Mao et al., Systema, Pan et al. zero-shot).
2. Simple baselines (no-change, mean response, CRISPR-informed mean,
   PCA-ridge linear) are competitive with or beat deep models on
   control-referenced metrics; in K562/RPE1 transfer a linear model on
   perturbation data won (Ahlmann-Eltze 2025; Wong 2025).
3. Perturbation responses decompose into template (27.8%), conserved
   perturbation (29.4%), interaction (23.5%) and noise (19.3%) of reproducible
   variance across K562/RPE1/HepG2/Jurkat (Molina and Zhang 2026, [P]). The
   conserved component transfers to unseen lines (MLP r = 0.25 to 0.39); the
   interaction did not transfer from basal state, DepMap priors, or source
   responses, but became learnable once 30% of target-context perturbations
   were measured. COMPASS reports the same qualitative picture with a
   different decomposition (0.66 vs 0.22 [S]).
4. Context-specific responses are biologically real (Nadig 2025: 241
   perturbations with cell-type-specific large effects; Pan 2026: ISR and p53
   axis differences).
5. Metric choice changes conclusions; control-referenced Pearson-delta and MSE
   are confounded by systematic variation (Systema, Heidari, Agarwal-Bisht,
   Miller). Ceiling-normalised, discrimination, and DE-recovery metrics are
   preferred, and Arc 2026 uses them.
6. Most perturbations in public screens are unreliable at split-half level
   (65% of 7,170; Wang et al. 2026), so evaluation must be stratified by
   measurement reliability.
7. Arc 2026 is strictly zero-shot on context. State reports both an
   underrepresented-context task (30% of target perturbations in training)
   and a zero-shot context task; only the latter matches the Arc setting.
8. Control expression captures most of the cell-line template but under 40%
   of the template-removed residual (top 100 PCs), so basal state alone is
   insufficient for perturbation-specific prediction (Molina and Zhang).
9. Context-specific graph priors from control cells (PrePR-CT) and context
   covariate embeddings (CPA, biolord) already exist.
10. Response-aligned simple models (Ridge/MLP on DepMap PCA) match or beat
    MORPH, STATE and scGPT on unseen perturbations and unseen lines
    (Molina and Zhang); simple baselines remain the bar.

## Open questions

Revised research questions for this project (2026-09-06):

- **A. Identifiability.** How much context-specific perturbation response
  (gamma) is identifiable under truly zero-shot context shift? Molina and
  Zhang put 23.5% of reproducible variance in gamma but recovered none of it
  zero-shot with basal expression, DepMap, or source responses. What is the
  noise-corrected ceiling, per perturbation and per context, and does any
  of it correlate with anything observable in the unseen context?
- **B. Richer priors.** Can richer context-conditioned biological priors
  (pathway/TF activity from controls, lineage, regulatory-state proxies,
  coessentiality restricted to lineage-matched DepMap lines, single-cell
  heterogeneity of controls, perturbation-efficiency proxies) recover any
  predictable portion of gamma beyond the information sources already tested?
  This is a high-risk hypothesis; the prior evidence is negative.
- **C. Transferability triage.** When gamma cannot be predicted accurately,
  can we identify perturbation/context pairs for which the conserved effect
  beta is sufficiently transferable that a beta-only prediction is useful,
  and flag the rest for experimental measurement?

Further open questions:

1. Does the additive decomposition (Molina/Zhang) and the multiplicative
   decomposition (COMPASS) attribute the same perturbations to "large
   interaction"? Which is the better statistical model for count data?
2. Which basal-state features, if any, carry information about gamma? Molina
   and Zhang tested control expression PCs and fine-tuning on target
   controls; pathway-level, TF-activity, and single-cell-distribution
   features remain untested in their framework.
3. Do gene priors predict *which* perturbations have large interaction terms
   (perturbation-level transferability), as distinct from predicting the
   interaction vector itself? COMPASS's beta-from-STRING result hints yes.
4. Is context diversity really the scaling axis (only one dataset supports
   it)?
5. Do calibrated metrics rescue deep models (Miller et al.) or not
   (Ahlmann-Eltze)? Unresolved.
6. What are Arc's six cell lines, and how far are they from
   K562/RPE1/HepG2/Jurkat/HCT116/HEK293T in basal state?
7. How should 400-cell count distributions be generated so that Arc's
   Wilcoxon-based DE metrics are not dominated by sampling artefacts?
8. How does the Arc score decompose by component: what does a beta-only
   submission score on each of the six metrics, and how much headroom is
   left for gamma?

## Possible novelty opportunities

Provisional; not to be finalised until Phase 3 reproduces the decomposition
and quantifies the zero-shot ceiling on our own splits.

1. **Identifiability analysis of gamma under zero-shot context shift**
   (question A). Molina and Zhang report a single negative result; a
   per-perturbation, per-context, reliability-corrected analysis of how much
   of gamma is even in principle recoverable, and from what, would extend it.
2. **Testing richer context-conditioned priors against the negative result**
   (question B), framed as a high-risk hypothesis with a pre-registered
   stopping rule. A positive result would be significant precisely because
   the strongest current evidence is negative; a clean negative result across
   richer priors is itself informative and publishable as part of 1.
3. **Transferability prediction** (question C, new reliability direction):
   predict, for each perturbation/context pair, whether the response is
   likely conserved (small ||gamma||) or strongly context-dependent, from
   gene priors, source-context disagreement, and basal-state distance. If
   this works it makes beta-only predictions actionable even when gamma is
   unpredictable. Nearest prior art: COMPASS predicts the shared-response
   coefficient from STRING (R^2 about 0.35); Wang et al. classify
   perturbations by measurement reliability; PRESCRIBE handles unseen genes.
   No paper found that predicts context-dependence of a perturbation for an
   unseen context.
4. **A clean public zero-shot cross-context benchmark** with fixed splits,
   Arc-style ceiling-normalised metrics, reliability tiers, component-level
   attribution, and the full baseline ladder, plus prospective Arc results
   analysed by component.

## Ideas that are NOT novel

- The additive decomposition mu + alpha_c + beta_p + gamma_cp itself, and
  the finding that beta transfers while gamma does not (Molina and Zhang
  2026); related decompositions in COMPASS 2026, Pan 2026, Nadig 2025 TRADE.
- Mapping DepMap coessentiality onto responses with Ridge/MLP for unseen
  genes and unseen lines (Molina and Zhang 2026; MORPH).
- Conditioning on cell context via embeddings (CPA, biolord, bioLord-emCell,
  State, CellFlow).
- Context-specific graph priors built from control cells for unseen cell
  types (PrePR-CT 2024).
- Set-based transformers over control cells (State).
- Flow matching / OT for perturbation distributions (CellOT, CellFlow,
  PRiMeFlow, MultiFlow).
- Uncertainty for perturbation prediction in general (PRESCRIBE, CILANTRO-SL,
  Iversen et al.).
- Arguing that context diversity matters more than scale (Dibaeinia 2026).
- Showing that simple baselines are strong (Ahlmann-Eltze, Wong, VCC 2025).
- Using STRING/GO/coessentiality gene embeddings as perturbation
  representations (GEARS, COMPASS, Zinchenko et al.).

## Papers/models we must benchmark against

Baselines (mandatory, cheap):

- No-change (control mean of target context).
- Context mean-perturbation response (Arc's normalisation baseline).
- Conserved effect: mean delta for the same perturbation across source
  contexts; also the CRISPR-informed variant (target gene zeroed).
- Nearest-context transfer by basal similarity.
- PCA-ridge linear model (Ahlmann-Eltze) with context and perturbation terms.
- Additive context + perturbation model, then the low-rank interaction model.

Published methods (at least one or two, run under our split):

- Molina and Zhang's response-aligned Ridge and MLP on DepMap PCA (code
  public; this is now the primary reference method and the strongest known
  zero-shot cross-line baseline).
- COMPASS (related decomposition; code availability to check).
- State (few-shot in paper; re-run zero-shot with SE embeddings;
  non-commercial licence).
- GEARS and/or scGPT (widely used reference points; expect near-zero
  discrimination cross-context).
- bioLord-emCell (scPerturBench's context-conditioned method).
- Optionally TabICL as used in COMPASS, and PrePR-CT-style graph prior as an
  ablation rather than a competitor.

Benchmarks/metrics to reuse:

- Molina and Zhang's component attribution (project predictions onto beta
  and gamma) and perturbed-reference Pearson.
- cell-eval2 with `--preset vcc2026` for Arc metrics.
- scPerturBench context-generalisation splits for comparability.
- Wang et al. reliability labels (rho, phi) for stratification.
- Systema perturbed-centroid reference as a secondary check.

---

## Change note (2026-09-06, revision)

Corrections applied after review:

1. Molina and Zhang (bioRxiv 2026-07-24/27) exists and is the primary
   decomposition reference; the earlier claim that no such paper existed was
   wrong. COMPASS remains as related work.
2. Arc normalisation: 0 is the official mean-response baseline and 1 an
   experimental replicate anchor. A perturbation-specific conserved-effect
   model is not that baseline and can score above 0 with gamma = 0. Removed
   the statement that all score lives in the interaction.
3. Molina and Zhang's negative finding on gamma incorporated; research
   questions revised to A/B/C; direct gamma prediction reframed as high-risk.
4. Added the transferability-prediction direction.
5. State: distinguished the underrepresented-context (30%) task from the
   zero-shot context task.
6. Arc contexts described as six cell lines from different tissues of
   origin, not necessarily cancer lines.
