# Zero-Shot Virtual Cell Modeling Across Unseen Cellular Contexts

> **Status note (2026-09-25).** This is the original project specification
> (last revised 2026-09-09) and is kept as written. Several of its plans were
> overtaken by results: the Molina & Zhang reproduction gate in §62 is
> **closed** (their processed data is not public) and was replaced by the
> independent four-context gate, which passed; the interaction / pathway model
> of §16 and §36 is **terminated**; prior-based prediction of unseen
> perturbations failed externally; and no deep distributional generator (§19,
> §37) is justified. Current state: `README.md` and
> `reports/repository_state_notes.md`.

## A Research-Engineering Project in Scientific Machine Learning, Single-Cell Biology, and Human-Relevant Perturbation Modeling

---

# 1. Executive Summary

### Working project title

**Zero-Shot Virtual Cell Modeling Across Unseen Cellular Contexts**

### Potential paper title

**From Conserved Effects to Context-Specific Responses: Predicting Genetic Perturbations in Unseen Cellular Contexts**

Alternative:

**Learning Cell–Perturbation Interactions for Zero-Shot Virtual Cell Prediction**

### Core question

> **Can a machine-learning model predict how a genetic perturbation will behave differently in a cellular context it has never seen perturbed?**

More specifically:

> Can we explicitly model the **cell-context × perturbation interaction** that causes the same gene knockdown to produce different transcriptional effects in different cells?

This is a scientifically meaningful unsolved problem.

Recent work shows that perturbation responses contain several distinct components: a general response, a cell-context-specific component, a perturbation effect shared across cell types, and a **cell×perturbation interaction**. Importantly, the context-specific interaction represents substantial biological variation and remains poorly recovered by current models.

### Why this matters

Scientists can perform CRISPR perturbation experiments in some cellular systems, but many biologically important contexts are difficult, expensive, or impossible to perturb comprehensively.

Examples include:

* rare human cell populations;
* primary human cells that do not expand well outside the body;
* disease states difficult to reproduce in culture;
* tissue-specific cellular states dependent on their natural environment.

A useful virtual-cell model could help scientists use experiments that **have already been performed** to intelligently prioritize experiments in biological contexts where experimentation is scarce. Arc explicitly identifies this cross-context generalization problem as the motivation for its 2026 Virtual Cell Challenge.

### External validation opportunity

The **2026 Arc Institute Virtual Cell Challenge** provides an unusually valuable prospective evaluation.

Arc experimentally performed CRISPRi Perturb-seq across six new cellular contexts.

Participants receive only:

1. unperturbed control cells describing each context; and
2. the identities of genes that were knocked down.

The actual post-perturbation expression measurements are hidden.

Three contexts are used for validation; three completely different contexts are reserved for final testing.

Therefore:

> your model predicts biology first, and an experiment you never saw determines whether you were right.

That is much stronger scientific validation than creating your own random train/test split.

---

# 2. What You Are Actually Trying to Learn

This project should serve two purposes simultaneously.

## Scientific purpose

Investigate how genetic perturbation effects transfer across biological contexts.

## Personal training purpose

Develop the technical background expected of an **ML Research Engineer working in computational biology**.

By completing this seriously, you should understand:

### Biology

* DNA and genes;
* transcription;
* gene expression;
* gene regulation;
* cell state;
* CRISPR;
* CRISPR interference;
* single-cell RNA sequencing;
* Perturb-seq;
* differential expression;
* pathways and regulatory networks.

### Mathematics/statistics

* high-dimensional data;
* probability distributions;
* dimensionality reduction;
* regression;
* uncertainty;
* distribution shift;
* hypothesis testing;
* differential expression;
* bootstrap confidence intervals;
* experimental noise;
* calibration;
* possibly low-rank modeling.

### ML

* representation learning;
* Transformers/set models;
* conditional generative modeling;
* transfer learning;
* domain generalization;
* multimodal representation;
* OOD generalization;
* ablation studies;
* uncertainty estimation.

### Engineering

* PyTorch;
* AnnData/HDF5;
* Scanpy;
* large sparse matrices;
* GPU training;
* configuration systems;
* reproducible experiments;
* scientific data pipelines;
* model evaluation;
* unit/integration testing;
* experiment tracking.

This is an unusually good bridge from strong mathematics/CS into genuine computational biology.

---

# 3. Biological Background

## 3.1 Gene expression

Different cells contain essentially the same genome but behave differently because different genes are active at different levels.

For a cell, we can approximately represent its transcriptomic state as

$$
x =
[x_1,x_2,\ldots,x_G]
$$

where

$$
x_g
$$

represents the measured RNA abundance associated with gene \(g\).

With roughly tens of thousands of measured genes, one cell becomes an extremely high-dimensional vector.

---

# 4. What Is a Genetic Perturbation?

Suppose Gene \(p\) participates in some biological process.

Researchers can intentionally reduce its activity and observe what changes.

That intervention is a **perturbation**.

If we compare:

$$
X^{control}
$$

against:

$$
X^{perturbed}
$$

then approximately

$$
\Delta_{p}
=
E[X^{perturbed}]
-
E[X^{control}]
$$

describes the transcriptional effect of perturbing gene \(p\).

But the important complication is that:

$$
\Delta_{p}
$$

is usually **not universal**.

Instead,

$$
\Delta_{c,p}
$$

depends on the cellular context \(c\).

The same gene can participate in different pathways or regulatory environments in different cells.

That is the heart of this project.

---

# 5. CRISPRi

The Arc challenge uses **CRISPR interference**, or CRISPRi.

Instead of cutting DNA, CRISPRi is generally used to reduce transcription of a targeted gene.

Conceptually:

$$
\text{Gene activity}
\rightarrow
\text{CRISPRi}
\rightarrow
\text{lower gene expression}
$$

Scientists then ask:

> What else changed because that gene was suppressed?

This provides causal information that observational expression alone cannot provide.

Correlation tells us:

> Genes A and B are associated.

Perturbation can help ask:

> If I deliberately suppress A, does B change?

This distinction between **observation and intervention** is fundamental.

---

# 6. Perturb-seq

Perturb-seq combines:

**CRISPR perturbation**

with

**single-cell RNA sequencing.**

The experimental workflow is approximately:

$$
\text{cells}
$$

$$
\downarrow
$$

$$
\text{genetically perturb different genes}
$$

$$
\downarrow
$$

$$
\text{measure RNA in thousands/millions of individual cells}
$$

$$
\downarrow
$$

$$
\text{identify each cell's perturbation}
$$

$$
\downarrow
$$

$$
\text{learn what each intervention changed}
$$

This creates datasets containing enormous numbers of:

$$
(\text{context},\text{perturbation},\text{expression})
$$

observations.

The scPerturb resource provides harmonized public single-cell perturbation datasets specifically to support computational modeling and benchmarking.

---

# 7. What Is a Virtual Cell?

A useful virtual cell would approximately learn:

$$
f(x,p) \rightarrow y
$$

where:

* \(x\) = starting cellular state;
* \(p\) = intervention;
* \(y\) = resulting cellular state.

The ambitious long-term idea is:

> Give the computer a biological cell state and an intervention, and predict what the real biological experiment would produce.

Today's models are nowhere close to simulating an entire cell mechanistically.

The practical current problem is narrower:

> Can we accurately predict **transcriptomic responses to perturbation?**

That is already scientifically useful and extremely difficult.

Large models such as **State** now train on more than 100 million perturbed cells and explicitly target cross-context perturbation prediction, but generalizing reliably to genuinely unseen contexts remains a core challenge.

---

# 8. Why Cross-Context Generalization Is Hard

Imagine learning a perturbation in:

$$
\text{Cell A}.
$$

You then want to predict it in:

$$
\text{Cell B}.
$$

A naïve assumption is:

$$
\Delta_{A,p}
\approx
\Delta_{B,p}.
$$

Sometimes this works.

Sometimes it fails badly.

Why?

Because biological context changes:

* gene expression;
* transcription factors;
* signaling;
* regulatory networks;
* metabolic state;
* pathway activity;
* chromatin state;
* gene dependencies.

Recent work argues that virtual-cell modeling should not assume a single context-invariant perturbation function; context diversity itself is crucial for generalization.

A major Nature Methods benchmark evaluating **27 methods on 29 datasets** similarly concluded that current perturbation models struggle under diverse unseen-context conditions.

---

# 9. The Key Scientific Decomposition

Recent 2026 work gives us an extremely useful formulation.

Represent the perturbation response for:

* cell context \(c\),
* perturbation \(p\),
* gene \(g\)

as:

$$
\delta_{c,p,g}
=
\mu_g
+
\alpha_{c,g}
+
\beta_{p,g}
+
\gamma_{c,p,g}.
$$

Interpretation:

### \(\mu_g\): global response

Expression behavior broadly shared across experiments.

### \(\alpha_{c,g}\): cell-context response

General perturbation behavior characteristic of that cellular context.

### \(\beta_{p,g}\): conserved perturbation effect

The part of perturbing \(p\) that transfers across cellular contexts.

### \(\gamma_{c,p,g}\): cell × perturbation interaction

The biologically interesting part:

> **How does this specific cell context respond differently to this specific perturbation?**

Recent analysis estimates substantial reproducible variation in both the conserved perturbation component and the interaction component. Critically, current cross-cell-line predictors recover some of the conserved component but fail to recover the target-specific interaction well.

That is your opening.

---

# 10. Primary Research Question

Instead of asking:

> Can I build a bigger neural network?

ask:

> **Can information about the basal regulatory state of an unseen cell context and biological information about the perturbed gene predict the context-specific interaction component of the perturbation response?**

Formally:

$$
\hat{\gamma}_{c,p}
=
h(z_c,z_p)
$$

where:

$$
z_c
=
\text{representation of the target cellular context}
$$

and

$$
z_p
=
\text{representation of the perturbed gene}.
$$

Your prediction becomes:

$$
\hat{\delta}_{c,p}
=
\hat{\mu}
+
\hat{\alpha}_c
+
\hat{\beta}_p
+
\hat{\gamma}_{c,p}.
$$

This architecture encodes an explicit biological hypothesis.

---

# 11. Central Hypothesis

### H1 — Conserved perturbation effects exist

Some transcriptional consequences of perturbing a gene transfer across cellular contexts.

### H2 — Context-specific effects contain meaningful additional signal

The interaction

$$
\gamma_{c,p}
$$

contains biologically important information not explained by the conserved perturbation response.

### H3 — Basal context contains partial information about the interaction

The unperturbed state of the new context contains information about:

* active pathways;
* regulatory state;
* available transcription factors;
* gene dependencies.

This information should help predict how perturbation \(p\) changes in that context.

### H4 — Biological priors improve interaction prediction

Gene-function or gene-network information should help identify how perturbations operate.

Possible sources include:

* gene coessentiality;
* regulatory networks;
* pathway membership;
* protein/gene embeddings.

DepMap provides large public CRISPR dependency datasets that can serve as one source of gene-function information.

### H5 — Some predictions are fundamentally more transferable than others

The model should ideally identify conditions under which extrapolation is unreliable.

Uncertainty is therefore a useful **secondary research question**, but not the primary novelty claim.

---

# 12. Novelty Position

Do **not** claim that:

> response decomposition itself is new.

It isn't.

Do not claim:

> conditioning on cell context is new.

It isn't.

Cell-type-specific graph priors have already been used to generalize perturbation-response models to unseen cell types.

Do not build the paper solely around:

> uncertainty estimation.

Reliability-aware perturbation modeling is already beginning to appear in 2026 work.

### The provisional novelty opportunity is narrower:

> **Explicitly model the target-context-specific interaction component that current cross-context models fail to recover, using representations of basal cellular state × perturbation biology, and evaluate whether this improves zero-shot biological generalization.**

The second contribution could be:

> **Determine when the interaction component is predictable and when the model should defer to physical experimentation.**

This novelty hypothesis must still survive a serious literature review.

Do not freeze the paper claim until that review is complete.

---

# 13. Project Architecture

Think of the model as four components.

---

## Component A — Context encoder

Input:

$$
X_c^{control}
$$

the unperturbed cells describing context \(c\).

Produce:

$$
z_c.
$$

Potential representations:

### Baseline

Pseudobulk mean expression.

### Better

Principal components.

### Better still

Pathway-activity vectors.

### ML version

Pretrained single-cell embedding.

### Research version

Set encoder operating across a population of control cells.

State itself uses permutation-invariant Transformer-style modeling across sets of cells, providing a useful reference architecture.

---

# 14. Component B — Perturbation Encoder

Input:

gene \(p\).

Possible information:

### Identity embedding

Learn one embedding per perturbed gene.

Weak for unseen genes.

### Expression-derived representation

How the gene behaves across contexts.

### Functional network embedding

Use:

* coessentiality;
* protein interactions;
* pathway membership;
* gene ontology;
* regulatory networks.

A biologically useful perturbation representation should tell the model something about:

> what biological system this gene participates in.

---

# 15. Component C — Conserved Response Model

Learn:

$$
\hat{\beta}_p.
$$

This captures what tends to happen when gene \(p\) is perturbed across contexts.

Simple baseline:

$$
\hat{\beta}_p
=
\frac{1}{|C|}
\sum_c
\delta_{c,p}.
$$

Then compare against learned models.

A surprising amount of strong scientific ML comes from understanding when such simple baselines work.

---

# 16. Component D — Context-Specific Interaction Model

This is the main research component.

Learn:

$$
\hat{\gamma}_{c,p}
=
h(z_c,z_p).
$$

Possible first architecture:

$$
u_c = W_c z_c
$$

$$
u_p = W_p z_p
$$

then:

$$
u_{c,p}
=
u_c \odot u_p
$$

and:

$$
\hat{\gamma}_{c,p}
=
W_o u_{c,p}.
$$

This is essentially a low-rank multiplicative interaction.

Advantages:

* mathematically interpretable;
* computationally manageable;
* explicitly models interaction;
* easy to compare against additive alternatives.

Only after establishing this should you try:

* cross-attention;
* hypernetworks;
* conditional flows;
* larger Transformers.

Do not start with maximal architectural complexity.

---

# 17. Distributional Prediction

The Arc Challenge requires more than a pseudobulk vector.

Your model must submit predicted **single-cell raw count profiles**.

For the current validation bundle, the official tooling expects all three contexts, the official perturbation panel, **400 cells per perturbation**, all **18,533 genes**, and raw nonnegative integer counts. Current documentation describes a complete validation submission as 300 perturbations × 400 cells × 3 contexts = 360,000 predicted cells.

Therefore you eventually need:

$$
p(X^{perturbed}|c,p)
$$

rather than only:

$$
E[X^{perturbed}|c,p].
$$

---

# 18. First Distributional Baseline

Do not immediately train an enormous diffusion model.

Start with control-cell transport.

Take a real control cell:

$$
x_i^{control}.
$$

Predict a perturbation shift:

$$
\hat{\Delta}_{c,p}.
$$

Construct a transformed expression profile:

$$
\tilde{x}_i
=
T(x_i^{control},\hat{\Delta}_{c,p}).
$$

Then convert predicted expression proportions into count distributions while preserving realistic:

* library sizes;
* sparsity;
* cell-to-cell variability.

Possible count model:

$$
X_{ig}
\sim
\text{NegativeBinomial}(\mu_{ig},\theta_g).
$$

Later compare this with a generative distribution model.

---

# 19. Advanced Distributional Model

Once baselines work, investigate:

### Conditional flow matching

The 2025 Virtual Cell Challenge's Generalist Prize was won by an Altos Labs team using a flow-matching generative approach to model heterogeneous cellular responses.

### Conditional VAE

$$
z\sim q(z|X)
$$

conditioned on:

$$
(c,p).
$$

### Conditional Transformer / State-style model

Predict distributions across cell sets.

The important research discipline is:

> Don't use the most sophisticated model because it sounds impressive.

Use it only when a simpler model demonstrably fails to capture something biologically important.

---

# 20. Public Data Strategy

You should not train only for the Arc leaderboard.

Build a research dataset that supports honest cross-context evaluation.

## Dataset family 1 — Replogle / Nadig CRISPRi datasets

Recent cross-context work uses four CRISPRi Perturb-seq datasets from:

* K562;
* RPE1;
* HepG2;
* Jurkat.

These provide overlapping perturbations across different cell lines and are therefore highly valuable for your research question.

---

## Dataset family 2 — scPerturb

Use scPerturb to add additional harmonized perturbation datasets where appropriate.

Be careful:

different datasets can differ in:

* sequencing technology;
* guide design;
* preprocessing;
* perturbation strength;
* cell number;
* read depth.

Dataset heterogeneity is not merely inconvenience.

It is part of the scientific problem.

---

## Dataset family 3 — Arc 2025

Arc has released the complete 2025 Challenge dataset, including held-out perturbation responses, making it available for additional development and benchmarking.

Use this as another controlled benchmark.

---

## Dataset family 4 — DepMap

Use DepMap carefully for:

* gene dependencies;
* coessentiality;
* cell-context functional information.

This should be considered a **biological prior**, not direct ground truth for transcriptomic response.

---

# 21. Locked Research Evaluation

This is essential if you want a paper.

Create your own **leave-one-context-out** benchmark.

For example:

Train:

$$
\text{K562},\text{RPE1},\text{HepG2}
$$

Test:

$$
\text{Jurkat}.
$$

Repeat:

$$
C_1,C_2,C_3
\rightarrow C_4.
$$

Then rotate the target.

The model must never see perturbation responses in the held-out context.

It may only see its **control/unperturbed cells**.

That reproduces the Arc scientific problem.

---

# 22. Baselines You Must Beat

Your credibility depends heavily on your baselines.

## Baseline 0 — No perturbation

Predict target control state.

## Baseline 1 — Mean perturbation response

Predict generic perturbation behavior.

## Baseline 2 — Conserved perturbation effect

Average perturbation response across source contexts.

## Baseline 3 — Nearest context

Find the most biologically similar training context and transfer its response.

## Baseline 4 — Linear model

Context + perturbation additive model.

## Baseline 5 — Interaction linear/MLP model

Simple context×perturbation interaction.

## Baseline 6 — Existing specialized method

Include at least one respected modern perturbation model.

## Baseline 7 — State

If computationally feasible, State is an especially relevant baseline because it explicitly targets cross-context perturbation prediction and Arc provides public implementation resources.

If your sophisticated method doesn't beat the simple baselines:

report it.

That is science.

---

# 23. Metrics

Do not report only one metric.

The field has learned that different metrics measure fundamentally different biological properties. Recent work explicitly argues that evaluation choices can substantially change apparent conclusions about perturbation models.

The 2026 Arc scorer uses six metrics:

### Perturbation discrimination

Does the predicted response resemble the correct perturbation rather than other perturbations?

### Expression accuracy

How accurate is overall expression?

### Differential-expression log-fold-change accuracy

Did the model capture the magnitude of gene-level changes?

### Differential-expression direction fidelity

Did genes change in the correct direction?

### Differential-expression direction reach

How much of the true directional signal was recovered?

### Differential-expression significance overlap

Did predicted and observed sets of significant DE genes overlap?

The six normalized scores are averaged equally for the overall Arc score. Arc normalizes them so that 0 corresponds approximately to a simple context-mean reference and 1 to split-half experimental reproducibility, although most metrics are not mathematically capped at 1.

---

# 24. Research-Specific Metrics

For the paper, additionally report:

### Conserved response recovery

$$
corr(\hat{\beta}_p,\beta_p).
$$

### Interaction recovery

$$
corr(\hat{\gamma}_{c,p},\gamma_{c,p}).
$$

This may be the most important metric for your scientific thesis.

### Pathway recovery

Do predicted transcriptional changes recover the correct biological pathways?

### Top-DE recovery

For the strongest \(k\) differential genes:

$$
Precision@k
$$

$$
Recall@k.
$$

### Direction accuracy

For genuine DE genes:

$$
P(
\text{sign}(\hat{\Delta}_g)
=
\text{sign}(\Delta_g)
).
$$

These are biologically interpretable.

---

# 25. Secondary Reliability Study

Once the predictor exists, ask:

> When should we trust it?

Possible uncertainty features:

* ensemble variance;
* source-context disagreement;
* distance from target context to training contexts;
* perturbation support;
* network novelty;
* model residual variance.

Define:

$$
u_{c,p}
$$

as predicted unreliability.

Then evaluate whether:

$$
u_{c,p}
$$

correlates with realized prediction error.

Use:

### Risk–coverage curves

Keep only the most confident \(q\%\) of predictions.

If reliability is useful:

$$
Risk(20\%)
<
Risk(50\%)
<
Risk(100\%).
$$

But remember:

this is a **secondary contribution**, because reliability-focused work already exists in the area.

---

# 26. Biological Analysis

This is what turns the work from ML benchmarking into computational biology.

After obtaining predictions, identify:

### Highly transferable perturbations

Perturbations for which:

$$
\gamma_{c,p}\approx 0.
$$

Interpretation:

effects are relatively conserved.

### Strongly context-specific perturbations

Large:

$$
\|\gamma_{c,p}\|.
$$

Then ask:

* Which pathways do these genes belong to?
* Are they lineage-specific?
* Are they tied to regulatory programs active only in certain contexts?
* Do they have context-specific DepMap dependencies?

### Failure cases

Find perturbations where:

$$
\hat{\beta}
$$

is correct but:

$$
\hat{\gamma}
$$

is wrong.

These may be your scientifically most informative examples.

---

# 27. Ablation Studies

A serious paper needs ablations.

Remove one component at a time.

### A

No context information.

### B

Context mean expression only.

### C

Context embedding.

### D

Perturbation ID only.

### E

Gene network prior.

### F

Additive context + perturbation.

### G

Explicit interaction.

### H

Interaction + uncertainty.

Then answer:

> What actually caused the improvement?

---

# 28. Preventing Leakage and Self-Deception

Biological ML is extremely vulnerable to leakage.

Create hard rules.

### Rule 1

No perturbation responses from the target context during zero-shot evaluation.

### Rule 2

Freeze test contexts.

### Rule 3

Never tune hyperparameters using final test results.

### Rule 4

Record every leaderboard submission.

### Rule 5

Do not repeatedly optimize specifically against Arc validation until the model merely memorizes leaderboard quirks.

### Rule 6

Maintain a separate public benchmark untouched by Arc submissions.

### Rule 7

Document every preprocessing decision.

This dramatically increases scientific credibility.

---

# 29. Challenge Track vs. Research Track

Run the project as **two related tracks**.

## Track A — Arc Challenge

Objective:

maximize strong biological generalization across Arc's hidden data.

Deadline:

**November 5, 2026 at 11:59 PM UTC.**

Final test contexts are released **October 22, 2026**.

## Track B — Research Paper

Objective:

answer the broader scientific question.

This should exist independently of leaderboard ranking.

Even if you rank poorly on Arc, a strong scientific finding could still become valuable research.

Conversely:

ranking highly without understanding why your model works is less useful academically.

---

# 30. Engineering Repository

Build this like real research software.

```text
virtual-cell/
│
├── README.md
├── pyproject.toml
├── configs/
│
├── data/
│   ├── raw/
│   ├── processed/
│   └── splits/
│
├── src/
│   ├── data/
│   ├── preprocessing/
│   ├── biology/
│   ├── embeddings/
│   ├── models/
│   │   ├── baselines/
│   │   ├── context/
│   │   ├── perturbation/
│   │   ├── interaction/
│   │   └── generative/
│   ├── uncertainty/
│   ├── evaluation/
│   └── visualization/
│
├── experiments/
├── notebooks/
├── tests/
├── scripts/
│
├── paper/
│   ├── main.tex
│   ├── figures/
│   └── tables/
│
└── reports/
```

Notebooks should be used for exploration.

Core results should run from scripts.

---

# 31. Technology Stack

### Core

Python

PyTorch

NumPy

SciPy

pandas

### Single-cell

Scanpy

AnnData

scvi-tools where useful

### ML

PyTorch Lightning or clean custom loops

scikit-learn for baselines

### Experiment management

Weights & Biases, MLflow, or a simple reproducible local equivalent

### Configuration

Hydra or structured YAML

### Challenge

Arc's `vcc` command-line tooling

Arc's current CLI supports dataset download, validation, packaging, submission, and score retrieval.

---

# 32. September 1–November 5 Execution Plan

## Phase 1 — September 1–7

### Learn the biology and reproduce the task

Understand:

* CRISPRi;
* Perturb-seq;
* scRNA-seq;
* counts;
* normalization;
* pseudobulk;
* differential expression.

Read:

1. Arc 2026 Challenge description.
2. State.
3. Nature Methods perturbation benchmark.
4. 2026 response-decomposition paper.
5. Virtual Cells Need Context, Not Just Scale.

Deliverable:

**3–5 page literature memo.**

At the end write:

> What exactly appears unsolved?

Do not begin architecture work until you can answer that in your own words.

---

# 33. Phase 2 — September 8–14

### Data and baseline infrastructure

Download and preprocess:

* Arc validation controls;
* Replogle/Nadig;
* selected scPerturb data;
* Arc 2025;
* optional DepMap priors.

Build:

* AnnData loaders;
* gene intersection;
* pseudobulk pipeline;
* context splits;
* differential-expression analysis.

Implement:

* control baseline;
* mean perturbation baseline;
* conserved-effect baseline;
* nearest-context baseline;
* linear baseline.

Deliverable:

**Baseline benchmark table.**

No deep learning yet.

---

# 34. Phase 3 — September 15–21

### Establish the scientific gap

Implement response decomposition:

$$
\delta
=
\mu+\alpha+\beta+\gamma.
$$

Reproduce the qualitative finding:

> conserved effects transfer better than context-specific interactions.

Analyze:

* variance explained by each component;
* which perturbations have largest interaction magnitude;
* which contexts are hardest.

Deliverable:

**Figure 1 candidate for paper: anatomy of cross-context perturbation response.**

If you cannot reproduce the phenomenon, rethink the project.

---

# 35. Phase 4 — September 22–28

### Build context and perturbation representations

Context:

* pseudobulk expression;
* PCA;
* pathway scores;
* pretrained cell embeddings.

Perturbation:

* learned ID;
* gene expression features;
* coessentiality;
* biological pathway embeddings.

Perform representation ablations.

Deliverable:

**representation benchmark.**

---

# 36. Phase 5 — September 29–October 5

### Interaction model

Implement:

$$
\hat{\gamma}_{c,p}
=
h(z_c,z_p).
$$

Start simple:

* bilinear;
* low-rank multiplicative;
* MLP.

Compare against:

$$
\hat{\gamma}=0.
$$

Primary research metric:

$$
corr(\hat{\gamma},\gamma).
$$

Deliverable:

**first evidence supporting or rejecting the central hypothesis.**

---

# 37. Phase 6 — October 6–12

### Distributional generation

Move from pseudobulk response to single-cell distributions.

Baseline:

control-cell transport + count resampling.

Advanced:

conditional generative model if necessary.

Evaluate:

* mean expression;
* variance;
* sparsity;
* library size;
* DE behavior;
* Arc-style metrics.

Deliverable:

**valid end-to-end `.vcc` submission pipeline.**

---

# 38. Phase 7 — October 13–21

### Lock the research model

Run:

* complete ablations;
* bootstrap confidence intervals;
* multiple seeds;
* all public test contexts;
* biological interpretation;
* error analysis.

Add uncertainty if it produces real information.

Do not endlessly change the architecture because of Arc leaderboard movements.

By October 21 produce:

**Model Freeze v1.0**

with:

* commit hash;
* dataset versions;
* configuration;
* checkpoint;
* predictions;
* known limitations.

---

# 39. Phase 8 — October 22–November 5

### Prospective Arc final evaluation

Arc releases the three final cellular contexts on **October 22**.

They are labeled separately from the validation contexts, and the final evaluation uses a different hidden test panel; current tooling emphasizes that validation and final scores should not be directly compared.

Run:

$$
\text{control cells}
\rightarrow
z_c
\rightarrow
\hat{\delta}_{c,p}
\rightarrow
\text{single-cell predictions}.
$$

No hidden ground truth is available.

Perform only predeclared model-selection procedures.

Final submission:

**November 5.**

---

# 40. Post-Challenge — November/December

Once Arc releases results:

### Analyze

* hidden generalization;
* component-level success;
* failures;
* uncertainty;
* biological findings.

### Write

Research manuscript.

### Release

* repository;
* model;
* benchmark;
* reproduction script;
* figures;
* technical report.

Before releasing a preprint, verify the current Arc rules/FAQ concerning challenge publication and participant obligations. The Challenge FAQ explicitly includes guidance for researchers and students interested in publishing their results.

---

# 41. Core Deliverables

By the end, you should have **six artifacts**.

## Deliverable 1 — Research repository

Polished public GitHub repository.

Not:

> three notebooks.

A reproducible research codebase.

---

## Deliverable 2 — Scientific benchmark

Cross-context perturbation benchmark showing:

* baselines;
* your model;
* ablations;
* component recovery;
* confidence intervals.

---

## Deliverable 3 — Arc Challenge submission

Externally evaluated against experimental measurements you never observed.

---

## Deliverable 4 — Technical report

Approximately 8–12 pages.

Structure:

1. Motivation
2. Biology
3. Related work
4. Problem formulation
5. Method
6. Datasets
7. Evaluation
8. Results
9. Biological analysis
10. Human impact
11. Limitations
12. Future work

---

## Deliverable 5 — Research demo

Keep it small.

Example:

Select:

**cell context**

and:

**perturbed gene**

Then visualize:

* predicted DE genes;
* predicted expression changes;
* conserved response;
* context-specific response;
* model confidence;
* relevant pathways.

The demo should communicate the science.

Do not let it become the project.

---

## Deliverable 6 — Potential paper

Only pursue the paper if the project produces a defensible scientific contribution.

---

# 42. What Would Make This Paper-Worthy?

I would want at least three of these five.

### 1. Novel method

A new interaction-modeling approach that isn't merely architectural remixing.

### 2. Strong empirical result

Significant improvement in cross-context prediction across multiple datasets.

### 3. Interaction recovery

Evidence that your model predicts:

$$
\gamma_{c,p}
$$

better than current methods.

This is especially valuable because recent work identifies this component as a major failure point.

### 4. Prospective evidence

Strong performance on Arc's hidden experimental contexts.

### 5. Biological insight

Discovery of patterns explaining:

* which perturbations transfer;
* which are context-specific;
* why models fail.

The strongest paper would contain all five.

---

# 43. Potential Paper Contribution

A plausible contribution statement might eventually become:

> We show that zero-shot perturbation prediction can be improved by explicitly separating conserved perturbation effects from cell-context-specific interactions. We introduce a model that conditions interaction effects on basal cellular state and perturbation biology, improving recovery of context-dependent transcriptional responses across multiple held-out cell lines. Independent evaluation on previously unseen Arc Institute CRISPRi experiments demonstrates [RESULT]. Analysis identifies biological regimes in which perturbation effects transfer reliably and regimes where experimental measurement remains necessary.

Do **not** write `[RESULT]` until reality gives it to you.

---

# 44. Paper Figure Plan

### Figure 1 — Problem

Same perturbation.

Different cells.

Different responses.

Visualize:

$$
\beta_p
$$

versus:

$$
\gamma_{c,p}.
$$

### Figure 2 — Method

Control cell population

*

perturbation representation

↓

conserved effect

*

context-specific interaction

↓

predicted cell population.

### Figure 3 — Main benchmark

Methods versus unseen contexts.

### Figure 4 — Interaction recovery

How well each method recovers:

$$
\gamma.
$$

### Figure 5 — Biological interpretation

Transferable versus context-specific pathways.

### Figure 6 — Arc prospective result

Hidden experimental evaluation.

If Figure 6 is strong, it may become the most powerful figure in the paper.

---

# 45. Human Impact

Be precise.

Do not say:

> This model cures cancer.

Do not say:

> Virtual cells replace experiments.

The credible impact is:

### Today

Genetic perturbation experiments are expensive and context-specific.

Researchers cannot experimentally test:

$$
\text{every gene}
\times
\text{every cell state}
\times
\text{every disease}
\times
\text{every intervention}.
$$

### Computational opportunity

Use existing experiments to prioritize which experiments are likely to be most informative.

### Potential effect

$$
\text{huge perturbation space}
$$

↓

$$
\text{computational screening}
$$

↓

$$
\text{smaller experimental shortlist}
$$

↓

$$
\text{better targeted biological experiments}
$$

↓

$$
\text{better understanding of disease mechanisms}
$$

↓

potentially

$$
\text{better therapeutic targets}.
$$

Arc specifically frames cross-context virtual-cell prediction as useful when relevant human cellular states are difficult, expensive, sparse, or impossible to perturb directly.

---

# 46. A Second Human-Impact Dimension: Knowing When Not to Trust AI

This is important.

A virtual-cell prediction can influence what experiments scientists pursue.

Bad predictions therefore waste:

* money;
* biological material;
* scientist time;
* potentially rare samples.

A scientifically responsible model should ideally say:

> I have strong evidence that this perturbation transfers.

or:

> This context appears sufficiently different that you should measure it experimentally.

That is the human purpose of the secondary reliability work.

The objective is not:

> eliminate experiments.

It is:

> **use computation to determine where experiments are most valuable.**

---

# 47. Ethical / Scientific Limits

The report should explicitly state:

### Transcriptome is not the entire cell

RNA measurements do not fully capture:

* proteins;
* protein modifications;
* metabolites;
* spatial structure;
* morphology;
* temporal dynamics.

### Cell lines are imperfect

Cultured cell lines do not fully reproduce human tissue.

### Prediction is not therapeutic validation

Even perfect transcriptomic prediction does not establish that a treatment is:

* safe;
* effective;
* clinically useful.

### Biological uncertainty remains

Models should augment experiments rather than be marketed as replacements for physical evidence.

This honesty makes your work more credible.

---

# 48. What This Could Look Like on Your Résumé

Do not write a bullet before obtaining results.

But the structure could eventually be:

> **Developed an interaction-aware virtual-cell model for zero-shot CRISPRi perturbation prediction across unseen cellular contexts; integrated single-cell state representations and biological gene priors to model context-specific transcriptional responses across X public Perturb-seq datasets.**

Then, if supported:

> **Improved context-specific perturbation-response recovery by X% over strong linear/ML baselines and externally evaluated predictions against hidden 2026 Arc Institute Perturb-seq experiments.**

If challenge placement is strong:

> **Placed X/Y teams in the 2026 Arc Virtual Cell Challenge on experimentally held-out cellular contexts.**

If you write a preprint:

> **First-author preprint: “[TITLE],” evaluating context-specific generalization of virtual-cell models across X perturbation datasets and prospective Arc Institute experiments.**

Every number must be real.

---

# 49. Why This Is Strong for ML Research Engineer Internships

A recruiter or researcher could infer that you understand:

### ML

* high-dimensional modeling;
* Transformers;
* generative modeling;
* domain shift;
* representation learning.

### Math/statistics

* decomposition;
* regression;
* uncertainty;
* evaluation;
* experimental variability.

### Biology

* CRISPR;
* gene regulation;
* Perturb-seq;
* single-cell transcriptomics.

### Scientific reasoning

* baselines;
* leakage;
* controls;
* causal interventions;
* biological validation.

### Engineering

* large sparse data;
* reproducibility;
* GPU pipelines;
* production-quality research software.

That is almost exactly the profile you're trying to develop for **ML Research Engineering in biology**.

---

# 50. How to Use It to Work With Professors

Do not initially email:

> “I'm interested in AI and biology. Can I join your lab?”

Instead, once you have meaningful progress, send:

### Artifact 1

One-page research summary.

### Artifact 2

GitHub repository.

### Artifact 3

Two or three figures showing your result.

### Artifact 4

Short explanation:

> I'm investigating why virtual-cell models fail to generalize genetic perturbations across cellular contexts. Recent work suggests models recover conserved perturbation effects but struggle with cell×perturbation interactions. I've built [method], evaluated it on [datasets], and am testing it prospectively through Arc's 2026 hidden Perturb-seq benchmark. I'm particularly interested in understanding [specific scientific question].

Then ask:

> Would you be open to discussing whether this connects with work in your group?

That conversation is dramatically different from generic undergraduate outreach.

You're arriving with a **research question**.

---

# 51. What You Should Ask a Professor For

Not immediately:

> Can you publish my paper?

Instead:

> Is my biological framing correct?

> Is this actually novel?

> Which baseline am I missing?

> Does \(\gamma_{c,p}\) correspond to something biologically meaningful?

> What biological analysis would make this useful rather than simply computational?

> Is there a relevant experimental dataset I should use?

Eventually, if the collaboration becomes substantive, authorship and supervision can emerge naturally.

---

# 52. Internship Outreach

The same strategy applies to AI-biotech companies.

Instead of:

> I'm a Stanford student interested in biotech.

You can say:

> I'm working on zero-shot prediction of genetic perturbations across unseen cellular contexts. I've been studying when perturbation effects transfer across cells and explicitly modeling cell×perturbation interactions. My system is being evaluated on Arc's hidden 2026 Perturb-seq experiments. I'm interested in research engineering for biological foundation models / perturbation modeling.

Now you're talking their language.

---

# 53. What NOT to Spend Time On

Avoid:

### Fancy frontend

One week maximum at the end.

### Training a gigantic foundation model from scratch

You do not have to compete with laboratories training on hundreds of millions of cells.

### Endless leaderboard optimization

The research insight matters more.

### Ten datasets with inconsistent preprocessing before your first baseline

Start focused.

### Twenty architectures

Good science beats architecture tourism.

### Branding

The project does not need a clever startup name.

The problem is interesting enough.

---

# 54. Minimum Successful Outcome

Even if the research hypothesis fails:

You complete:

* serious literature review;
* Perturb-seq pipeline;
* rigorous baseline suite;
* multi-context benchmark;
* Arc submission;
* failure analysis;
* public research repository;
* technical report.

That is already valuable.

You will have actually learned computational biology.

---

# 55. Strong Outcome

Your interaction model improves:

* cross-context prediction;
* DE recovery;
* interaction recovery;

across multiple public cell contexts.

Arc hidden performance is competitive.

You produce:

* polished repository;
* technical report;
* research presentation.

Very strong résumé project.

---

# 56. Exceptional Outcome

You establish a reproducible scientific finding such as:

> basal regulatory state × perturbation functional representation predicts meaningful context-specific response components that existing methods systematically miss.

The finding survives:

* multiple datasets;
* strong baselines;
* ablations;
* hidden Arc evaluation.

Then:

**write the paper.**

Possible path:

1. faculty feedback;
2. manuscript;
3. bioRxiv;
4. computational-biology / ML-for-science venue;
5. internship and research outreach.

---

# 57. The Reading Sequence

Do not read 50 papers before coding.

Start with these conceptual anchors:

### 1. Arc 2026 Virtual Cell Challenge

Understand the exact practical task.

### 2. State

Understand modern large-scale virtual-cell modeling.

### 3. Nature Methods perturbation benchmark

Understand why evaluation/generalization are difficult.

### 4. Virtual Cells Need Context, Not Just Scale

Understand the context-generalization argument.

### 5. Perturbation Response Decomposition

This is the most important paper for your specific research question.

### 6. PrePR-CT

Understand how biological context-specific graph priors have already been used, so you don't accidentally claim that idea as novel.

### 7. Principled Evaluation of Perturbation Models

Understand how evaluation choices affect scientific conclusions.

Then expand the literature review based on your exact method.

---

# 58. Your Research Notebook

Maintain a single document with four sections.

## Things I know

Facts supported by literature/experiments.

## Things I suspect

Hypotheses.

## Things I don't understand

Questions.

## Things my experiments disproved

This fourth section is important.

Research isn't about making your original idea look correct.

It is about discovering what survives contact with evidence.

---

# 59. The Question You Should Be Able to Answer at the End

Not:

> What is a virtual cell?

You should be able to explain:

> **Why is cross-context perturbation prediction difficult?**

> **Which components of a perturbation response actually transfer?**

> **What information is needed to predict context-dependent effects?**

> **Which current methods succeed and fail?**

> **How should virtual-cell predictions be evaluated biologically?**

> **When should a scientist trust a computational prediction versus run the experiment?**

If you can answer those deeply after this project, you'll have moved from:

> student interested in biotech

to:

> someone beginning to think like a computational biology researcher.

---

# 60. Final Project Definition

### Domain

Computational biology / functional genomics.

### Technical field

Scientific machine learning.

### Biological data

Single-cell transcriptomics.

### Intervention

CRISPRi gene perturbation.

### ML problem

Zero-shot domain/context generalization.

### Scientific question

Can the context-specific component of a perturbation response be predicted from basal cellular state and perturbation biology?

### Engineering challenge

Generate realistic post-perturbation single-cell count distributions at scale.

### External benchmark

2026 Arc Institute Virtual Cell Challenge.

### Human impact

More intelligently prioritize biological experiments in difficult-to-measure human cellular contexts.

### Career signal

ML Research Engineering for Biology.

### Paper opportunity

**Real, but contingent on scientific novelty and results—not guaranteed simply by participating in the competition.**

---

# 61. What This Project Should Become

The ideal final package is:

> **Research paper / technical report**

*

> **reproducible ML repository**

*

> **externally evaluated Arc Challenge result**

*

> **biological analysis**

*

> **small scientific demo**

*

> **clear understanding of single-cell biology**

*

> **a credible story for why you want to pursue ML research engineering in biology.**

If you execute that well, this should not read on your résumé as:

> “Participated in a Kaggle-style competition.”

It should read as:

> **“Conducted independent scientific ML research on an unsolved problem in computational biology, developed a new modeling approach, evaluated it rigorously across public datasets and prospectively against hidden experimental measurements, and investigated when computational predictions can meaningfully reduce biological search.”**

That is the standard I would set for the project.

---

# 62. Framing Addendum (2026-09-06)

This addendum revises sections 9 through 12 after reading Molina and Zhang,
"Perturbation response decomposition enables biologically aligned
generalization to unseen perturbations and cellular contexts," bioRxiv
10.64898/2026.07.24.740459 (July 2026). Full notes and sources are in
`reports/literature_notes.md`.

### What section 9 refers to

The decomposition in section 9 is Molina and Zhang's. On K562, RPE1, HepG2
and Jurkat they attribute reproducible variance as: template (mu + alpha)
27.8%, conserved perturbation beta 29.4%, interaction gamma 23.5%, noise
19.3%. COMPASS (Liang and Singh, August 2026) is a separate related
decomposition.

### The negative result that changes section 10

Molina and Zhang report that gamma could not be predicted in held-out cell
lines without measured perturbations from that context, using basal
expression, DepMap coessentiality, source-context responses, or fine-tuning
on target controls, across linear, MLP, MORPH and STATE. Beta did transfer
(MLP r = 0.25 to 0.39). Direct prediction of gamma from basal state and gene
priors is therefore a **high-risk hypothesis**, not an assumed result.

### Revised research questions (replaces the single question in section 10)

A. How much context-specific perturbation response is identifiable under
truly zero-shot context shift?

B. Can richer context-conditioned biological priors recover any predictable
portion of gamma beyond existing approaches?

C. When gamma cannot be predicted accurately, can we identify
perturbation/context pairs for which the conserved effect beta is
sufficiently transferable?

### Added direction (extends section 25)

Predict whether a perturbation is likely to be conserved across contexts or
strongly context-dependent (transferability prediction). Useful even if exact
gamma prediction remains out of reach.

### Corrections to sections 17 and 23

Arc normalises each metric so that 0 is the official mean-response baseline
and 1 an experimental replicate anchor. A perturbation-specific
conserved-effect model is not that baseline and can score above 0 even with
gamma = 0. The six contexts are described as cell lines from different
tissues of origin; do not assume they are all cancer lines until official
metadata says so. State reports both an underrepresented-context task (30%
of target perturbations in training) and a separate zero-shot context task;
cite the right one.

### Revised Phase 3 gate (validation gate, set 2026-09-06)

Reproduction of Molina & Zhang is a hard gate. No novel model is implemented
until every item below is reproduced on our own pipeline:

1. the decomposition delta = mu + alpha_c + beta_p + gamma_cp and its
   variance attribution;
2. the split-half reliability (noise) correction;
3. leave-one-context-out evaluation;
4. the beta-only transfer baseline;
5. response-aligned Ridge and MLP baselines.

Each must be run with documented dataset versions (source, accession or
download date, checksum) and frozen splits committed as files. Report the
zero-shot ceiling for gamma alongside these numbers.

If our numbers materially disagree with the paper, stop and investigate the
discrepancy (preprocessing, gene set, split definition, metric normalisation)
rather than proceeding to modelling.
