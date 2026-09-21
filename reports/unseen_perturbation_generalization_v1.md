# Unseen-perturbation generalization — v1

**Unseen `beta_p` is predictable on the internal benchmark and that result does
not survive a real context shift.**

On the four-context tensor, predicting a perturbation's conserved effect from
biological priors alone — with no measurement of that knockdown anywhere —
works: STRING-network k-neighbour transfer reaches `r = +0.51` and removes 25%
of the response energy, comfortably beating basal expression, pathways and
DepMap. Held out on both axes at once it degrades to `r = +0.24`, still real.

Then it is tested on `arch1`, a genuinely external context whose perturbations
appear in no source. There it collapses to **`r = +0.06`, leaving *more* error
than predicting zero** (1.24 against 1.00). The same collapse hits the context
main effect independently: even with a **perfect** scalar, 89% of `arch1`'s
main effect is unrecoverable from the sources, against 32–51% internally.

The internal benchmark was measuring transfer between four screens that share a
study design, a protocol and a gene panel. Arc's contexts are further apart
from each other (basal `r` 0.75–0.84) than ours are (0.89–0.93), so the
external number is the one that should be believed.

**A nonlinear model is not warranted.** The failure is not capacity-limited.

**No gamma model, no deep model, no Arc submission, no generative modelling.**
All prior conclusions stand unmodified; Arc bridge v1 is frozen.

Date: 2026-09-20
Reproduce:
`uv run python scripts/audit_public_datasets.py`
`uv run python scripts/run_unseen_perturbation.py`
`uv run python scripts/run_arch1_external.py`
Outputs: `outputs/unseen_perturbation_v1/`, `data/splits/arc_target_support_v1.csv`

---

## 1. Final seven-dataset Arc-target coverage

Recomputed **from the downloaded files**, not from the earlier remote audit.
`arch1`, `kaden25rpe1` and `wessels23` were acquired this phase; byte size and
upstream MD5 verified for all three.

| dataset | cells | genes | perts | median cells/pert | perts ≥400 cells | Arc genes | Arc targets |
|---|---|---|---|---|---|---|---|
| `replogle22k562` | 308,646 | 8,563 | 1,971 | 125 | 64 | 7,938 (42.8%) | **0** |
| `replogle22rpe1` | 240,774 | 8,749 | 2,016 | 82 | 37 | 8,259 (44.6%) | **0** |
| `nadig25hepg2` | 133,757 | 9,623 | 1,818 | 55 | 13 | 9,024 (48.7%) | **0** |
| `nadig25jurkat` | 258,202 | 8,881 | 2,137 | 89 | 42 | 8,284 (44.7%) | **0** |
| `arch1` | 221,273 | 18,020 | 150 | 1,045 | 126 | 18,017 (97.2%) | **13** |
| `kaden25rpe1` | 850,225 | 28,907 | 1,836 | 400 | 919 | 16,828 (90.8%) | **80** |
| `wessels23` | 28,490 | 16,775 | 157 | 151 | 11 | 15,439 (83.3%) | **0** |
| **union** | | | | | | **18,367 (99.1%)** | **86 (28.7%)** |

**Every figure matches the values reported in `arc_bridge_v1.md` exactly** —
per-dataset and union, genes and targets, including the 7 targets in ≥2
datasets. No material difference, so no stop condition was triggered.

Two facts only the local read could supply:

- **All seven files are log-normalized, not raw counts** (`arch1` values run
  0.10–4.65, `kaden25rpe1` 0.43–5.75). The four already in use are the same, so
  the frozen protocol applies to the new files unchanged. It also means the
  count-space requirement of an Arc submission is not satisfied by any public
  file and must come from the generators in `virtual_cell.arc.generate`.
- **`kaden25rpe1` sits exactly in Arc's cell-count regime**: median 400
  cells/perturbation, 919 perturbations with ≥400 cells.

## 2. Tier 0/1/2 counts

`data/splits/arc_target_support_v1.csv` — one row per Arc target, with contexts
observed, datasets observed, per-dataset cell counts, output-space presence and
tier. Classification uses **identifier presence only**; no response data is read,
so no outcome can influence it.

| tier | definition | count | share |
|---|---|---|---|
| **TIER 2** | perturbed in ≥2 public contexts | **7** | 2.3% |
| **TIER 1** | perturbed in exactly 1 public context | **79** | 26.3% |
| **TIER 0** | never perturbed in any public context | **214** | 71.3% |

All 7 Tier-2 targets are the `arch1` ∩ `kaden25rpe1` intersection: `HSBP1`,
`MTA1`, `RNF2`, `SMARCA5`, `STAT6`, `TARBP2`, `ZNF714`.

## 3. Feasible context-main-effect baseline

The VCC score anchors 0 at `m_c = mean_p delta[c,p] = mu + alpha_c`. For Arc's
A/B/C that is **hidden** — it summarises exactly the outcomes the challenge
withholds. The distinction is enforced in code: `oracle_main_effect` is an
evaluation target and is not reachable from `FEASIBLE_ESTIMATORS`, and a test
asserts that corrupting the target context's responses leaves every feasible
estimator bit-identical.

**Public leave-one-context-out** (unexplained fraction of `m_c`, lower better;
1.0 is what predicting the control scores):

| target context | E0 pooled | E1 basal-weighted | E2 basal-shrunk | direction `r` |
|---|---|---|---|---|
| K562 | 2.599 | 2.541 | **0.985** | 0.690 |
| RPE1 | 0.569 | **0.568** | 0.604 | 0.766 |
| HepG2 | 0.337 | **0.337** | 0.434 | 0.817 |
| Jurkat | 0.857 | 0.841 | **0.552** | 0.680 |

`E0`/`E1` recover the **direction** well (`r` 0.68–0.82) but the **magnitude**
badly — K562 overshoots 2.1×, RPE1 undershoots to 0.36×. The single scalar in
`E2` is what makes the estimator usable at all: it is the only one that beats
the zero response in all four folds.

**So no, target controls do not automatically deliver score 0.** Even `E2`
leaves 43–99% of the main effect unexplained.

### And on a genuinely external context it fails outright

On `arch1`:

| estimator | direction `r` | unexplained | norm ratio |
|---|---|---|---|
| E0 pooled | 0.329 | 18.63 | 4.54 |
| E1 basal-weighted | 0.331 | 18.11 | 4.48 |
| E2 basal-shrunk | 0.331 | **8.78** | 3.14 |

`arch1`'s responses are 3–5× smaller in norm than the four screens'
(median `‖delta_p‖` 1.29 against 3.50–6.30), so the pooled main effect
overshoots by 4.5×. But the problem is not only scale: **fitting the best
possible scalar with the oracle still leaves 89.0% unexplained**, because the
direction correlation is only 0.33. The same diagnostic run internally leaves
32–51%.

**The context main effect does not transfer to a distant context, and no
calibration fixes it.**

## 4. P1 — held-out perturbation, seen context

Perturbations are withheld from **every** context at once, so no measurement of
them exists in training. `mu` is recomputed from training perturbations only,
which is the easy leak to miss. Median over 5 folds.

| prior family | U1 nearest | U2 k-NN | U3 ridge | best unexplained |
|---|---|---|---|---|
| **string** | +0.297 | **+0.509** | +0.518 | **0.753** |
| combined | +0.393 | +0.499 | +0.512 | 0.761 |
| depmap | +0.128 | +0.366 | +0.407 | 0.889 |
| reactome | +0.200 | +0.282 | +0.373 | 0.939 |
| basal | +0.074 | +0.107 | +0.131 | 0.988 |
| hallmark | −0.087 | −0.016 | +0.218 | 0.987 |

(`beta_pearson`; U0 zero is `r` undefined, unexplained 1.000 by construction.)

**Unseen `beta_p` is predictable here.** STRING k-NN removes 24.7% of the
response energy with no measurement of the knockdown anywhere.

## 5. P2 — held-out perturbation AND held-out context

| prior family | U2 k-NN `r` | U2 unexplained | response `r` | response unexplained |
|---|---|---|---|---|
| combined | +0.230 | 0.973 | +0.249 | **1.003** |
| **string** | **+0.235** | **0.971** | +0.248 | 1.006 |
| depmap | +0.170 | 0.992 | +0.218 | 1.013 |
| reactome | +0.144 | 0.997 | +0.249 | 1.029 |
| basal | +0.050 | 1.012 | +0.215 | 1.031 |
| **U0 zero** | — | 1.000 | +0.225 | 1.023 |

### The ceilings that make these numbers mean something

| reference | `r` | unexplained |
|---|---|---|
| `oracle_beta` (perfect conserved effect) | +0.651 | **0.587** |
| `seen_perturbation_transfer_scaled` (Tier 2) | +0.324 | **0.903** |
| `seen_perturbation_transfer_raw` | +0.324 | 1.043 |
| best unseen prediction (STRING k-NN) | +0.235 | 0.971 |
| zero | — | 1.000 |

The residual of `oracle_beta` is pure `gamma` — 58.7% of the energy, which the
frozen conclusion says is not zero-shot predictable. So the whole available
budget for *any* beta-predicting model in P2 is 41.3%, of which scale-calibrated
Tier-2 transfer captures 9.7% and Tier-0 prediction captures **2.9%**.

Note `seen_perturbation_transfer_raw` scores 1.043 — *worse than zero*. The
frozen scale calibration is not a refinement, it is what makes conserved
transfer beat doing nothing at all.

## 6. Performance by prior family

**STRING network > DepMap > Reactome > basal > Hallmark**, consistently in both
regimes.

Answering question 4 directly: **yes, network, DepMap and pathway information
all beat simple target-gene basal expression.** Basal removes 1.2% of P1 energy;
STRING removes 24.7%, twenty times more. Hallmark is worse than useless —
74% of Arc targets are in none of its retained sets, and `U1` on Hallmark scores
`r = −0.087`.

Combining families adds nothing over STRING alone (0.761 against 0.753), so the
network is carrying essentially all the signal.

Coverage, which partly explains the ranking:

| family | features | leakage | coverage, training perts | coverage, Arc targets |
|---|---|---|---|---|
| basal | 7 | none | 0.848 | 0.677 |
| hallmark | 39 | none | 0.294 | 0.227 |
| reactome | 587 | none | 0.764 | 0.583 |
| string | 1,564 | none | 0.868 | 0.627 |
| depmap | 7 | **partial** | 0.942 | 1.000 |

DepMap is classified `partial` because gene effect **is itself a perturbation
outcome**. It is usable only because the readout (viability, not transcriptome),
the assay (knockout, not CRISPRi) and the cell panel all differ — and the panel
is the part that cannot be verified, since Arc's contexts are unidentified. The
mitigation is enforced in code: `depmap_features` **refuses** to emit per-line
profiles, so only summaries pooled over all DepMap lines are ever read and no
context-specific quantity can enter.

## 7. Performance versus support distance

**30.8% of Tier-0 Arc targets have zero STRING similarity to any training
perturbation** (Tier 1: 55.7%; Hallmark: 74.3%; Reactome: 39.3%). For those
genes the representation is not merely weak, it is empty.

Where the prior does cover them, Tier-0 Arc targets sit slightly *outside* the
training support: median max STRING similarity 0.403 against 0.485 for a
held-out training gene.

DepMap is the opposite failure — 0% zero-support, but median max similarity
0.9997 for every group. Seven pooled summary statistics cannot tell genes apart,
so its coverage is nominal.

## 8. Is unseen `beta_p` prediction scientifically viable?

**Internally yes; externally no, and the external answer is the one that
counts.**

`arch1` — 100 perturbations present in no source context, 6,209 shared genes,
a context from no source study. The estimator, its `k`, its ridge penalty and
the prior family were **fixed on the internal benchmark before `arch1` was
scored once**.

| estimator | beta `r` | beta unexplained | response `r` | response unexplained |
|---|---|---|---|---|
| U0 zero | — | 1.000 | +0.044 | 2.230 |
| U1 nearest | −0.017 | 5.918 | +0.006 | 5.866 |
| U2 k-NN | +0.057 | 1.704 | +0.033 | 1.302 |
| U3 ridge | +0.059 | **1.241** | +0.049 | 1.590 |

**No estimator beats predicting zero.** `r` falls from +0.235 internally to
+0.059, and every estimator leaves more error than doing nothing.

The contrast that isolates the cause — the 17 `arch1` perturbations that *are*
measured in the sources, scored on the same context and the same gene axis:

| | `r` | unexplained |
|---|---|---|
| Tier-2 scale-calibrated transfer (s = 0.517) | **+0.302** | **0.966** |
| Tier-0 prediction from STRING | +0.057 | 1.704 |

**Direct measurement of the perturbation still transfers to `arch1`. Predicting
it from priors does not.** The context shift is survivable when you have
measured the knockdown somewhere; it is not survivable when you are inferring
the knockdown from the gene's annotations.

## 9. Uncertainty estimator for Tier-0 genes

Source agreement is undefined for a Tier-0 gene — there are no source responses
to agree. Two candidates were tested.

**Support distance fails.** Spearman against per-perturbation accuracy: basal
**−0.254** (backwards), STRING −0.055 (nil), combined +0.017 (nil). Hallmark
scores highest (+0.279) and is the worst predictor, so the ranking is anti-
correlated with usefulness.

**Neighbour agreement works.** Rather than asking how close the gene is to the
training set, it asks whether the neighbours being averaged agree with *each
other* — the natural Tier-0 analogue of source agreement.

| prior family | support distance (P2) | **neighbour agreement (P2)** |
|---|---|---|
| depmap | +0.080 | **+0.502** |
| combined | +0.015 | **+0.466** |
| reactome | +0.189 | +0.426 |
| string | −0.073 | +0.260 |
| basal | −0.246 | +0.179 |
| hallmark | +0.251 | −0.039 |

It holds up externally: on `arch1`'s 100 unseen perturbations, Spearman
**+0.314**.

**Recommendation: neighbour agreement is the Tier-0 confidence statistic.** It
needs no target measurement, it is defined wherever the prior is, and it is the
one thing in this study that transferred to the external benchmark.

## 10. Simplest justified hybrid Arc strategy

Given section 8, the honest tiering is narrower than the phase plan anticipated:

| tier | n | prediction | confidence |
|---|---|---|---|
| **TIER 2** | 7 | scale-calibrated conserved transfer, `s` fitted inner-LOO on sources | raw source agreement |
| **TIER 1** | 79 | single-source response, scale-calibrated; **no prior blending** | raw source agreement (2 sources minimum is not met — treat as weak) |
| **TIER 0** | 214 | **the feasible context main effect alone** | neighbour agreement, reported but not acted on |

All tiers then pass through a count-space generator (`virtual_cell.arc.generate`).

**Tier 0 gets no prior-based `beta`.** The phase plan proposed one; the external
benchmark says it would add error rather than remove it. Adding it would be
optimising against the internal benchmark, which is precisely the mistake the
`arch1` result exposes.

The consequence for the Arc score is stark and should be stated plainly. Per
`arc_bridge_v1.md`, 0 is the mean-response baseline and a control-emitting
submission scores −0.31. This strategy puts:

- 214 targets (71%) at approximately the feasible main effect — and section 3
  shows that estimate is poor on a distant context, so **at or below 0**;
- 79 targets at weak single-source transfer;
- 7 targets at the validated Tier-2 predictor.

There is no honest path to a competitive Arc score from public data on this
evidence.

## 11. Is a nonlinear model warranted next?

**No.**

The evidence against is specific, not a general caution:

1. **The failure is not capacity.** Ridge already reaches the k-NN ceiling
   internally (0.826 against 0.761), and both collapse identically on `arch1`. A
   model with more capacity fits the internal benchmark better and would look
   *more* promising on exactly the benchmark that proved misleading.
2. **The budget is small even if it worked.** `oracle_beta` caps any
   beta-predicting model at 41.3% energy reduction in P2, and 58.7% is `gamma`,
   which is permanently closed.
3. **Two independent components failed together on `arch1`** — the context main
   effect (89% unexplained with a perfect scalar) and unseen `beta_p`. A better
   `beta` model does not address the first.
4. **30.8% of Tier-0 Arc targets have no STRING representation at all.** No
   architecture predicts from an empty feature vector.

What *would* change the answer, in order of value:

1. **Perturbation data on Arc-like targets in Arc-like contexts.** The binding
   constraint is 0/300 overlap with the four research contexts and 214/300 with
   no data anywhere. This is a data problem.
2. **A second external context** to confirm the `arch1` collapse is general and
   not specific to hESC biology or that study's CRISPRi efficiency. `arch1` is
   currently a single external observation, and that is the main caveat on this
   report's central claim.
3. **A representation with real coverage of Tier-0 genes** — sequence- or
   structure-derived embeddings would at least be defined for all 300, unlike
   STRING and Reactome.

## Validation

Two-axis nested protocol, with both axes held out in the hardest regime. Any
projection is fitted on training perturbations only. Leakage is not argued, it
is tested: for every estimator, replacing the outer context's responses with
noise, replacing the held-out perturbations' responses in every context with
noise, and doing both at once each leave predictions **bit-identical** — plus a
guard test confirming that corrupting the *readable* block does change them, so
the leakage tests cannot pass vacuously.

`uv run pytest` — 405 passed. `ruff check`, `ruff format --check`, `uv build`
clean. All 14 freeze manifests verified.
