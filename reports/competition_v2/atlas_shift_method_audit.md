# AtlasShift — method audit (competition track)

**Scope.** A line-by-line audit of the public AtlasShift solution, read from source,
not from its README or from the phase brief. Competition track only. Nothing here
changes a frozen research conclusion (see §9).

| field | value |
|---|---|
| repository | https://github.com/kaipengm2/Virtual-Cell-Challenge-2026 |
| audited commit | `d24ce4fdae0cd7cba3ba29546cd8737094eae98a` (2026-09-23) |
| license | MIT, "Copyright (c) 2026 Virtual Cell Challenge 2026 solution contributors" |
| vendored copy | `third_party/atlasshift/` (byte-identical; provenance `third_party/README.md`) |
| files audited | `README.md`, `sources.json`, `requirements.txt`, `prepare.py` (598 lines), `model.py` (444), `predict.py` (228), `pack.py` (132), `compact.py` (249) |
| reported score | 0.1545618019, public leaderboard #82 on 2026-09-10 (README; **a claim, not reproduced by us** — a C0 bundle was built but not submitted, see the source audit) |
| runtime pins | Python 3.13, numpy 2.5.2, scipy 1.18.1, pandas 3.0.5, anndata 0.13.3.post0, h5py 3.16.0, numba 0.67.0; `slafdb` for X-Atlas |

Attribution: every formula below is the upstream authors' work. Our code only
*runs* theirs (`scripts/competition_v2/run_atlasshift_c0.sh`) and measures it.

---

## 1. One-paragraph summary

AtlasShift predicts each target-context perturbation as **the target context's own
control profile, multiplied by a per-gene fold change borrowed from direct
measurements of the same knockdown in other cell lines**. The fold change is a
weighted average across four single-cell atlases (K562 genome-wide Perturb-seq, X-Atlas
HCT116, X-Atlas HEK293T, the complete public VCC 2025 H1 data) plus published CD4
T-cell DE statistics. Each source's effect is **centred on that source's own mean
perturbation response** before fusion, then **shrunk** (amplitude 0.6 / 0.3), clipped at
±3, and applied to the target control mean. Genes whose TSS lies within 5 kb of the
target's TSS are additionally **capped** (a CRISPRi promoter-bleed prior). Integer
counts are produced by a deterministic dual-moment fit on a smoothed control
template. There is **no learned model, no context adaptation, and no context mean
shift**.

## 2. Source response representation (`prepare.py`, `model.load_xatlas`)

Every single-cell source (K562, HCT116, HEK293T, H1) is reduced to per-target sufficient
statistics on a unioned gene axis:

* `target_count_sums[t, g]` — raw counts summed over the target's cells;
* `target_mean_cpm[t, g]` — the **mean of per-cell CPM** (cell-weighted);
* `matched_control_probability[t, g]` — control composition matched to the target's
  batches (K562: `gem_group`, weighted by the target's library mass per batch; X-Atlas:
  `sample`, same weighting; H1: the control of the same 2025 split);
* `matched_control_mean_cpm[t, g]` — the matched control's mean per-cell CPM (weighted by
  the target's cell count per batch);
* `n_cells[t]`.

K562 reads the full raw genome-wide file (`K562_gwps_raw_singlecell_01.h5ad`), keeping
only the 300 official targets + the 150 H1-training targets (`source_targets`). X-Atlas
keeps the same target set and at most **250 control cells per sample** (first by
`cell_integer_id`), and requires ≥ 99 % of the expected expression records to be found.
H1 concatenates train / validation / test (300 targets; each target in exactly one
split) with **its own split's** control as the matched control.

CD4 (`prepare_cd4`) is different in kind: it is the publisher's **DE table**
(`log_fc`, `adj_p_value`, `lfcSE`) per culture condition, averaged across rows for the
same (condition, target), with a strict quality flag: ≥ 2 guides, not a single-guide
estimate, on-target knockdown significant, no distal off-target flag, not low target
expression — **all** rows must pass.

**Pseudocount shrinkage.** `load_xatlas(prior_counts=100000)` adds 100,000 pseudo-counts
of the matched control composition to each target's count sum before normalising, and
blends mean CPM toward the control by `fraction = N / (N + 100000)` where `N` is the
target's total counts. A target with 10^5 total UMIs is therefore pulled halfway to
its control; with 400 cells at 20k UMIs (8·10^6) the shrinkage is ≈ 1 %. It is an
empirical-Bayes-like **reliability shrinkage by sequencing mass**, not a learned
parameter.

## 3. Effect in two spaces (`aligned_effect`)

Two effect representations are computed for every (target, gene):

| space | formula | used for |
|---|---|---|
| `log2fc` | `log2((mean_cpm_pert + 1) / (mean_cpm_ctrl + 1))` | the per-cell mean-CPM moment |
| `bulk_delta` | `log1p(5e4 · p_pert) − log1p(5e4 · p_ctrl)` | the pseudobulk moment (the exact PDS / MSE space, `TS = 5e4`) |

A target is used from a source only if `n_cells ≥ 20` (`minimum_cells`), and only genes the
source measures (`measured_genes`) count toward that source's coverage.

**Centering (`common_subtract = 1`).** For each source, the **mean effect over all the
source's own targets** (excluding each target's own gene) is subtracted from every
target's effect. This removes the source's generic perturbation response (its
"context main effect", in our decomposition language `mu + alpha_source`) and keeps only
the target-specific part — the same move as our `centred_source_beta`, applied per
atlas. Note the centering set is the source's *retained* targets (≈ the 300 official +
150 H1-train targets), not its whole screen.

## 4. Source fusion and weights (`fuse_source_centered`, `add_cd4_family`)

For each official target and gene:

```
effect = Σ_s w_s · 1[covered_s] · effect_s  /  Σ_s w_s · 1[covered_s]
```

with **w = K562 2, HCT116 1, HEK293T 1, H1 2** (`predict.WEIGHTS`). A source that did
not measure the target contributes nothing to numerator or denominator, so a target
measured in one source gets that source's effect unshrunk by fusion.

CD4 is then folded in as a further family with weight **0.5** (per culture condition
the CD4 log2FC is centred on the source-wide mean across targets, `center_scope="source"`,
conditions averaged). In `bulk_delta` space the CD4 log2FC is converted through the
target context's own control composition first.

No weight, amplitude, or threshold is fitted anywhere in the released code; they are
constants. The repository contains no record of how they were chosen (plausibly on
the public leaderboard — that is unknowable from source and is recorded as a risk:
constants chosen against the A/B/C leaderboard would transfer to D/E/F only as far as
A/B/C resemble D/E/F).

## 5. Context adaptation, amplitude, clipping, conversion (`desired_mean`)

There is **no context adaptation** beyond using the target context's own control as
the multiplicative base. Per context, per space:

```
effect'  = clip(amplitude · effect, −3, +3)
log2fc:      p_target = ctrl_mean_cpm · 2^effect'            amplitude 0.6
bulk_delta:  p_target = expm1(max(log1p(5e4 · ctrl_bulk) + effect', 0))   amplitude 0.3
```

each renormalised to sum to one. The context's **mean perturbation response is never
added**: a target with no source measurement receives effect 0, i.e. its prediction
equals the context control (modulo the promoter cap). This is the opposite design
choice to V1, which added a (shrunk) context main effect `m_hat` to every target and a
perturbation-specific term only to 86.

## 6. Promoter-neighbour prior (`prepare_promoters`, `pairs`, `apply_promoter_prior`)

* **Source of pairs.** GENCODE v47 `gene` records (`gencode.v47.annotation.gtf.gz`,
  SHA-256 `df11938c…`). TSS = start on `+`, end on `−`. Gene symbols duplicated in the GTF
  are dropped (`~duplicated(keep=False)`).
* **Rule.** For each official target, every panel gene on the same chromosome with
  `|TSS_target − TSS_neighbour| ≤ 5000 bp` (excluding the target itself). A `divergent`
  flag is recorded but **not used**.
* **Formula.** `remaining(d) = 0.15 + 0.85 · clip(log10(max(d, 500) / 500), 0, 1)`: 15 % of
  control at ≤ 500 bp, rising log-linearly to 100 % at 5 kb. The neighbour's predicted
  composition is capped at `ctrl[j] · remaining(d)` (only lowered, never raised), then the
  row is renormalised.
* **Biology.** CRISPRi (dCas9-KRAB) deposits H3K9me3 over ~1–2 kb around the guide, so
  genes sharing a bidirectional promoter or with a TSS within a few kb of the target
  are silenced along with it. This is a direct, mechanistic, *assay* effect — not a
  statement about gene-regulatory biology — and it is measurable in every CRISPRi screen.
* **Leakage.** None: it reads the genome annotation and the official target / gene lists
  only.
* Counts affected on the Arc panel: 75 pairs, 74 targets, 75 neighbour genes, 40 within 500 bp (`atlas_shift_source_audit.md` §3). Public ablation: main report §7.

This prior is **categorically different from our failed STRING / DepMap / pathway
priors**. Those predicted a *never-measured perturbation's* downstream transcriptional
effect from functional annotation. The promoter prior predicts only a
chromatin-proximity side effect of the CRISPRi assay on a handful of neighbours.

## 7. Integer cell construction (`control_template`, `dual_moment_counts`)

1. Draw `400 × 4` control cells without replacement (seed `20260910 + context index`),
   sort by library size, and average their **compositions in groups of 4** → 400 smoothed
   template cells; each template cell's depth = the rounded mean of its 4 donors' depths.
2. Iteratively (≤ 100 steps, 2 · 10^-4 tolerance) rescale the template so that (a) the
   mean per-cell composition equals the `log2fc`-space target and (b) the depth-weighted
   composition equals the `bulk_delta`-space target — two moments at once, via a
   per-gene exponential tilt along relative depth.
3. Floor `x · depth`, distribute each cell's residual counts by systematic sampling
   (seeded per `context:target`), then **repair column totals** so every gene's
   pseudobulk total equals its expectation to the nearest integer.

Consequences: pseudobulk and mean-CPM moments are hit almost exactly (deterministic
up to one count per gene), cell depths equal the template's, and per-cell variance is
the variance of 4-cell-averaged controls — **lower than real cells**. With 400 cells per
group, averaging raises Wilcoxon power for any shifted gene: a small mean shift
applied to low-noise cells is called significant more readily than the same shift in
real cells. That is directly relevant to the DE yield members (§ comparison report).

## 8. Packaging (`compact.py`, `pack.py`)

`compact.py` losslessly permutes observations to (context, cell, target) order and
writes CSC with gzip, verifying row hashes, sums, nnz and 256 dense rows. `pack.py`
calls `vcc.prep.run_prep` from the installed `vcc-cli 0.2.0` with memory-mapped
arrays — the official checks, no bypass.

## 9. Relation to our frozen research results

* AtlasShift's gain comes from **expanded direct evidence**: the same knockdown measured
  in more cell lines. Our research asked whether a perturbation's effect can be
  predicted from annotation priors with **no** direct measurement anywhere; that
  answer (no, on two external datasets) is untouched by anything AtlasShift does.
* Its centred-per-source fusion is the conserved-effect (`beta`) half of the
  decomposition we studied; it makes no claim about the context × perturbation
  interaction (`gamma`), which it ignores.
* The promoter prior is an assay artefact prior, not a biological unseen-effect prior.
