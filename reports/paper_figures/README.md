# Paper figures (frozen publication set)

Four main figures and five Extended Data figures for the manuscript in `paper/`. Every panel is drawn from a small
table in `data/figure_sources/paper/` that is built only from frozen artefacts. No experiment was re-run; no
estimator, threshold, preregistered result or numerical value was changed; no frozen file was written.

* Claims, evidentiary status and constraints: [`figure_design_spec.md`](figure_design_spec.md)
* Layout, rendering and integrity checks: [`figure_qc.md`](figure_qc.md) (machine output: `figure_qc_auto.json`,
  `validation.json`)
* Captions: `fig*_caption.md`, `ext_fig*_caption.md`
* Per-figure manuscript notes: [`../../paper/figure_notes.md`](../../paper/figure_notes.md)

## Rebuild and validate

```bash
bash scripts/paper_figures/run_all.sh
```

This takes about 25 s. It runs four steps:

1. **Source builders:** each one checks every frozen input against the manifest or freeze file that pins it and
   stops on a mismatch.
2. **Plot scripts.**
3. **`qc_figures.py`:** font floor, clipping, text–text overlap and text-on-data checks.
4. **`validate_figures.py`:**
   * full re-hash of the N1/N4, N3, N5 and N6 output manifests and the decomposition/N5 freeze files;
   * source-table and frozen-input hashes against their provenance sidecars;
   * figure SVG/PNG/PDF and source hashes against their manifests;
   * captions present.

Shared style and semantics: `scripts/paper_figures/style.py`. Provenance and saving: `scripts/paper_figures/_common.py`.

## Evidentiary status labels

| label | meaning |
|---|---|
| **Confirmatory** | preregistered primary endpoint or preregistered control, evaluated by its frozen rule |
| **Prereg. replication** | the same frozen protocol re-run on new contexts (N3-A) |
| **Prereg. quantity, descriptive** | a quantity computed under a preregistered protocol, read without a preregistered test |
| **Exploratory** | defined after results existed (labelled in the panel and caption) |
| **Background** | replication of published analyses (Molina & Zhang decomposition) |
| **Descriptive aggregate** | cross-context median; no inferential uncertainty exists |
| **Schematic** | design illustration; no data |

## Main figures

All figures are 183 mm wide. Outputs are SVG (editable text), 600 dpi PNG and PDF (embedded TrueType fonts), with a
6 pt minimum text size.

### Fig. 1 — Conserved and context-specific perturbation effects differ in reproducibility

`fig1_transfer_components.*` · source `fig1_components.csv` · upstream `outputs/four_context_v1/summary.json`,
`outputs/four_context_sensitivity/variant_table.csv`

| panel | content | status |
|---|---|---|
| A | zero-shot vs target-calibration information boundary | Schematic |
| B | δ = μ + α + β + γ, following Molina & Zhang (2026) | Schematic / Background |
| C | noise-corrected energy shares (β 30.1 %, γ 21.0 %) with 21-variant ranges | Background (predeclared decomposition) |
| D | split-half reproducibility (β 80.8 %, γ 49.5 %) on the same rows | Background |

### Fig. 2 — Context-specific perturbation effects require substantially more target data to estimate

`fig2_sample_complexity.*` · sources `fig2_curves.csv`, `fig2_gain_fraction.csv` · upstream `outputs/n3/n3a/`,
`outputs/n3/n3_decision.json`

| panel | content | status |
|---|---|---|
| A | template + scale recovery (normalised to each context's own k_ref gain) vs k | Prereg. replication (k_T50 construct); median line descriptive |
| B | absolute reliable γ⊥ recovery (M3) vs k | Prereg. replication; median line descriptive |
| C | six held-out contexts, perturbation-bootstrap 95 % CIs | Prereg. replication |
| D | G(k) at k = 20/50/100, γ⊥ vs template + scale; C3 bar | ratio point confirmatory at k = 20 (C3); draw intervals and median row descriptive |

### Fig. 3 — Broader source panels improve zero-shot transfer across most target contexts

`fig3_source_count.*` · source `fig3_source_count.csv` (k = 0 rows) · upstream `outputs/n3/n3b/`

| panel | content | status |
|---|---|---|
| A | source count m × target budget k design (this figure: k = 0) | Schematic |
| B | zero-shot template-removed full response (M1) vs m, per target + median | Prereg. quantity, descriptive; median descriptive |
| C | per-target m = 2 → m = 5 change (the same values as B) | Prereg. quantity, descriptive |

No own-frame or anchor-based (k > 0) analysis appears in the main figure; those are in ED 5.

### Fig. 4 — Source choice strongly influences perturbation transfer

`fig4_source_compatibility.*` · sources `fig4_compatibility.csv`, `fig4_robustness.csv`,
`fig4_per_perturbation.csv` · upstream `outputs/n5/`

| panel | content | status |
|---|---|---|
| A | C_S (noise-corrected latent cosine) per source with 95 % CIs; median source reliability | Confirmatory (primary C_S; others preregistered secondary) |
| B | GWPS − RPE1 advantage under every preregistered control, by estimand | Confirmatory |
| C | ECDF of per-perturbation r difference | Prereg., descriptive |
| D | same cell line × same study grid; N6 cell marked gate-failed | Schematic (design facts from N5/N6) |

## Extended Data

| figure | panels | status | sources | upstream |
|---|---|---|---|---|
| **ED 1** `ext_fig1_robustness.*` — decomposition robustness | A energy shares, 21 variants · B β/γ reproducibility · C γ before/after correction · D reliability vs depth | Background (predeclared battery) | `ext1_variants.csv`, `ext1_depth.csv` | `outputs/four_context_sensitivity/` |
| **ED 2** `ext_fig2_n1_controls.*` — validity controls for Fig. 2 | A anchor-permutation null · B shared-control check · C single-anchor harm · D N1 vs N3-A | Confirmatory failure checks; D intervals descriptive | `ext2_*.csv` | `outputs/n3/n3a/`, `outputs/n1_n4/n1/` |
| **ED 3** `ext_fig3_agreement_correction.*` — agreement correction | A v1 vs corrected ρ · B partial Spearman (T2) · C T1 vs baselines | A exploratory; B–C confirmatory (N4 WEAK) | `ext3_*.csv` | `data/figure_sources/n1_n4/`, `outputs/n3/n4/` |
| **ED 4** `ext_fig4_n6_reliability.*` — N6 reliability gate | A pooled reliability vs gate · B per-perturbation reliability · C positive control · D gate checklist; VIPerturb compatibility omitted | Gate confirmatory (Outcome D); two pooled GWPS values exploratory | `ext4_*.csv` | `outputs/n6/` |
| **ED 5** `ext_fig5_n3_reference_frame.*` — source count with anchors | **A preregistered Q1 Δ(20), fixed frame (primary) vs own frame** · B own-frame calibration gain vs m (former Fig. 3C) · C fixed vs own frame curves · D frame-free full response, zero-shot vs 20 anchors · E own-frame k × m heatmaps (former Fig. 3D) | A fixed frame confirmatory (criterion fired); A own frame, B, C own-frame curve and E exploratory; D prereg. quantities, descriptive | `ext5_delta20.csv`, `ext5_reference_frame.csv`, `fig3_source_count.csv` | `outputs/n3/n3b/`, `outputs/n3/n3_decision.json` |

All Extended Data source tables are built by `build_ext_sources.py`.

## Colour semantics (identical in every figure)

| meaning | colour |
|---|---|
| template / scale | teal `#2A9D78` (squares) |
| conserved β | blue `#1F6FB4` (circles) |
| context-specific γ, γ⊥ | vermilion `#D55E00` (diamonds) |
| target-calibrated / emphasised | near-black `#1A1A1A` (filled) |
| zero-shot / source-only / other sources | grey `#8F8F8F` (open) |
| failed gate, non-interpretable | light grey `#BDBDBD`, hatched |
| γ magnitude heatmap | white → vermilion → brown (single hue) |

## Notes

* File names (`fig4_source_compatibility`, `fig4_compatibility.csv`, `fig3_source_count`) predate the final titles.
  They are kept so that manifests stay stable. "Compatibility" in a file name refers to the N5 protocol's name for
  C_S, not a biological claim.
* `scripts/paper_figures/*.py` carry a file-level `ruff: noqa: E501` for label strings; all other lint rules pass.
* The original figure pipeline (`scripts/figures/`, `data/figure_sources/*.csv`, `reports/figures/`) is untouched.
