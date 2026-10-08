# Figure QC

Date 2026-10-08 (final freeze). All checks after the final rebuild: `bash scripts/paper_figures/run_all.sh`.

## How it was checked

1. **Automated layout check** (`scripts/paper_figures/qc_figures.py` → `figure_qc_auto.json`). It rebuilds every
   figure and reports min/max font size, text beyond the canvas, and text–text overlaps (> 2 px in both directions).
   It also flags text sitting on data, meaning any data-line vertex or marker inside a text box. For annotations
   it uses the text box only, not the leader line. It does not test line segments between vertices.
2. **Visual inspection** of every PNG at ~1/3 scale (about printed double-column size) after every edit. Text-on-data
   collisions found this way were fixed:
   * Fig 2: A annotation, D legend;
   * Fig 3: C overplotting, the Δ column, B legend replaced by direct labels;
   * Fig 4: B/C labels, the D cell overflow;
   * ED 1, ED 2, ED 4 and ED 5 legends.
3. **Independent renders.**
   * Every PDF was rasterised with poppler `pdftoppm` at 100 dpi (print scale) and every SVG with QuickLook; all
     render.
   * PDF fonts are embedded TrueType: Arial regular/bold/italic. STIX General is used for the ⊥ glyph, and
     DejaVu Sans as a fallback for a few symbols (▼, ◀, →) in two figures.
   * SVG text stays editable text.
4. **Grayscale** check of Figs 1, 2 and 4. Every critical distinction has a second channel:
   * marker shape (square / circle / diamond for template / β / γ);
   * fill (zero-shot open vs calibrated filled);
   * line style (full vs depth-matched ECDF);
   * direct labels.
5. **Integrity (`validate_figures.py` → `validation.json`).**
   * Full re-hash of every entry in `outputs/{n1_n4,n3,n5,n6}/manifest_sha256.txt`,
     `data/provenance/scperteval/{decomposition_phase,canonical_v1}_freeze.txt` and
     `data/provenance/research_v3/n5_freeze_sha256.txt`: 160 entries, 0 mismatches, 0 missing.
   * All 22 source tables match their provenance sidecars (table hash and every frozen-input hash).
   * All 9 figure manifests match (SVG/PNG/PDF and source hashes); every caption is present.
   * `git status` shows only the three new directories; no tracked or frozen file changed.
   * Every manifest-pinned input re-hashed equal at build time.
   * Derived quantities that also exist in frozen files were asserted equal: k_T50, the C3 ratios, the N3-B
     own-frame medians, the N5 fractions 0.858 / 0.705, and the depth-experiment medians.
   * The ED 5 fixed-frame Δ(20) values are copied from `n3_decision.json`.
6. **Lint:** `ruff check` and `ruff format --check` pass. E501 is exempted per file for label strings.

## Per figure

| figure | size (mm) | PNG px (600 dpi) | min / max font (pt) | clipped | text overlaps / text on data | SVG / PDF | colour semantics | uncertainty shown |
|---|---|---|---|---|---|---|---|---|
| Fig 1 | 183 × 108 | 4320 × 2550 | 6.0 / 12.5 | 0 | 0 / 0 | ok / ok | consistent | split-half s.d. (below marker size, stated); 21-variant ranges labelled as ranges |
| Fig 2 | 183 × 141 | 4320 × 3330 | 6.0 / 8.5 | 0 | 0 / 0 | ok / ok | consistent | C: perturbation-bootstrap 95 % CI; D: draw intervals; A/B and D median row descriptive (stated) |
| Fig 3 | 183 × 66 | 4320 × 1554 | 6.0 / 8.5 | 0 | 0 / 0 | ok / ok | consistent | B, C: range over 5 cell-split repeats; median descriptive |
| Fig 4 | 183 × 124 | 4320 × 2940 | 6.0 / 8.5 | 0 | 0 / 0 | ok / ok | consistent | A/B: 95 % bootstrap CIs; C: full distribution |
| ED 1 | 183 × 99 | 4320 × 2340 | 6.0 / 8.5 | 0 | 0 / 0 | ok / ok | consistent | D: 95 % CI over pairs; A–C every variant shown |
| ED 2 | 183 × 117 | 4320 × 2760 | 6.0 / 8.5 | 0 | 0 / 0 | ok / ok | consistent | draw intervals (A, D); C band = range of contexts |
| ED 3 | 183 × 109 | 4320 × 2580 | 6.0 / 8.5 | 0 | 0 / 0 | ok / ok | consistent | B, C bootstrap CIs; A exploratory point values |
| ED 4 | 183 × 107 | 4320 × 2520 | 6.0 / 8.5 | 0 | 0 / 0 | ok / ok | consistent; VIPerturb compatibility omitted | C bootstrap CIs; A single statistics |
| ED 5 | 183 × 222 | 4320 × 5244 | 6.0 / 8.5 | 0 | 0 / 0 | ok / ok | consistent | A, B: 95 % draw intervals; C band: draw interval; D repeat range / draw intervals; E none (stated) |

## Remaining limitations

* **Different normalisations in Fig 2A and 2B.** A is a fraction of each context's own k_ref gain; B is absolute.
  This is deliberate and is stated in the panel subtitles and in the first sentence of the caption. D is the
  like-for-like comparison.
* **No inferential uncertainty on cross-context medians.** None exists in the frozen analysis; the medians are
  labelled descriptive.
* **Fig 3 has no inferential interval** for the m-trend; the repeat ranges show cell-split variation only.
* **ED 5 is a full-page figure (222 mm).** It is within full-page limits but dense, and its heatmap cells carry no
  uncertainty by design (intervals for k = 20/50 are in panel B).
* **Schematics are drawn programmatically.** They use one geometry system, but a final pass by a professional
  illustrator could still refine Fig 1A.
