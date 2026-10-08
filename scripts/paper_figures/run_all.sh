#!/usr/bin/env bash
# Rebuild every paper-figure source table and figure from frozen artefacts, then run layout QC.
set -euo pipefail
cd "$(dirname "$0")"
for b in build_fig1_sources build_fig2_sources build_fig3_sources build_fig4_sources build_ext_sources; do
  uv run python "$b.py"
done
for p in plot_fig1 plot_fig2 plot_fig3 plot_fig4 plot_ext_fig1 plot_ext_fig2 plot_ext_fig3 plot_ext_fig4 plot_ext_fig5; do
  uv run python "$p.py"
done
uv run python qc_figures.py
uv run python validate_figures.py
