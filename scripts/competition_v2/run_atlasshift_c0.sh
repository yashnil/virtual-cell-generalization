#!/usr/bin/env bash
# C0_ATLASSHIFT_REPRODUCTION: the exact command sequence used in the competition-v2 phase.
#
# Runs the vendored AtlasShift code (third_party/atlasshift, commit d24ce4f, MIT) UNMODIFIED,
# in its own pinned Python 3.13 environment. The only deviation from upstream is
# --atlas-root pointing at our local, revision-pinned X-Atlas/Orion copy instead of the
# live Hugging Face dataset. It packages a .vcc and NEVER submits it.
#
# Prerequisites: scripts/competition_v2/download_sources.py has completed (checksums
# verified), and `vcc` (vcc-cli 0.2.0) is on PATH.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
UP="$ROOT/third_party/atlasshift"
PY="$UP/.venv/bin/python"
C0="$ROOT/outputs/competition_v2/atlasshift_c0"
RAW="$ROOT/data/raw/competition_v2"

if [ ! -x "$PY" ]; then
  uv venv --python 3.13 "$UP/.venv"
  uv pip install --python "$PY" -r "$UP/requirements.txt" slafdb
fi

mkdir -p "$C0/data" "$C0/scratch"
for f in gene_names.csv pert_counts.csv context_A.h5ad context_B.h5ad context_C.h5ad; do
  ln -sf "$ROOT/data/raw/arc2026/controls/$f" "$C0/data/$f"
done

cd "$C0/data"
for src in promoters h1 cd4 hct hek k562; do
  "$PY" "$UP/prepare.py" --data-dir . --raw-dir "$RAW" --source "$src" \
    --atlas-root "$RAW/xatlas_orion/data" > "../prepare_$src.log" 2>&1
done
"$PY" "$UP/predict.py" --data-dir . --output ../prediction.h5ad > ../predict.log 2>&1
cd "$C0"
"$PY" "$UP/compact.py" prediction.h5ad prediction_compact.h5ad > compact.log 2>&1
"$PY" "$UP/pack.py" prediction_compact.h5ad --data-dir data \
  --output c0_atlasshift_reproduction_val.vcc --scratch-dir scratch > pack.log 2>&1

cd "$ROOT"
uv run python scripts/competition_v2/validate_c0_bundle.py
echo "C0 packaged at $C0/c0_atlasshift_reproduction_val.vcc — NOT submitted."
