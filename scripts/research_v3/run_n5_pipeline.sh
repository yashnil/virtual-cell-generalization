#!/bin/bash
# Frozen N5 pipeline (reports/n5_protocol.md §8). Refuses to run if the protocol or
# N5 code differs from the recorded digest.
set -euo pipefail
cd "$(dirname "$0")/../.."
DIGEST=data/provenance/research_v3/n5_protocol_digest.txt
grep -E '^[0-9a-f]{64}  ' "$DIGEST" | shasum -a 256 -c - || { echo "N5 digest mismatch: refusing to run"; exit 2; }
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 VECLIB_MAXIMUM_THREADS=2
mkdir -p outputs/n5
uv run python scripts/research_v3/build_n5_gwps.py
uv run python scripts/research_v3/run_n5.py
uv run python scripts/research_v3/analyse_n5.py
(cd outputs/n5 && find . -type f ! -name manifest_sha256.txt | sort | xargs shasum -a 256 > manifest_sha256.txt)
echo "N5 complete"
