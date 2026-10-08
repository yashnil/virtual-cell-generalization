#!/bin/bash
# Frozen N6 pipeline (reports/n6_protocol.md §8). Refuses to run if the protocol,
# code, axes or converted artefact differ from the recorded digest.
set -euo pipefail
cd "$(dirname "$0")/../.."
DIGEST=data/provenance/research_v3/n6_protocol_digest.txt
grep -E '^[0-9a-f]{64}  ' "$DIGEST" | shasum -a 256 -c - || { echo "N6 digest mismatch: refusing to run"; exit 2; }
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 VECLIB_MAXIMUM_THREADS=2
mkdir -p outputs/n6
(cd scripts/research_v3 && uv run python build_n6.py)
uv run python scripts/research_v3/run_n6.py
uv run python scripts/research_v3/analyse_n6.py
(cd outputs/n6 && find . -type f ! -name manifest_sha256.txt | sort | xargs shasum -a 256 > manifest_sha256.txt)
echo "N6 complete"
