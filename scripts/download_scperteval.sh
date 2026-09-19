#!/usr/bin/env bash
# Download the four audited scPertEval datasets. Resumable and idempotent:
# a file whose size already matches the audited byte count is left untouched.
set -uo pipefail

BASE="https://storage.googleapis.com/scperteval/processed"
DEST="${1:-data/raw/scperteval}"
mkdir -p "$DEST"

# dataset:expected_bytes  (from reports/scperteval_four_context_data_spec.md)
FILES=(
  "replogle22k562:2430512332"
  "replogle22rpe1:1877364555"
  "nadig25hepg2:1236448196"
  "nadig25jurkat:2004474709"
)

status=0
for entry in "${FILES[@]}"; do
  ds="${entry%%:*}"; want="${entry##*:}"
  f="$DEST/${ds}_processed_complete.h5ad"
  if [ -f "$f" ]; then
    have=$(wc -c < "$f" | tr -d ' ')
    if [ "$have" = "$want" ]; then
      echo "SKIP  $ds already complete ($have bytes)"
      continue
    fi
    echo "RESUME $ds at $have / $want bytes"
  else
    echo "START  $ds ($want bytes)"
  fi
  curl -fSL --retry 8 --retry-delay 5 --retry-all-errors -C - \
       --connect-timeout 30 -o "$f" "$BASE/${ds}_processed_complete.h5ad"
  rc=$?
  have=$(wc -c < "$f" 2>/dev/null | tr -d ' ' || echo 0)
  if [ "$rc" != "0" ] || [ "$have" != "$want" ]; then
    echo "FAIL  $ds rc=$rc size=$have want=$want"; status=1
  else
    echo "OK    $ds $have bytes"
  fi
done
exit $status
