#!/usr/bin/env bash
set -euo pipefail

ANCHOR="${1:?Usage: bash expression.sh ANCHOR SAMPLE_ORDER_STRING [OUTPUT_PREFIX] [extra expression.py args...]}"
SAMPLE_ORDER="${2:?Usage: bash expression.sh ANCHOR SAMPLE_ORDER_STRING [OUTPUT_PREFIX] [extra expression.py args...]}"
OUTPUT_PREFIX="${3:-results/260714-00-3ants-challenge/filter1/noCluster/target1/2000-clusters/expression_${ANCHOR}}"
if [[ $# -ge 3 ]]; then
  shift 3
else
  shift 2
fi

SATC="${SATC:-results/260714-00-3ants-challenge/filter1/noCluster/target1/2000-clusters/all_satc.filtered.dump}"
if [[ ! -s "$SATC" ]]; then
  ALT="results/260714-00-3ants-challenge/filter1/noCluster/target1/2000-clusters/all_satc_merged.txt"
  if [[ -s "$ALT" ]]; then
    SATC="$ALT"
  fi
fi

PYTHON="${PYTHON:-python3}"

"$PYTHON" expression.py \
  --anchor "$ANCHOR" \
  --satc "$SATC" \
  --sample_order "$SAMPLE_ORDER" \
  --output_prefix "$OUTPUT_PREFIX" \
  "$@"
