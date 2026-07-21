#!/usr/bin/env bash
#SBATCH --job-name=anchor_expression
#SBATCH --partition=horence,owners,normal
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=logs/anchor_expression.%j.out
#SBATCH --error=logs/anchor_expression.%j.err

set -euo pipefail

FLASH_CONDA="${FLASH_CONDA:-/oak/stanford/groups/horence/chester/dabs_ref/miniforge3}"
FLASH_PY_ENV="${FLASH_PY_ENV:-$FLASH_CONDA/envs/biopython_env-R}"
if [[ -f "$FLASH_CONDA/bin/activate" ]]; then
  source "$FLASH_CONDA/bin/activate" "$FLASH_PY_ENV"
fi
PYTHON="${PYTHON:-python}"

join_csv() { local IFS=,; echo "$*"; }

tc7_samples() {
  local rows=(
    "02:GTAGAG:ATCACG:AGTCAA" "04:GGTAGC:CGATGT:AGTTCC" "06:ATGAGC:TTAGGC:ATGTCA"
    "08:CAAAAG:TGACCA:CCGTCC" "10:CAACTA:GTCCGC:ACAGTG" "12:CACGAT:GTGAAA:GCCAAT"
    "14:CACTCA:GTGGCC:CAGATC" "16:CAGGCG:GTTTCG:ACTTGA" "18:CATGGC:GATCAG:CGTACG"
    "20:CATTTT:TAGCTT:GAGTGG" "22:CCAACA:GGCTAC:ACTGAT" "24:TAATCG:CTTGTA:ATTCCT"
  )
  local out=() n a b c
  for row in "${rows[@]}"; do
    IFS=: read -r n a b c <<<"$row"
    out+=("TC7-${n}A_${a}" "TC7-${n}B_${b}" "TC7-${n}C_${c}")
  done
  join_csv "${out[@]}"
}

manipulation_samples() {
  local out=() rep stage stage_id prefix tissue
  for rep in {1..5}; do out+=("0_CH${rep}" "0_CT${rep}"); done
  for stage in "1:BLT25:OLT25" "2:BLM:OLM" "3:BDM:ODM"; do
    IFS=: read -r stage_id b_prefix o_prefix <<<"$stage"
    for prefix in "$b_prefix" "$o_prefix"; do
      for tissue in H T; do
        for rep in {1..5}; do out+=("${stage_id}_${prefix}${tissue}${rep}"); done
      done
    done
  done
  join_csv "${out[@]}"
}

usage() {
  echo "Usage: sbatch expression.sh ANCHOR [tc7|manipulation|auto|SAMPLE_CSV] [OUTPUT_PREFIX] [extra expression.py args...]"
}

[[ "${1:-}" =~ ^(-h|--help)$ ]] && { usage; exit 0; }
ANCHOR="${1:-${ANCHOR:-}}"; [[ -n "$ANCHOR" ]] || { usage; exit 2; }
PRESET="${2:-${SAMPLE_PRESET:-tc7}}"
OUT="${3:-results/260714-00-3ants-challenge/filter1/noCluster/target1/2000-clusters/expression_${ANCHOR}_${PRESET}}"
if [[ $# -ge 3 ]]; then shift 3; elif [[ $# -ge 2 ]]; then shift 2; else shift 1; fi

SATC="${SATC:-results/260714-00-3ants-challenge/filter1/noCluster/target1/2000-clusters/all_satc.filtered.dump}"
[[ -s "$SATC" ]] || SATC="results/260714-00-3ants-challenge/filter1/noCluster/target1/2000-clusters/all_satc_merged.txt"

PLOT_MODE="${PLOT_MODE:-plain}"
ORDER_FLAG=()
case "$PRESET" in
  tc7) SAMPLE_ORDER="$(tc7_samples)"; PLOT_MODE="${PLOT_MODE_OVERRIDE:-tc7}"; ORDER_FLAG=(--only_ordered_samples) ;;
  manipulation|infection) SAMPLE_ORDER="$(manipulation_samples)"; PLOT_MODE="${PLOT_MODE_OVERRIDE:-manipulation}"; ORDER_FLAG=(--only_ordered_samples) ;;
  auto|natural) SAMPLE_ORDER=auto; PLOT_MODE="${PLOT_MODE_OVERRIDE:-plain}" ;;
  *) SAMPLE_ORDER="$PRESET"; PLOT_MODE="${PLOT_MODE_OVERRIDE:-$PLOT_MODE}"; [[ "${ONLY_ORDERED_SAMPLES:-1}" == 1 ]] && ORDER_FLAG=(--only_ordered_samples) ;;
esac

mkdir -p "$(dirname "$OUT")" logs
"$PYTHON" expression.py \
  --anchor "$ANCHOR" --satc "$SATC" --sample_order "$SAMPLE_ORDER" \
  --output_prefix "$OUT" --plot_mode "$PLOT_MODE" \
  --grid_cols "${GRID_COLS:-3}" --grid_rows "${GRID_ROWS:-3}" \
  "${ORDER_FLAG[@]}" "$@"
