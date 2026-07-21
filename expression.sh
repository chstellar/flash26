#!/bin/bash
#
#SBATCH --job-name=anchor_expression
#SBATCH --partition=horence
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=chesteryu@stanford.edu

set -euo pipefail

if command -v ml >/dev/null 2>&1; then
  ml purge
elif command -v module >/dev/null 2>&1; then
  module purge
fi
unset PYTHONPATH PYTHONHOME

FLASH_CONDA="${FLASH_CONDA:-/oak/stanford/groups/horence/chester/dabs_ref/miniforge3}"
if [[ -x "$FLASH_CONDA/bin/conda" ]]; then
  eval "$("$FLASH_CONDA/bin/conda" shell.bash hook)"
fi
if command -v mamba >/dev/null 2>&1; then
  eval "$(mamba shell hook --shell bash)"
  mamba activate "${FLASH_PY_ENV:-biopython_env-R}"
else
  conda activate "${FLASH_PY_ENV:-biopython_env-R}"
fi
PYTHON="${PYTHON:-python}"
export MPLBACKEND="${MPLBACKEND:-Agg}"

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
  echo "Usage: sbatch expression.sh ANCHOR_OR_ANCHOR_TXT [tc7|manipulation|auto|SAMPLE_CSV] [OUTPUT_PREFIX] [extra expression.py args...]"
  echo "If ANCHOR_OR_ANCHOR_TXT is a file, it should contain one anchor sequence per row."
}

[[ "${1:-}" =~ ^(-h|--help)$ ]] && { usage; exit 0; }
ANCHOR_INPUT="${1:-${ANCHOR:-}}"; [[ -n "$ANCHOR_INPUT" ]] || { usage; exit 2; }
PRESET="${2:-${SAMPLE_PRESET:-tc7}}"
if [[ -f "$ANCHOR_INPUT" ]]; then
  DEFAULT_OUT="results/260714-00-3ants-challenge/filter1/noCluster/target1/2000-clusters/expression_$(basename "$ANCHOR_INPUT" .txt)_${PRESET}"
else
  DEFAULT_OUT="results/260714-00-3ants-challenge/filter1/noCluster/target1/2000-clusters/expression_${ANCHOR_INPUT}_${PRESET}"
fi
OUT="${3:-$DEFAULT_OUT}"
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

safe_anchor_name() {
  printf '%s' "$1" | tr -c 'A-Za-z0-9._-' '_'
}

run_anchor() {
  local anchor="$1"
  local out_prefix="$2"
  shift 2
  mkdir -p "$(dirname "$out_prefix")" logs
  "$PYTHON" expression.py \
    --anchor "$anchor" --satc "$SATC" --sample_order "$SAMPLE_ORDER" \
    --output_prefix "$out_prefix" --plot_mode "$PLOT_MODE" \
    --grid_cols "${GRID_COLS:-3}" --grid_rows "${GRID_ROWS:-3}" \
    "${ORDER_FLAG[@]}" "$@"
}

if [[ -f "$ANCHOR_INPUT" ]]; then
  while IFS= read -r anchor || [[ -n "$anchor" ]]; do
    anchor="${anchor//$'\r'/}"
    anchor="${anchor%%[[:space:]]*}"
    [[ -n "$anchor" ]] || continue
    [[ "$anchor" == \#* ]] && continue
    run_anchor "$anchor" "${OUT}_$(safe_anchor_name "$anchor")" "$@"
  done < "$ANCHOR_INPUT"
else
  run_anchor "$ANCHOR_INPUT" "$OUT" "$@"
fi
