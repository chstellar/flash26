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

TC7_SAMPLES=(
  TC7-02A_GTAGAG
  TC7-02B_ATCACG
  TC7-02C_AGTCAA
  TC7-04A_GGTAGC
  TC7-04B_CGATGT
  TC7-04C_AGTTCC
  TC7-06A_ATGAGC
  TC7-06B_TTAGGC
  TC7-06C_ATGTCA
  TC7-08A_CAAAAG
  TC7-08B_TGACCA
  TC7-08C_CCGTCC
  TC7-10A_CAACTA
  TC7-10B_GTCCGC
  TC7-10C_ACAGTG
  TC7-12A_CACGAT
  TC7-12B_GTGAAA
  TC7-12C_GCCAAT
  TC7-14A_CACTCA
  TC7-14B_GTGGCC
  TC7-14C_CAGATC
  TC7-16A_CAGGCG
  TC7-16B_GTTTCG
  TC7-16C_ACTTGA
  TC7-18A_CATGGC
  TC7-18B_GATCAG
  TC7-18C_CGTACG
  TC7-20A_CATTTT
  TC7-20B_TAGCTT
  TC7-20C_GAGTGG
  TC7-22A_CCAACA
  TC7-22B_GGCTAC
  TC7-22C_ACTGAT
  TC7-24A_TAATCG
  TC7-24B_CTTGTA
  TC7-24C_ATTCCT
)

MANIPULATION_SAMPLES=(
  0_CH1
  0_CH2
  0_CH3
  0_CH4
  0_CH5
  0_CT1
  0_CT2
  0_CT3
  0_CT4
  0_CT5
  1_BLT25H1
  1_BLT25H2
  1_BLT25H3
  1_BLT25H4
  1_BLT25H5
  1_BLT25T1
  1_BLT25T2
  1_BLT25T3
  1_BLT25T4
  1_BLT25T5
  1_OLT25H1
  1_OLT25H2
  1_OLT25H3
  1_OLT25H4
  1_OLT25H5
  1_OLT25T1
  1_OLT25T2
  1_OLT25T3
  1_OLT25T4
  1_OLT25T5
  2_BLMH1
  2_BLMH2
  2_BLMH3
  2_BLMH4
  2_BLMH5
  2_BLMT1
  2_BLMT2
  2_BLMT3
  2_BLMT4
  2_BLMT5
  2_OLMH1
  2_OLMH2
  2_OLMH3
  2_OLMH4
  2_OLMH5
  2_OLMT1
  2_OLMT2
  2_OLMT3
  2_OLMT4
  2_OLMT5
  3_BDMH1
  3_BDMH2
  3_BDMH3
  3_BDMH4
  3_BDMH5
  3_BDMT1
  3_BDMT2
  3_BDMT3
  3_BDMT4
  3_BDMT5
  3_ODMH1
  3_ODMH2
  3_ODMH3
  3_ODMH4
  3_ODMH5
  3_ODMT1
  3_ODMT2
  3_ODMT3
  3_ODMT4
  3_ODMT5
)

join_by_comma() {
  local IFS=,
  echo "$*"
}

usage() {
  cat <<'EOF'
Usage:
  bash expression.sh ANCHOR [tc7|manipulation|auto|SAMPLE_ORDER_STRING] [OUTPUT_PREFIX] [extra expression.py args...]
  sbatch expression.sh ANCHOR [tc7|manipulation|auto|SAMPLE_ORDER_STRING] [OUTPUT_PREFIX] [extra expression.py args...]

Examples:
  sbatch expression.sh AAAAAAGAAAAGAATACAAAACACAAGGGGA tc7
  sbatch expression.sh AAAAAAGAAAAGAATACAAAACACAAGGGGA manipulation
  bash expression.sh AAAAAAGAAAAGAATACAAAACACAAGGGGA "sampleA,sampleB,sampleC" out/custom --plot_mode plain --only_ordered_samples
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

ANCHOR="${1:-${ANCHOR:-}}"
if [[ -z "$ANCHOR" ]]; then
  usage
  exit 2
fi

PRESET_OR_ORDER="${2:-${SAMPLE_PRESET:-tc7}}"
OUTPUT_PREFIX="${3:-results/260714-00-3ants-challenge/filter1/noCluster/target1/2000-clusters/expression_${ANCHOR}_${PRESET_OR_ORDER}}"
if [[ $# -ge 3 ]]; then
  shift 3
elif [[ $# -ge 2 ]]; then
  shift 2
else
  shift 1
fi

SATC="${SATC:-results/260714-00-3ants-challenge/filter1/noCluster/target1/2000-clusters/all_satc.filtered.dump}"
if [[ ! -s "$SATC" ]]; then
  ALT="results/260714-00-3ants-challenge/filter1/noCluster/target1/2000-clusters/all_satc_merged.txt"
  if [[ -s "$ALT" ]]; then
    SATC="$ALT"
  fi
fi

PYTHON="${PYTHON:-python3}"
PLOT_MODE="${PLOT_MODE:-plain}"
ORDER_FLAG=()

case "$PRESET_OR_ORDER" in
  tc7)
    SAMPLE_ORDER="$(join_by_comma "${TC7_SAMPLES[@]}")"
    PLOT_MODE="${PLOT_MODE_OVERRIDE:-tc7}"
    ORDER_FLAG=(--only_ordered_samples)
    ;;
  manipulation|infection)
    SAMPLE_ORDER="$(join_by_comma "${MANIPULATION_SAMPLES[@]}")"
    PLOT_MODE="${PLOT_MODE_OVERRIDE:-manipulation}"
    ORDER_FLAG=(--only_ordered_samples)
    ;;
  auto|natural)
    SAMPLE_ORDER="auto"
    PLOT_MODE="${PLOT_MODE_OVERRIDE:-plain}"
    ;;
  *)
    SAMPLE_ORDER="$PRESET_OR_ORDER"
    PLOT_MODE="${PLOT_MODE_OVERRIDE:-$PLOT_MODE}"
    if [[ "${ONLY_ORDERED_SAMPLES:-1}" == "1" ]]; then
      ORDER_FLAG=(--only_ordered_samples)
    fi
    ;;
esac

mkdir -p "$(dirname "$OUTPUT_PREFIX")" logs

"$PYTHON" expression.py \
  --anchor "$ANCHOR" \
  --satc "$SATC" \
  --sample_order "$SAMPLE_ORDER" \
  --output_prefix "$OUTPUT_PREFIX" \
  --plot_mode "$PLOT_MODE" \
  --grid_cols "${GRID_COLS:-3}" \
  --grid_rows "${GRID_ROWS:-3}" \
  "${ORDER_FLAG[@]}" \
  "$@"
