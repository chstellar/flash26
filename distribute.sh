#!/bin/bash
#
#SBATCH --job-name=extendor_distribution
#SBATCH --partition=horence
#SBATCH --time=0-01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
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
elif command -v conda >/dev/null 2>&1; then
  conda activate "${FLASH_PY_ENV:-biopython_env-R}"
fi

PYTHON="${PYTHON:-python}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---- Edit these defaults once, then run: sbatch distribute.sh ----
# The SATC table is hardcoded to the 260714 run requested here.
PROJECT_DIR="${PROJECT_DIR:-$SCRIPT_DIR/flash26}"
SATC="${SATC:-$PROJECT_DIR/results/260714-00-3ants-challenge/filter1/noCluster/target1/2000-clusters/all_satc_merged.txt}"
ANCHOR_SATC="${ANCHOR_SATC:-$PROJECT_DIR/results/260714-00-3ants-challenge/filter1/noCluster/target1/2000-clusters/anchor_satc_merged.txt}"

# Put your extendor list/table and sample partition sheet here, or override them
# at submit time, for example:
#   sbatch --export=ALL,INPUT_TSV=my_extendors.tsv,PARTITION_SHEET=my_partitions.csv distribute.sh
INPUT_TSV="${INPUT_TSV:-extendors.tsv}"
PARTITION_SHEET="${PARTITION_SHEET:-sample_partitions.tsv}"
OUTPUT_TSV="${OUTPUT_TSV:-extendor_partition_distribution.tsv}"

# Column settings. These can be header names or 1-based column indices.
EXTENDOR_COL="${EXTENDOR_COL:-1}"
SAMPLE_COL="${SAMPLE_COL:-1}"
MAJOR_COL="${MAJOR_COL:-2}"
MINOR_COL="${MINOR_COL:-3}"

# Header and delimiter handling.
INPUT_HAS_HEADER="${INPUT_HAS_HEADER:-auto}"
PARTITION_HAS_HEADER="${PARTITION_HAS_HEADER:-auto}"
PARTITION_DELIMITER="${PARTITION_DELIMITER:-auto}"
SATC_HAS_HEADER="${SATC_HAS_HEADER:-auto}"

# Existing code in this repo treats an extendor as target+anchor.
EXTENDOR_ORDER="${EXTENDOR_ORDER:-target-anchor}"

# Optional:
#   LONG_OUTPUT=extendor_partition_distribution.long.tsv
LONG_OUTPUT="${LONG_OUTPUT:-}"

usage() {
  echo "Usage: sbatch distribute.sh [extra distribute.py args...]"
  echo
  echo "This wrapper is configured at the top of distribute.sh."
  echo "Default input/output:"
  echo "  INPUT_TSV=$INPUT_TSV"
  echo "  PARTITION_SHEET=$PARTITION_SHEET"
  echo "  OUTPUT_TSV=$OUTPUT_TSV"
  echo "  SATC=$SATC"
  echo
  echo "Common optional environment variables:"
  echo "  EXTENDOR_COL    extendor column name or 1-based index; default 1"
  echo "  SAMPLE_COL      sample column in partition sheet; default 1"
  echo "  MAJOR_COL       major partition column; default 2"
  echo "  MINOR_COL       minor partition column; default 3"
  echo "  PARTITION_HAS_HEADER auto|yes|no; default auto"
  echo "  PARTITION_DELIMITER  auto|tab|comma; default auto"
  echo "  LONG_OUTPUT     optional long-format output TSV"
  echo
  echo "Example:"
  echo "  sbatch --export=ALL,INPUT_TSV=extendors.tsv,PARTITION_SHEET=sample_partitions.csv,OUTPUT_TSV=dist.tsv distribute.sh"
}

[[ "${1:-}" =~ ^(-h|--help)$ ]] && { usage; exit 0; }

args=(
  --input "$INPUT_TSV"
  --partition_tsv "$PARTITION_SHEET"
  --output "$OUTPUT_TSV"
  --satc "$SATC"
  --anchor_satc "$ANCHOR_SATC"
  --extendor_col "$EXTENDOR_COL"
  --sample_col "$SAMPLE_COL"
  --major_col "$MAJOR_COL"
  --minor_col "$MINOR_COL"
  --input_has_header "$INPUT_HAS_HEADER"
  --partition_has_header "$PARTITION_HAS_HEADER"
  --partition_delimiter "$PARTITION_DELIMITER"
  --satc_has_header "$SATC_HAS_HEADER"
  --extendor_order "$EXTENDOR_ORDER"
)

if [[ -n "$LONG_OUTPUT" ]]; then
  args+=(--long_output "$LONG_OUTPUT")
fi

"$PYTHON" "$SCRIPT_DIR/distribute.py" "${args[@]}" "$@"
