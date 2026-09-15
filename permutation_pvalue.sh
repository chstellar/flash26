#!/bin/bash
#
#SBATCH --partition=horence
#SBATCH --time=0-01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=chesteryu@stanford.edu

set -euo pipefail

ml purge

FLASH_CONDA="${FLASH_CONDA:-/oak/stanford/groups/horence/chester/dabs_ref/miniforge3}"
if [[ -x "$FLASH_CONDA/bin/conda" ]]; then
  eval "$("$FLASH_CONDA/bin/conda" shell.bash hook)"
fi
if command -v mamba >/dev/null 2>&1; then
  eval "$(mamba shell hook --shell bash)"
  mamba activate "${FLASH_PY_ENV:-adelie_env}"
elif command -v conda >/dev/null 2>&1; then
  conda activate "${FLASH_PY_ENV:-adelie_env}"
fi

PROJECT_DIR="${PROJECT_DIR:-/scratch/users/jiamuyu/proj_botryllus/flash}"
PYTHON="${PYTHON:-python}"

# Override these with sbatch --export=ALL,INPUT_DIR=...,METADATA_CATEGORIES=...
INPUT_DIR="${INPUT_DIR:-${PROJECT_DIR}/results}"
METADATA_CATEGORIES="${METADATA_CATEGORIES:-}"
MATRIX="${MATRIX:-test}"                  # test, train, or both
PERMUTATIONS="${PERMUTATIONS:-10000}"
SEED="${SEED:-42}"
OUTPUT_TSV="${OUTPUT_TSV:-}"
PDF="${PDF:-}"
SIDECAR="${SIDECAR:-}"

usage() {
  echo "Usage: sbatch permutation_pvalue.sh [extra permutation_pvalue.py args...]"
  echo
  echo "Required:"
  echo "  INPUT_DIR              directory containing one *_confusion_matrices.pdf"
  echo "  METADATA_CATEGORIES    comma-separated metadata category name(s)"
  echo
  echo "Optional environment variables:"
  echo "  MATRIX=$MATRIX"
  echo "  PERMUTATIONS=$PERMUTATIONS"
  echo "  SEED=$SEED"
  echo "  OUTPUT_TSV=$OUTPUT_TSV"
  echo "  PDF=$PDF"
  echo "  SIDECAR=$SIDECAR"
  echo
  echo "Example:"
  echo "  sbatch --export=ALL,INPUT_DIR=/path/to/results,METADATA_CATEGORIES=species,PERMUTATIONS=100000 permutation_pvalue.sh"
}

[[ "${1:-}" =~ ^(-h|--help)$ ]] && { usage; exit 0; }

if [[ -z "$METADATA_CATEGORIES" ]]; then
  echo "ERROR: set METADATA_CATEGORIES to one or more comma-separated metadata columns." >&2
  usage >&2
  exit 2
fi

args=(
  --input_dir "$INPUT_DIR"
  --metadata_category "$METADATA_CATEGORIES"
  --matrix "$MATRIX"
  --permutations "$PERMUTATIONS"
  --seed "$SEED"
)

if [[ -n "$OUTPUT_TSV" ]]; then
  args+=(--output "$OUTPUT_TSV")
fi
if [[ -n "$PDF" ]]; then
  args+=(--pdf "$PDF")
fi
if [[ -n "$SIDECAR" ]]; then
  args+=(--sidecar "$SIDECAR")
fi

"$PYTHON" "$PROJECT_DIR/permutation_pvalue.py" "${args[@]}" "$@"
