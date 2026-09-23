#!/bin/bash
#
#SBATCH --partition=horence
#SBATCH --time=0-01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=chesteryu@stanford.edu

set -euo pipefail

ml purge 2>/dev/null || true
FLASH_CONDA="${FLASH_CONDA:-/oak/stanford/groups/horence/chester/dabs_ref/miniforge3}"
eval "$("$FLASH_CONDA/bin/conda" shell.bash hook)"
conda activate adelie_env

PROJECT_DIR="${PROJECT_DIR:-/scratch/users/jiamuyu/proj_botryllus/flash/}"
PYTHON="${PYTHON:-python}"

# ant
# INPUT_DIR="${INPUT_DIR:-${PROJECT_DIR}/results/260819-00-cfloridanus-fungus/filter1/noCluster/hyena/normalized}"
# METADATA_CATEGORIES="${METADATA_CATEGORIES:-fungus_species,tissue}"

# wolbachia
INPUT_DIR="${INPUT_DIR:-${PROJECT_DIR}/results/260826-01-2flies-wolbachia/filter1/shiftDist-levFilter/hyena/normalized}"
METADATA_CATEGORIES="${METADATA_CATEGORIES:-infection_status}"

MATRIX="${MATRIX:-test}"                  # test, train, or both
PERMUTATIONS="${PERMUTATIONS:-100000}"
METHOD="${METHOD:-exact}"                 # auto, exact, or monte_carlo
MAX_EXACT_STATES="${MAX_EXACT_STATES:-2000000}"
SEED="${SEED:-42}"
OUTPUT_TSV="${OUTPUT_TSV:-${INPUT_DIR}/permutation.tsv}"
PDF="${PDF:-*confusion_matrices*pdf}"
SIDECAR="${SIDECAR:-*confusion_matrices*csv}"

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
  echo "  METHOD=$METHOD"
  echo "  MAX_EXACT_STATES=$MAX_EXACT_STATES"
  echo "  SEED=$SEED"
  echo "  OUTPUT_TSV=$OUTPUT_TSV"
  echo "  PDF=$PDF"
  echo "  SIDECAR=$SIDECAR"
  echo
  echo "Example:"
  echo "  sbatch --export=ALL,INPUT_DIR=/path/to/results,METADATA_CATEGORIES=species,PERMUTATIONS=100000 permutation.sh"
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
  --method "$METHOD"
  --max_exact_states "$MAX_EXACT_STATES"
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

"$PYTHON" "$PROJECT_DIR/permutation.py" "${args[@]}" "$@"
