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
eval "$("$FLASH_CONDA/bin/conda" shell.bash hook)"
conda activate adelie_env

PROJECT_DIR="${PROJECT_DIR:-/scratch/users/jiamuyu/proj_botryllus/flash/}"
PYTHON="${PYTHON:-python}"

# Override these with sbatch --export=ALL,INPUT_DIR=...,METADATA_CATEGORIES=...
INPUT_DIR="${INPUT_DIR:-${PROJECT_DIR}/results/260819-00-cfloridanus-fungus/filter1/noCluster/hyena/normalized}"
METADATA_CATEGORIES="${METADATA_CATEGORIES:-fungus_species,tissue}"
MATRIX="${MATRIX:-test}"                  # test, train, or both
PERMUTATIONS="${PERMUTATIONS:-100000}"
SEED="${SEED:-42}"
OUTPUT_TSV="${OUTPUT_TSV:-${INPUT_DIR}/permutation.tsv}"
# PDF="${PDF:-${INPUT_DIR}/*confusion_matrices.pdf}"
PDF="${PDF:-${INPUT_DIR}/260819-00-cfloridanus-fungus_hyena_adelie_results_top2000_target1_k41_s41_trainProp0.8_confusion_matrices.pdf}"
SIDECAR="${SIDECAR:-${INPUT_DIR}/260819-00-cfloridanus-fungus_hyena_adelie_results_top2000_target1_k41_s41_trainProp0.8_confusion_matrices.csv}"

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

python "$PROJECT_DIR/permutation.py" "${args[@]}" "$@"
