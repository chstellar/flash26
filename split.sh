#!/bin/bash
#
# Temporary BLAST replot launcher. Parameters are intentionally kept here.
#
#SBATCH --job-name=split_blast_replot
#SBATCH --partition=horence
#SBATCH --time=0-04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
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
RSCRIPT="${RSCRIPT:-Rscript}"

# Main tunables.
PROJECT_DIR="${PROJECT_DIR:-/scratch/users/jiamuyu/proj_botryllus/flash}"
RESULTS_DIR="${RESULTS_DIR:-${PROJECT_DIR}/results/260826-01-2flies-wolbachia}"
TEMP_DIR="${TEMP_DIR:-${PROJECT_DIR}/results}"
DATASET_TABLE="${DATASET_TABLE:-${PROJECT_DIR}/dataset_table.csv}"
METADATA_FILE="${METADATA_FILE:-/scratch/users/jiamuyu/proj_botryllus/splash2/260826_00_2flies-wolbachia/metadata.csv}"

METADATA_CATEGORIES="${METADATA_CATEGORIES:-infectant,infection_status}"
SPLIT_METADATA_COL="${SPLIT_METADATA_COL:-sra_study}"
NUM_HITS="${NUM_HITS:-10}"

PLOT_SCRIPT="${PLOT_SCRIPT:-${PROJECT_DIR}/src/annotation/blast_code/plot_blast_annotations_each_feature_split_species.R}"
OUTPUT_DIR="${OUTPUT_DIR:-${RESULTS_DIR}/split_replots}"
WORK_DIR="${WORK_DIR:-${OUTPUT_DIR}/filtered_inputs}"
FILE_GLOB="${FILE_GLOB:-**/*_nonzero_coefficients_blastp_annotated.tsv}"
DRY_RUN="${DRY_RUN:-0}"
PRODUCTS="${PRODUCTS:-0}"

usage() {
  echo "Usage: sbatch split.sh [extra split.py args...]"
  echo
  echo "Edit parameters near the top of split.sh or override them with --export."
  echo "Defaults:"
  echo "  RESULTS_DIR=$RESULTS_DIR"
  echo "  METADATA_CATEGORIES=$METADATA_CATEGORIES"
  echo "  SPLIT_METADATA_COL=$SPLIT_METADATA_COL"
  echo "  METADATA_FILE=$METADATA_FILE"
  echo "  OUTPUT_DIR=$OUTPUT_DIR"
  echo
  echo "Example:"
  echo "  sbatch --export=ALL,SPLIT_METADATA_COL=sra_study,NUM_HITS=20 split.sh"
}

[[ "${1:-}" =~ ^(-h|--help)$ ]] && { usage; exit 0; }

cd "$PROJECT_DIR"

args=(
  --results_dir "$RESULTS_DIR"
  --metadata_categories "$METADATA_CATEGORIES"
  --split_metadata_col "$SPLIT_METADATA_COL"
  --metadata_file "$METADATA_FILE"
  --dataset_table "$DATASET_TABLE"
  --temp_dir "$TEMP_DIR"
  --plot_script "$PLOT_SCRIPT"
  --rscript "$RSCRIPT"
  --num_hits "$NUM_HITS"
  --output_dir "$OUTPUT_DIR"
  --work_dir "$WORK_DIR"
  --file_glob "$FILE_GLOB"
)

if [[ "$DRY_RUN" == "1" ]]; then
  args+=(--dry_run)
fi
if [[ "$PRODUCTS" == "1" ]]; then
  args+=(--products)
fi

"$PYTHON" "$PROJECT_DIR/split.py" "${args[@]}" "$@"
