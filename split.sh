#!/bin/bash
#SBATCH --partition=horence
#SBATCH --time=0-04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=chesteryu@stanford.edu

set -euo pipefail

ml purge
eval "$(/oak/stanford/groups/horence/chester/dabs_ref/miniforge3/bin/conda shell.bash hook)" 
eval "$(mamba shell hook --shell bash)"
mamba activate biopython_env-R

PROJECT_DIR="${PROJECT_DIR:-/scratch/users/jiamuyu/proj_botryllus/flash}"
RESULTS_DIR="${RESULTS_DIR:-${PROJECT_DIR}/results/260826-01-2flies-wolbachia/filter1/shiftDist-levFilter/hyena/normalized/}"
DATASET_TABLE="${DATASET_TABLE:-${PROJECT_DIR}/dataset_table.csv}"
METADATA_FILE="${METADATA_FILE:-/scratch/users/jiamuyu/proj_botryllus/splash2/260826_00_2flies-wolbachia/metadata.csv}"

METADATA_CATEGORIES="${METADATA_CATEGORIES:-infectant,infection_status}"
SPLIT_METADATA_COL="${SPLIT_METADATA_COL:-sra_study}"
NUM_HITS="${NUM_HITS:-10}"

PLOT_SCRIPT="${PLOT_SCRIPT:-${PROJECT_DIR}/src/annotation/blast_code/plot_blast_annotations_each_feature_split_species.R}"
FILE_GLOB="${FILE_GLOB:-*summary_compactor.tsv}"
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
  echo "  OUTPUT_DIR=$RESULTS_DIR"
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
  --temp_dir "$RESULTS_DIR"
  --plot_script "$PLOT_SCRIPT"
  --rscript Rscript
  --num_hits "$NUM_HITS"
  --output_dir "$RESULTS_DIR"
  --work_dir "$WORK_DIR"
  --file_glob "$FILE_GLOB"
)

if [[ "$DRY_RUN" == "1" ]]; then
  args+=(--dry_run)
fi
if [[ "$PRODUCTS" == "1" ]]; then
  args+=(--products)
fi

python "$PROJECT_DIR/split.py" "${args[@]}" "$@"
