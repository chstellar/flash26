#!/bin/bash
#
#SBATCH --partition=horence
#SBATCH --time=0-02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
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
  mamba activate "${FLASH_R_ENV:-default-R_env}"
elif command -v conda >/dev/null 2>&1; then
  conda activate "${FLASH_R_ENV:-default-R_env}"
fi

PROJECT_DIR="${PROJECT_DIR:-/scratch/users/jiamuyu/proj_botryllus/flash}"

# Current manuscript subset example.
DATASET="${DATASET:-260819-00-cfloridanus-fungus-cp}"
RESULTS_DIR="${RESULTS_DIR:-${PROJECT_DIR}/results/${DATASET}/filter1/noCluster/hyena/normalized}"
METADATA="${METADATA:-/scratch/users/jiamuyu/proj_botryllus/splash2/260818_00_cfloridanus_fungus/metadata.csv}"
METADATA_COLUMN="${METADATA_COLUMN:-fungus_species}"
CLUSTERS="${CLUSTERS:-1512,1483}"
CLUSTER_LABEL="${CLUSTERS//,/_}"
OUTPUT="${OUTPUT:-${RESULTS_DIR}/${DATASET}_${METADATA_COLUMN}_clusters_${CLUSTER_LABEL}_blast_plots.pdf}"
NUM_HITS="${NUM_HITS:-10}"
COMPACTOR_SUMMARY="${COMPACTOR_SUMMARY:-${RESULTS_DIR}/260819-00-cfloridanus-fungus_hyena_adelie_results_top2000_target1_k41_s41_trainProp0.8_nonzero_coefficients_blast_annotated_plots_summary_compactor.tsv}"
GENE_LABEL_SIZE="${GENE_LABEL_SIZE:-5.5}"
HIST_LABEL_SCALE="${HIST_LABEL_SCALE:-1.25}"
AXIS_TITLE_SIZE="${AXIS_TITLE_SIZE:-16}"
AXIS_TEXT_SIZE="${AXIS_TEXT_SIZE:-12}"
PLOT_TITLE_SIZE="${PLOT_TITLE_SIZE:-16}"
PLOT_SUBTITLE_SIZE="${PLOT_SUBTITLE_SIZE:-11}"
CLASS_ORDER="${CLASS_ORDER:-none,bbassiana,ocamponoti-floridani}"
CLASS_COLORS="${CLASS_COLORS:-none:#7F7F7F,bbassiana:#0072B2,ocamponoti-floridani:#D55E00}"

usage() {
  echo "Usage: sbatch blot.sh [extra blot.R args...]"
  echo
  echo "Hardcoded defaults:"
  echo "  PROJECT_DIR=$PROJECT_DIR"
  echo "  RESULTS_DIR=$RESULTS_DIR"
  echo "  METADATA=$METADATA"
  echo "  METADATA_COLUMN=$METADATA_COLUMN"
  echo "  CLUSTERS=$CLUSTERS"
  echo "  OUTPUT=$OUTPUT"
  echo "  COMPACTOR_SUMMARY=$COMPACTOR_SUMMARY"
  echo "  GENE_LABEL_SIZE=$GENE_LABEL_SIZE"
  echo "  HIST_LABEL_SCALE=$HIST_LABEL_SCALE"
  echo "  AXIS_TITLE_SIZE=$AXIS_TITLE_SIZE"
  echo "  AXIS_TEXT_SIZE=$AXIS_TEXT_SIZE"
  echo "  PLOT_TITLE_SIZE=$PLOT_TITLE_SIZE"
  echo "  PLOT_SUBTITLE_SIZE=$PLOT_SUBTITLE_SIZE"
  echo "  CLASS_ORDER=$CLASS_ORDER"
  echo "  CLASS_COLORS=$CLASS_COLORS"
  echo
  echo "Common overrides:"
  echo "  sbatch --export=ALL,METADATA_COLUMN=fungus_species,CLUSTERS=12\\,34,GENE_LABEL_SIZE=5,AXIS_TITLE_SIZE=18 blot.sh"
}

[[ "${1:-}" =~ ^(-h|--help)$ ]] && { usage; exit 0; }

Rscript --vanilla "$PROJECT_DIR/blot.R" \
  --project_dir "$PROJECT_DIR" \
  --results_dir "$RESULTS_DIR" \
  --metadata "$METADATA" \
  --metadata_column "$METADATA_COLUMN" \
  --clusters "$CLUSTERS" \
  --output "$OUTPUT" \
  --compactor_summary "$COMPACTOR_SUMMARY" \
  --num_hits "$NUM_HITS" \
  --label_size "$GENE_LABEL_SIZE" \
  --hist_label_scale "$HIST_LABEL_SCALE" \
  --axis_title_size "$AXIS_TITLE_SIZE" \
  --axis_text_size "$AXIS_TEXT_SIZE" \
  --plot_title_size "$PLOT_TITLE_SIZE" \
  --plot_subtitle_size "$PLOT_SUBTITLE_SIZE" \
  --class_order "$CLASS_ORDER" \
  --class_colors "$CLASS_COLORS" \
  "$@"
