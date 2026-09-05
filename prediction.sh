#!/bin/bash
#
#SBATCH --partition=horence
#SBATCH --time=0-02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=chesteryu@stanford.edu

set -euo pipefail

ml purge
FLASH_CONDA="${FLASH_CONDA:-/oak/stanford/groups/horence/chester/dabs_ref/miniforge3}"
eval "$("$FLASH_CONDA/bin/conda" shell.bash hook)"
conda activate "${FLASH_ADELIE_ENV:-adelie_env}"

PROJECT_DIR="${PROJECT_DIR:-/scratch/users/jiamuyu/proj_botryllus/flash}"

# Edit these defaults or override them through sbatch --export.
DATASET="${DATASET:-260903-01-bschlosseri-age-a10t31}"
RESULTS_DIR="${RESULTS_DIR:-${PROJECT_DIR}/results/${DATASET}/filter1/noCluster/hyena/genomes/normalized}"
SIDECAR="${SIDECAR:-${RESULTS_DIR}/${DATASET}_hyena_adelie_genomes_results_top2000_k41_s41_trainProp0.8_confusion_matrices.tsv}"
METADATA_CATEGORY="${METADATA_CATEGORY:-age_in_days}"
MATRIX="${MATRIX:-both}"                 # train, test, or both
PLOT_TYPE="${PLOT_TYPE:-auto}"           # auto, regression, or confusion
OUTPUT="${OUTPUT:-${RESULTS_DIR}/${DATASET}_${METADATA_CATEGORY}_${MATRIX}_prediction_plots.pdf}"

# Regression point colors. Leave COLOR_COLUMN empty for one fixed point color.
COLOR_COLUMN="${COLOR_COLUMN-electric_shock}"
DOT_COLORS="${DOT_COLORS-no:#0072B2,yes:#D55E00}"
POINT_COLOR="${POINT_COLOR:-#2F6F9F}"
MISSING_COLOR="${MISSING_COLOR:-#A6A6A6}"
REGRESSION_CMAP="${REGRESSION_CMAP:-viridis}"

# Confusion-matrix class order and axis-label colors.
CLASS_ORDER="${CLASS_ORDER-no,yes}"
CLASS_COLORS="${CLASS_COLORS-no:#7F7F7F,yes:#D55E00}"
CONFUSION_CMAP="${CONFUSION_CMAP:-viridis}"

# Figure and text controls.
WIDTH="${WIDTH:-7.2}"
HEIGHT="${HEIGHT:-6.6}"
POINT_SIZE="${POINT_SIZE:-42}"
POINT_ALPHA="${POINT_ALPHA:-0.82}"
CELL_FONT_SIZE="${CELL_FONT_SIZE:-12}"
AXIS_TITLE_SIZE="${AXIS_TITLE_SIZE:-13}"
AXIS_TEXT_SIZE="${AXIS_TEXT_SIZE:-11}"
TITLE_SIZE="${TITLE_SIZE:-14}"
OTHER_TEXT_SIZE="${OTHER_TEXT_SIZE:-10}"
LEGEND_TEXT_SIZE="${LEGEND_TEXT_SIZE:-9}"

usage() {
  echo "Usage: sbatch prediction.sh [extra prediction.py args...]"
  echo
  echo "Hardcoded defaults:"
  echo "  PROJECT_DIR=$PROJECT_DIR"
  echo "  SIDECAR=$SIDECAR"
  echo "  METADATA_CATEGORY=$METADATA_CATEGORY"
  echo "  MATRIX=$MATRIX"
  echo "  PLOT_TYPE=$PLOT_TYPE"
  echo "  OUTPUT=$OUTPUT"
  echo "  COLOR_COLUMN=$COLOR_COLUMN"
  echo "  DOT_COLORS=$DOT_COLORS"
  echo "  CLASS_ORDER=$CLASS_ORDER"
  echo "  CLASS_COLORS=$CLASS_COLORS"
  echo "  CELL_FONT_SIZE=$CELL_FONT_SIZE"
  echo "  AXIS_TITLE_SIZE=$AXIS_TITLE_SIZE"
  echo "  AXIS_TEXT_SIZE=$AXIS_TEXT_SIZE"
  echo "  TITLE_SIZE=$TITLE_SIZE"
  echo "  OTHER_TEXT_SIZE=$OTHER_TEXT_SIZE"
  echo
  echo "Examples:"
  echo "  sbatch --export=ALL,METADATA_CATEGORY=age_in_days,MATRIX=test,COLOR_COLUMN=electric_shock prediction.sh"
  echo "  sbatch --export=ALL,METADATA_CATEGORY=infection_status,MATRIX=both,PLOT_TYPE=confusion,CLASS_ORDER=uninfected\,infected prediction.sh"
}

[[ "${1:-}" =~ ^(-h|--help)$ ]] && { usage; exit 0; }

python "$PROJECT_DIR/prediction.py" \
  --input "$SIDECAR" \
  --metadata-category "$METADATA_CATEGORY" \
  --matrix "$MATRIX" \
  --plot-type "$PLOT_TYPE" \
  --output "$OUTPUT" \
  --color-column "$COLOR_COLUMN" \
  --dot-colors "$DOT_COLORS" \
  --point-color "$POINT_COLOR" \
  --missing-color "$MISSING_COLOR" \
  --regression-cmap "$REGRESSION_CMAP" \
  --class-order "$CLASS_ORDER" \
  --class-colors "$CLASS_COLORS" \
  --confusion-cmap "$CONFUSION_CMAP" \
  --width "$WIDTH" \
  --height "$HEIGHT" \
  --point-size "$POINT_SIZE" \
  --point-alpha "$POINT_ALPHA" \
  --cell-font-size "$CELL_FONT_SIZE" \
  --axis-title-font-size "$AXIS_TITLE_SIZE" \
  --axis-text-font-size "$AXIS_TEXT_SIZE" \
  --title-font-size "$TITLE_SIZE" \
  --other-font-size "$OTHER_TEXT_SIZE" \
  --legend-font-size "$LEGEND_TEXT_SIZE" \
  "$@"
