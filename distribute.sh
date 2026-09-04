#!/bin/bash
#
#SBATCH --partition=horence
#SBATCH --time=0-01:00:00
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
  mamba activate "${FLASH_PY_ENV:-biopython_env-R}"
elif command -v conda >/dev/null 2>&1; then
  conda activate "${FLASH_PY_ENV:-biopython_env-R}"
fi

PYTHON="${PYTHON:-python}"

# PROJECT_DIR="${PROJECT_DIR:-/scratch/users/jiamuyu/proj_botryllus/flash}"
# RESULTS_DIR="${PROJECT_DIR}/results/260826-01-2flies-wolbachia/filter1/shiftDist-levFilter/hyena/normalized"
# PARTITION_SHEET="${PARTITION_SHEET:-/scratch/users/jiamuyu/proj_botryllus/splash2/260826_00_2flies-wolbachia/partition.species.csv}"
# # cut -d',' -f1,2,7 metadata.csv > partition.csv
# SATC=$(ls ${RESULTS_DIR}/../../target1/*clusters/all_satc_merged.txt 2>/dev/null | head -1)
# INPUT_TSV="${INPUT_TSV:-${RESULTS_DIR}/infectant.tsv}"
# # grep "fungus_species" ${RESULTS_DIR}/*summary_compactor.tsv > $INPUT_TSV
# grep "Nedd8" ${RESULTS_DIR}/*summary_compactor.tsv | grep -v "residual" > $INPUT_TSV
# OUTPUT_TSV="${OUTPUT_TSV:-${RESULTS_DIR}/distribution.tsv}"
# HEATMAP_PDF="${HEATMAP_PDF:-${RESULTS_DIR}/distribution_heatmaps.pdf}"


PROJECT_DIR="${PROJECT_DIR:-/scratch/users/jiamuyu/proj_botryllus/flash}"
RESULTS_DIR="${PROJECT_DIR}/results/260819-00-cfloridanus-fungus-cp/filter1/noCluster/hyena/normalized"
PARTITION_SHEET="${PARTITION_SHEET:-/scratch/users/jiamuyu/proj_botryllus/splash2/260818_00_cfloridanus_fungus/partition.csv}"
# cut -d',' -f1,2,7 metadata.csv > partition.csv
SATC=$(ls ${RESULTS_DIR}/../../target1/*clusters/all_satc_merged.txt 2>/dev/null | head -1)
INPUT_TSV="${INPUT_TSV:-${RESULTS_DIR}/manuscript.tsv}"
# grep "fungus_species" ${RESULTS_DIR}/*summary_compactor.tsv > $INPUT_TSV
OUTPUT_TSV="${OUTPUT_TSV:-${RESULTS_DIR}/distribution.tsv}"
HEATMAP_PDF="${HEATMAP_PDF:-${RESULTS_DIR}/distribution_heatmaps.pdf}"

# Column settings. These can be header names or 1-based column indices.
# EXTENDOR_COL="${EXTENDOR_COL:-4}"
# ANNOTATION_COL="${ANNOTATION_COL:-20}"
# SAMPLE_COL="${SAMPLE_COL:-1}"
# MAJOR_COL="${MAJOR_COL:-3}"
# MINOR_COL="${MINOR_COL:-2}"
# MAJOR_ORDER="${MAJOR_ORDER:-no_wolbachia,reduced_wolbachia,wolbachia}"
# MINOR_ORDER="${MINOR_ORDER:-SRP179717,SRP243262,SRP295874,SRP423200}"
# MAJOR_DISPLAY_NAME="${MAJOR_DISPLAY_NAME:-wolbachia_infection}"
# MINOR_DISPLAY_NAME="${MINOR_DISPLAY_NAME:-sra_study}"
# MAJOR_COLORS="${MAJOR_COLORS:-no_wolbachia:#7F7F7F,reduced_wolbachia:#9467BD,wolbachia:#0072B2}"
EXTENDOR_COL="${EXTENDOR_COL:-4}"
ANNOTATION_COL="${ANNOTATION_COL:-20}"
SAMPLE_COL="${SAMPLE_COL:-1}"
MAJOR_COL="${MAJOR_COL:-2}"
MINOR_COL="${MINOR_COL:-3}"
MAJOR_ORDER="${MAJOR_ORDER:-no_fungus,bbassiana,ocamponoti-floridani}"
MINOR_ORDER="${MINOR_ORDER:-manipulation,timecourse}"
MAJOR_DISPLAY_NAME="${MAJOR_DISPLAY_NAME:-fungus_species}"
MINOR_DISPLAY_NAME="${MINOR_DISPLAY_NAME:-dataset}"
MAJOR_COLORS="${MAJOR_COLORS:-no_fungus:#7F7F7F,not_given:#9467BD,bbassiana:#0072B2,ocamponoti-floridani:#D55E00}"
# no_fungus:#7F7F7F,not_given:#9467BD,bbassiana:#0072B2,ocamponoti-floridani:#D55E00

INPUT_HAS_HEADER="${INPUT_HAS_HEADER:-auto}"
PARTITION_HAS_HEADER="${PARTITION_HAS_HEADER:-auto}"
PARTITION_DELIMITER="${PARTITION_DELIMITER:-auto}"
SATC_HAS_HEADER="${SATC_HAS_HEADER:-auto}"
EXTENDOR_ORDER="${EXTENDOR_ORDER:-anchor-target}"
LONG_OUTPUT="${LONG_OUTPUT:-}"
HEATMAP_BASE_FONT_SIZE="${HEATMAP_BASE_FONT_SIZE:-14}"
HEATMAP_ENTRY_FONT_SIZE="${HEATMAP_ENTRY_FONT_SIZE:-18}"
HEATMAP_TOTAL_FONT_SIZE="${HEATMAP_TOTAL_FONT_SIZE:-14}"

usage() {
  echo "Usage: sbatch distribute.sh [extra distribute.py args...]"
  echo
  echo "This wrapper is configured at the top of distribute.sh."
  echo "Default input/output:"
  echo "  INPUT_TSV=$INPUT_TSV"
  echo "  PARTITION_SHEET=$PARTITION_SHEET"
  echo "  OUTPUT_TSV=$OUTPUT_TSV"
  echo "  HEATMAP_PDF=$HEATMAP_PDF"
  echo "  SATC=$SATC"
  echo
  echo "Common optional environment variables:"
  echo "  EXTENDOR_COL    extendor column name or 1-based index; default 4"
  echo "  ANNOTATION_COL  annotation column name or 1-based index; default 20"
  echo "  SAMPLE_COL      sample column in partition sheet; default 1"
  echo "  MAJOR_COL       major partition column; default 2"
  echo "  MINOR_COL       minor partition column; default 3"
  echo "  MAJOR_ORDER     comma-separated row order"
  echo "  MINOR_ORDER     comma-separated column order"
  echo "  MAJOR_COLORS    comma-separated row-label color map, label:#RRGGBB"
  echo "  HEATMAP_BASE_FONT_SIZE   label/title baseline font size; default 14"
  echo "  HEATMAP_ENTRY_FONT_SIZE  numbers inside main heatmap cells; default 18"
  echo "  HEATMAP_TOTAL_FONT_SIZE  numbers inside total cells; default 14"
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
  --extendor_col "$EXTENDOR_COL"
  --annotation_col "$ANNOTATION_COL"
  --sample_col "$SAMPLE_COL"
  --major_col "$MAJOR_COL"
  --minor_col "$MINOR_COL"
  --major_order "$MAJOR_ORDER"
  --minor_order "$MINOR_ORDER"
  --major_display_name "$MAJOR_DISPLAY_NAME"
  --minor_display_name "$MINOR_DISPLAY_NAME"
  --major_colors "$MAJOR_COLORS"
  --heatmap_base_font_size "$HEATMAP_BASE_FONT_SIZE"
  --heatmap_entry_font_size "$HEATMAP_ENTRY_FONT_SIZE"
  --heatmap_total_font_size "$HEATMAP_TOTAL_FONT_SIZE"
  --input_has_header "$INPUT_HAS_HEADER"
  --partition_has_header "$PARTITION_HAS_HEADER"
  --partition_delimiter "$PARTITION_DELIMITER"
  --satc_has_header "$SATC_HAS_HEADER"
  --extendor_order "$EXTENDOR_ORDER"
)

if [[ -n "$LONG_OUTPUT" ]]; then
  args+=(--long_output "$LONG_OUTPUT")
fi
if [[ -n "$HEATMAP_PDF" ]]; then
  args+=(--heatmap_pdf "$HEATMAP_PDF")
fi

"$PYTHON" "$PROJECT_DIR/distribute.py" "${args[@]}" "$@"
