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

# First example: C. floridanus fungus, infection_status, clusters 8143 and 5883.
DATASET="${DATASET:-260819-00-cfloridanus-fungus-cp}"
RESULTS_DIR="${RESULTS_DIR:-${PROJECT_DIR}/results/${DATASET}/filter1/noCluster/hyena/normalized}"
METADATA="${METADATA:-/scratch/users/jiamuyu/proj_botryllus/splash2/260818_00_cfloridanus_fungus/metadata.csv}"
METADATA_COLUMN="${METADATA_COLUMN:-infection_status}"
CLUSTERS="${CLUSTERS:-8143,5883}"
CLUSTER_LABEL="${CLUSTERS//,/_}"
OUTPUT="${OUTPUT:-${RESULTS_DIR}/${DATASET}_${METADATA_COLUMN}_clusters_${CLUSTER_LABEL}_blast_plots.pdf}"
NUM_HITS="${NUM_HITS:-10}"
COMPACTOR_SUMMARY="${COMPACTOR_SUMMARY:-${RESULTS_DIR}/260819-00-cfloridanus-fungus_hyena_adelie_results_top2000_target1_k41_s41_trainProp0.8_nonzero_coefficients_blast_annotated_plots_summary_compactor.tsv}"
CLUSTERS_FILE="${CLUSTERS_FILE:-${PROJECT_DIR}/results/${DATASET}/260819-00-cfloridanus-fungus/filter1/noCluster/260819-00-cfloridanus-fungus_sequences_per_cluster_top2000-clusters_target1_k41_s41_annotated.tsv}"
FEATHER_FILE="${FEATHER_FILE:-${PROJECT_DIR}/results/${DATASET}/260819-00-cfloridanus-fungus_hyena_top_variance_features_for_glmnet_filter1_noCluster_top2000_target1_k41_s41_normalized.feather}"
SAMPLE_SEQS="${SAMPLE_SEQS:-${PROJECT_DIR}/results/${DATASET}/260819-00-cfloridanus-fungus_prepared_sequences_filter1_noCluster_top2000_target1_k41_s41_sample_sequences.tsv}"
CLUSTER_LENGTH="${CLUSTER_LENGTH:-41}"

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
  echo "  CLUSTERS_FILE=$CLUSTERS_FILE"
  echo "  FEATHER_FILE=$FEATHER_FILE"
  echo "  SAMPLE_SEQS=$SAMPLE_SEQS"
  echo "  CLUSTER_LENGTH=$CLUSTER_LENGTH"
  echo
  echo "Common overrides:"
  echo "  sbatch --export=ALL,METADATA_COLUMN=fungus_species,CLUSTERS=12\\,34 blot.sh"
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
  --clusters_file "$CLUSTERS_FILE" \
  --feather_file "$FEATHER_FILE" \
  --sample_seqs "$SAMPLE_SEQS" \
  --cluster_length "$CLUSTER_LENGTH" \
  --num_hits "$NUM_HITS" \
  "$@"
