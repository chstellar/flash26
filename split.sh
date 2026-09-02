#!/bin/bash
#SBATCH --partition=horence
#SBATCH --time=0-01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=chesteryu@stanford.edu

set -euo pipefail

ml purge
eval "$(/oak/stanford/groups/horence/chester/dabs_ref/miniforge3/bin/conda shell.bash hook)"
eval "$(mamba shell hook --shell bash)"
mamba activate default-R_env

PROJECT_DIR="/scratch/users/jiamuyu/proj_botryllus/flash"
RESULTS_DIR="${PROJECT_DIR}/results/260826-01-2flies-wolbachia/filter1/shiftDist-levFilter/hyena/normalized"
METADATA_FILE="/scratch/users/jiamuyu/proj_botryllus/splash2/260826_00_2flies-wolbachia/metadata.csv"
RUN_PREFIX="${RESULTS_DIR}/260826-01-2flies-wolbachia_hyena_adelie_results_top20000_target1_k54_s54_trainProp0.8"

METADATA_CATEGORIES="infectant,infection_status"
SPLIT_METADATA_COL="sra_study"

cd "$PROJECT_DIR"

Rscript --vanilla split.R \
  --nonzero_annotations "${RUN_PREFIX}_nonzero_coefficients_blastp_annotated_compactor.tsv" \
  --clusters "${PROJECT_DIR}/results/260826-01-2flies-wolbachia/filter1/shiftDist-levFilter/260826-01-2flies-wolbachia_sequences_per_cluster_top20000-clusters_target1_k54_s54.tsv" \
  --feather_file "${PROJECT_DIR}/results/260826-01-2flies-wolbachia/260826-01-2flies-wolbachia_hyena_top_variance_features_for_glmnet_filter1_shiftDist-levFilter_top20000_target1_k54_s54_normalized.feather" \
  --sample_seqs "${PROJECT_DIR}/results/260826-01-2flies-wolbachia/260826-01-2flies-wolbachia_prepared_sequences_filter1_shiftDist-levFilter_top20000_target1_k54_s54_sample_sequences.tsv" \
  --metadata "$METADATA_FILE" \
  --compactor_summary "${RUN_PREFIX}_nonzero_coefficients_blast_annotated_plots_summary_compactor.tsv" \
  --output "${RUN_PREFIX}_nonzero_coefficients_blast_annotated_plots_compactor_split-by-${SPLIT_METADATA_COL}.pdf" \
  --num_hits 10 \
  --cluster_length 54 \
  --metadata_categories "$METADATA_CATEGORIES" \
  --split_metadata_col "$SPLIT_METADATA_COL"
