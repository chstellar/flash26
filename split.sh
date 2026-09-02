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

METADATA_CATEGORIES="infectant,infection_status"
SPLIT_METADATA_COL="sra_study"

cd "$PROJECT_DIR"

python split.py \
  --results_dir "$RESULTS_DIR" \
  --metadata_file "$METADATA_FILE" \
  --metadata_categories "$METADATA_CATEGORIES" \
  --split_metadata_col "$SPLIT_METADATA_COL"
