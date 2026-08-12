#!/bin/bash
#
#SBATCH -p horence
#SBATCH --time=0-00:04:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=chesteryu@stanford.edu

set -euo pipefail

ml purge
eval "$(/oak/stanford/groups/horence/chester/dabs_ref/miniforge3/bin/conda shell.bash hook)"
conda activate biopython_env

BASE_DIR="/scratch/users/jiamuyu/proj_botryllus/flash/results/260714-00-3ants-challenge/filter1/noCluster/hyena/normalized"
PREFIX="260714-00-3ants-challenge_hyena_adelie_results_top2000_target1_k41_s41_trainProp0.8"

INPUT="${BASE_DIR}/test.txt"
OUTPUT="${BASE_DIR}/compactor.txt"
SEED_ANNOTATIONS="${BASE_DIR}/${PREFIX}_compactor_fungus_regular_seed_annotations.tsv"
SELECTED="${BASE_DIR}/${PREFIX}_compactor_fungus_regular_selected.tsv"

python extendor_to_compactor.py \
  --input "${INPUT}" \
  --output "${OUTPUT}" \
  --seed_annotations "${SEED_ANNOTATIONS}" \
  --selected "${SELECTED}" \
  --anchor_len 10 \
  --target_len 31