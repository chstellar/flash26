#!/usr/bin/env bash
#SBATCH --job-name=extendor_to_compactor
#SBATCH --output=extendor_to_compactor_%j.out
#SBATCH --error=extendor_to_compactor_%j.err
#SBATCH --time=00:10:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=horence,owners,normal

set -euo pipefail

module load devel 2>/dev/null || true
module load math 2>/dev/null || true

CONDA_SH="/oak/stanford/groups/horence/chester/dabs_ref/miniforge3/etc/profile.d/conda.sh"
if [[ -f "${CONDA_SH}" ]]; then
  source "${CONDA_SH}"
  conda activate /oak/stanford/groups/horence/chester/dabs_ref/miniforge3/envs/biopython_env
fi

PYTHON="${PYTHON:-python}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

BASE_DIR="/scratch/users/jiamuyu/proj_botryllus/flash/results/260714-00-3ants-challenge/filter1/noCluster/hyena/normalized"
PREFIX="260714-00-3ants-challenge_hyena_adelie_results_top2000_target1_k41_s41_trainProp0.8"

INPUT="${BASE_DIR}/test.txt"
OUTPUT="${BASE_DIR}/compactor.txt"
SEED_ANNOTATIONS="${BASE_DIR}/${PREFIX}_compactor_fungus_regular_seed_annotations.tsv"
SELECTED="${BASE_DIR}/${PREFIX}_compactor_fungus_regular_selected.tsv"

"${PYTHON}" "${SCRIPT_DIR}/extendor_to_compactor.py" \
  --input "${INPUT}" \
  --output "${OUTPUT}" \
  --seed_annotations "${SEED_ANNOTATIONS}" \
  --selected "${SELECTED}" \
  --anchor_len 10 \
  --target_len 31
