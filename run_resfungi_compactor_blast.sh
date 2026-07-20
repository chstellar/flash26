#!/bin/bash
#
#SBATCH -p horence
#SBATCH --time=0-4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=chesteryu@stanford.edu

set -euo pipefail

# Run from the FLASH repository root even if sbatch is submitted elsewhere.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# The BLAST helper scripts call `python` and `Rscript` internally, so make sure
# those resolve to modern module-provided executables.
if command -v ml >/dev/null 2>&1; then
  ml purge || true
  ml python/3.9.0
  ml R/4.4.2 || ml R/4.3.2
else
  module --force purge || true
  module load python/3.9.0
  module load R/4.4.2 || module load R/4.3.2
fi

CONDA_BASE="${CONDA_BASE:-/oak/stanford/groups/horence/chester/dabs_ref/miniforge3}"
CONDA_ENV="${CONDA_ENV:-/oak/stanford/groups/horence/chester/dabs_ref/miniforge3/envs/biopython_env-R}"
if [ -s "$CONDA_BASE/bin/activate" ] && [ -d "$CONDA_ENV" ]; then
  # The FLASH BLAST scripts need Biopython, pandas, and R packages. This is
  # the same environment used by the Snakemake BLAST rules.
  eval "$("$CONDA_BASE/bin/conda" shell.bash hook)"
  if command -v mamba >/dev/null 2>&1; then
    eval "$(mamba shell hook --shell bash)"
    mamba activate "$CONDA_ENV"
  else
    conda activate "$CONDA_ENV"
  fi
else
  echo "WARNING: Could not find conda env $CONDA_ENV; falling back to module python/R." >&2
fi

echo "Using python: $(command -v python)"
python --version
echo "Using python3: $(command -v python3)"
python3 --version
echo "Using Rscript: $(command -v Rscript)"
Rscript --version

python - <<'PY'
import Bio
import pandas
print("Python dependencies OK: Bio, pandas")
PY

Rscript -e 'library(data.table); cat("R dependencies OK: data.table\n")'

# FLASH helper wrappers call bare `python`; on Sherlock that can still be
# /bin/python 2.7. Put a tiny shim first in PATH so wrappers get Python 3.
PYTHON_SHIM_DIR="${TMPDIR:-/tmp}/resfungi_python_shim_${USER:-user}_$$"
mkdir -p "$PYTHON_SHIM_DIR"
ln -sf "$(command -v python)" "$PYTHON_SHIM_DIR/python"
ln -sf "$(command -v python)" "$PYTHON_SHIM_DIR/python3"
export PATH="$PYTHON_SHIM_DIR:$PATH"
trap 'rm -rf "$PYTHON_SHIM_DIR"' EXIT

echo "After shim, python: $(command -v python)"
python --version

export ENTREZ_EMAIL="${ENTREZ_EMAIL:-v8514616@outlook.com}"

if [ -n "${RESFUNGI_THREADS:-}" ]; then
  THREADS="$RESFUNGI_THREADS"
elif [ -n "${SLURM_JOB_ID:-}" ]; then
  THREADS="${SLURM_CPUS_PER_TASK:-32}"
else
  THREADS=32
fi
echo "Using BLAST threads: $THREADS"
echo "Using BLAST modes: ${RESFUNGI_BLAST_MODES:-both}"
echo "Using reblast mode: ${RESFUNGI_REBLAST_MODE:-missing}"
echo "Using protein BLAST database: ${RESFUNGI_PROTEIN_DB:-refseq_protein}"

python /scratch/users/jiamuyu/proj_botryllus/flash/resfungi_compactor_blast.py \
  --threads "$THREADS" \
  --taxids "300111;102681;104421" \
  --translation_table 1 \
  --entrez_email "$ENTREZ_EMAIL" \
  --blast_modes "${RESFUNGI_BLAST_MODES:-both}" \
  --reblast_mode "${RESFUNGI_REBLAST_MODE:-missing}" \
  --protein_db "${RESFUNGI_PROTEIN_DB:-refseq_protein}" \
  "$@"
