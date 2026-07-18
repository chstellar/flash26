#!/bin/bash
#
#SBATCH -p horence
#SBATCH --time=0-12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
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

echo "Using python: $(command -v python)"
python --version
echo "Using python3: $(command -v python3)"
python3 --version
echo "Using Rscript: $(command -v Rscript)"
Rscript --version

# FLASH helper wrappers call bare `python`; on Sherlock that can still be
# /bin/python 2.7. Put a tiny shim first in PATH so wrappers get Python 3.
PYTHON_SHIM_DIR="${TMPDIR:-/tmp}/resfungi_python_shim_${USER:-user}_$$"
mkdir -p "$PYTHON_SHIM_DIR"
ln -sf "$(command -v python3)" "$PYTHON_SHIM_DIR/python"
ln -sf "$(command -v python3)" "$PYTHON_SHIM_DIR/python3"
export PATH="$PYTHON_SHIM_DIR:$PATH"
trap 'rm -rf "$PYTHON_SHIM_DIR"' EXIT

echo "After shim, python: $(command -v python)"
python --version

export ENTREZ_EMAIL="${ENTREZ_EMAIL:-v8514616@outlook.com}"

python3 resfungi_compactor_blast.py \
  --threads "${SLURM_CPUS_PER_TASK:-32}" \
  --taxids "300111;102681;104421" \
  --translation_table 1 \
  --entrez_email "$ENTREZ_EMAIL" \
  "$@"
