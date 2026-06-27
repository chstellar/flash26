#!/bin/bash
#
#SBATCH -p horence
#SBATCH --time=2-00:00:00 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=chesteryu@stanford.edu

ml purge
eval "$(/oak/stanford/groups/horence/chester/dabs_ref/miniforge3/bin/conda shell.bash hook)" 
eval "$(mamba shell hook --shell bash)"
mamba activate flash

SNAKEMAKE_FILE="${1:-260620_00}"
FORCE_RULE="${2:-}"

snakemake --unlock -s $SNAKEMAKE_FILE

# Build the snakemake command
SNAKEMAKE_CMD="snakemake --sdm conda --use-conda --conda-base-path /oak/stanford/groups/horence/chester/dabs_ref/miniforge3 --profile slurm_profile/"

# Add forcerun flag if a rule is specified
if [ -n "$FORCE_RULE" ]; then
    SNAKEMAKE_CMD="$SNAKEMAKE_CMD -R $FORCE_RULE"
fi

SNAKEMAKE_CMD="$SNAKEMAKE_CMD all_embeddings -s $SNAKEMAKE_FILE"

# Execute the command
eval $SNAKEMAKE_CMD

# snakemake --unlock -s $SNAKEMAKE_FILE
# snakemake --sdm conda --use-conda --conda-base-path /oak/stanford/groups/horence/chester/dabs_ref/miniforge3 --profile slurm_profile/ all_embeddings -s $SNAKEMAKE_FILE

# NUM_CORES=${SLURM_CPUS_PER_TASK:-1}
# snakemake --sdm conda --use-conda --conda-base-path /oak/stanford/groups/horence/chester/dabs_ref/miniforge3 --profile slurm_profile/ -j $NUM_CORES all_embeddings