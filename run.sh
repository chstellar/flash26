#!/bin/bash
#
#SBATCH -p horence
#SBATCH --time=0-12:00:00 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=chesteryu@stanford.edu

eval "$(/oak/stanford/groups/horence/chester/dabs_ref/miniforge3/bin/conda shell.bash hook)" 
eval "$(mamba shell hook --shell bash)"
mamba activate flash

NUM_CORES=${SLURM_CPUS_PER_TASK:-1}

snakemake --sdm conda -j $NUM_CORES all_embeddings
# snakemake all_embeddings