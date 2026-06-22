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

snakemake --unlock -s 260620_00
snakemake --sdm conda --use-conda --conda-base-path /oak/stanford/groups/horence/chester/dabs_ref/miniforge3 --profile slurm_profile/ all_embeddings -s 260620_00

# snakemake --unlock -s Snakefile26062100
# snakemake --sdm conda --use-conda --conda-base-path /oak/stanford/groups/horence/chester/dabs_ref/miniforge3 --profile slurm_profile/ all_embeddings -s Snakefile26062100

# NUM_CORES=${SLURM_CPUS_PER_TASK:-1}
# snakemake --sdm conda --use-conda --conda-base-path /oak/stanford/groups/horence/chester/dabs_ref/miniforge3 --profile slurm_profile/ -j $NUM_CORES all_embeddings