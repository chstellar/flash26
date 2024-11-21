#!/bin/bash
#
#SBATCH --job-name=metaSPLASH_snakemake
#SBATCH -p horence,normal
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --time=24:00:00

source /home/users/dcotter1/.bashrc
mamba activate snakemake
snakemake --profile slurm_profile/ -R prepare_data_for_glmnet_ohe all_ohe