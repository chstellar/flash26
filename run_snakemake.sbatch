#!/bin/bash
#
#SBATCH --job-name=metaSPLASH_snakemakes
#SBATCH -p horence,normal
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --time=24:00:00

mamba activate snakemake
snakemake --profile slurm_profile/