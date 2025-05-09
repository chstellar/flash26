#!/bin/bash
#
#SBATCH --job-name=metaSPLASH_snakemake
#SBATCH -p horence
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --time=72:00:00

source /home/users/dcotter1/.bashrc
mamba activate snakemake
snakemake -R plot_blast_features --profile slurm_profile/ all


ml system
ml rclone
rclone copy results/ gdrive:Salzman\ Lab/Shared\ Materials/MetaSPLASH_all_results/results_by_dataset/may_08_2025

