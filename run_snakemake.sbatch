#!/bin/bash
#
#SBATCH --job-name=metaSPLASH_snakemake
#SBATCH -p horence
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --time=72:00:00

source /home/users/dcotter1/.bashrc
mamba activate snakemake
snakemake --rerun-incomplete --resources msa=1 --profile slurm_profile/ -R plot_blast_features -R plot_blast_features_OHE -R plot_blast_heatmaps -R plot_blast_heatmaps_OHE -R merge_annotations -R merge_annotations_OHE all all_ohe


ml system
ml rclone
rclone copy results/ gdrive:Salzman\ Lab/Shared\ Materials/MetaSPLASH_all_results/results_by_dataset/may_09_2025

