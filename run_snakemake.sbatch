#!/bin/bash

source /home/users/dcotter1/.bashrc
mamba activate snakemake
snakemake --rerun-incomplete --resources msa=1 --profile slurm_profile/ -R plot_blast_features -R plot_blast_features_OHE -R plot_blast_heatmaps -R plot_blast_heatmaps_OHE -R merge_annotations -R merge_annotations_OHE all all_ohe

