#!/bin/bash

source ~/.bashrc
mamba activate snakemake
snakemake --rerun-incomplete --profile slurm_profile/ all_ohe

