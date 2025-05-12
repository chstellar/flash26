#!/bin/bash
NUM_ARRAY=50
export IMAGE="/oak/stanford/groups/horence/julias/splash_postprocessing/khoa_annotation/SPLASH_UTILS/feature_annotation/python-sequtils-blast-splash-pfam_latest.sif"
export REPO="/oak/stanford/groups/horence/khoa/scratch/repos/SPLASH_UTILS"
export BLAST_DB="/scratch/users/dcotter1/blast_db/"
export ENV="singularity run -B $REPO,$OAK,$BLAST_DB $IMAGE"

INPUT_FASTA=$1
SPLIT_TEMP_FOLDER=$2
BLAST_OUTPUT_FOLDER=$3
OUTPUT_FILE=$4
THREADS=$5

mkdir -p $SPLIT_TEMP_FOLDER
mkdir -p $BLAST_OUTPUT_FOLDER

$ENV python src/annotation/blast_code/run_blast.py --input $INPUT_FASTA --split_folder $SPLIT_TEMP_FOLDER --blast_folder $BLAST_OUTPUT_FOLDER --max_workers $THREADS

$ENV python src/annotation/blast_code/blast_features.py --blast_folder $BLAST_OUTPUT_FOLDER --output_file $OUTPUT_FILE --max_workers $THREADS
