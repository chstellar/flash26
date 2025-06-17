#!/bin/bash
NUM_ARRAY=50
export BIOPYTHON="/oak/stanford/groups/horence/dcotter1/envs/biopython/bin/activate"

INPUT_FASTA=$1
SPLIT_TEMP_FOLDER=$2
BLAST_OUTPUT_FOLDER=$3
OUTPUT_FILE=$4
THREADS=$5
TAXID=$6
TRANSLATION_TABLE=$7

if [[ -z $TAXID ]] ; then
  TAXID=0
fi

ml R/4.3.2
ml python/3.9.0
source $BIOPYTHON

mkdir -p $SPLIT_TEMP_FOLDER
mkdir -p $BLAST_OUTPUT_FOLDER

ml biology && ml ncbi-blast+ && python src/annotation/blast_code/run_blastp.py --input $INPUT_FASTA --split_folder $SPLIT_TEMP_FOLDER --blast_folder $BLAST_OUTPUT_FOLDER --max_workers $THREADS --taxid $TAXID --translation_table $TRANSLATION_TABLE

Rscript --vanilla src/annotation/blast_code/merge_blastp_with_GO.R --blast_folder $BLAST_OUTPUT_FOLDER --output_file $OUTPUT_FILE --max_workers $THREADS
