#!/bin/bash
INPUT_FASTA=$1
SPLIT_TEMP_FOLDER=$2
BLAST_OUTPUT_FOLDER=$3
OUTPUT_FILE=$4
THREADS=$5
TAXID=$6
ENTREZ_EMAIL=$7
TEMP_DIR=$8

if [[ -z $TAXID ]] ; then
  TAXID=0
fi

mkdir -p $SPLIT_TEMP_FOLDER
mkdir -p $BLAST_OUTPUT_FOLDER

python src/annotation/blast_code/run_blast.py --input $INPUT_FASTA --split_folder $SPLIT_TEMP_FOLDER --blast_folder $BLAST_OUTPUT_FOLDER --max_workers $THREADS --taxid $TAXID

python src/annotation/blast_code/blast_features.py --blast_folder $BLAST_OUTPUT_FOLDER --output_file $OUTPUT_FILE --max_workers $THREADS --entrez_email $ENTREZ_EMAIL --temp_dir $TEMP_DIR

rm -r $SPLIT_TEMP_FOLDER
rm -r $BLAST_OUTPUT_FOLDER