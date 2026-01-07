#!/bin/bash
INPUT_FASTA=$1
SPLIT_TEMP_FOLDER=$2
BLAST_OUTPUT_FOLDER=$3
OUTPUT_FILE=$4
THREADS=$5
TAXID=$6
ENTREZ_EMAIL=$7
TEMP_DIR=$8
LOCAL_BLAST_DB=$9

if [[ -z $TAXID ]] ; then
  TAXID=0
fi

# if there is no local blast db provided, it will be passed as 0
if [[ $LOCAL_BLAST_DB == "0" ]] ; then
  LOCAL_BLAST_DB=""
else
  LOCAL_BLAST_DB="--local_blast_db $LOCAL_BLAST_DB"
fi

mkdir -p $SPLIT_TEMP_FOLDER
mkdir -p $BLAST_OUTPUT_FOLDER

python src/annotation/blast_code/run_blast.py \
  --input $INPUT_FASTA \
  --split_folder $SPLIT_TEMP_FOLDER \
  --blast_folder $BLAST_OUTPUT_FOLDER \
  --max_workers $THREADS \
  --taxid $TAXID \
  $LOCAL_BLAST_DB # flag will be provided as --local_blast_db "/path/to/db" or will be empty

python src/annotation/blast_code/blast_features.py \
  --blast_folder $BLAST_OUTPUT_FOLDER \
  --output_file $OUTPUT_FILE \
  --entrez_email $ENTREZ_EMAIL \
  --temp_dir $TEMP_DIR

rm -r $SPLIT_TEMP_FOLDER
rm -r $BLAST_OUTPUT_FOLDER