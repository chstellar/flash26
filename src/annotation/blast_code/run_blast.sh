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
TOP_N_SEQUENCES_PER_CLUSTER=${10:-0}
BLAST_SELECTION_MODE=${11:-all}
COEFFICIENTS_FILE=${12:-}
NUM_PLOT_HITS=${13:-10}
SAMPLE_SEQUENCES=${14:-}
CLUSTER_LENGTH=${15:-0}

if [[ -z $TAXID ]] ; then
  TAXID=0
fi

# if there is no local blast db provided, it will be passed as 0
if [[ $LOCAL_BLAST_DB == "0" ]] ; then
  LOCAL_BLAST_DB=""
else
  LOCAL_BLAST_DB="--local_blast_db $LOCAL_BLAST_DB"
fi

# if there is no entrez email provided, it will be passed as 0
# throw error if email not provided
if [[ $ENTREZ_EMAIL == "0" ]] ; then
  echo "Error: Entrez email must be set in config.yml to run BLAST."
  exit 1
fi

mkdir -p $SPLIT_TEMP_FOLDER
mkdir -p $BLAST_OUTPUT_FOLDER

python src/annotation/blast_code/run_blast.py \
  --input $INPUT_FASTA \
  --split_folder $SPLIT_TEMP_FOLDER \
  --blast_folder $BLAST_OUTPUT_FOLDER \
  --max_workers $THREADS \
  --taxid "$TAXID" \
  --top_n_sequences_per_cluster "$TOP_N_SEQUENCES_PER_CLUSTER" \
  --blast_selection_mode "$BLAST_SELECTION_MODE" \
  --coefficients "$COEFFICIENTS_FILE" \
  --num_plot_hits "$NUM_PLOT_HITS" \
  --sample_sequences "$SAMPLE_SEQUENCES" \
  --cluster_length "$CLUSTER_LENGTH" \
  $LOCAL_BLAST_DB # flag will be provided as --local_blast_db "/path/to/db" or will be empty

python src/annotation/blast_code/blast_features.py \
  --blast_folder $BLAST_OUTPUT_FOLDER \
  --output_file $OUTPUT_FILE \
  --entrez_email $ENTREZ_EMAIL \
  --temp_dir $TEMP_DIR

rm -r $SPLIT_TEMP_FOLDER
rm -r $BLAST_OUTPUT_FOLDER
