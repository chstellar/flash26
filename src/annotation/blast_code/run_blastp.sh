#!/bin/bash

INPUT_FASTA=$1
SPLIT_TEMP_FOLDER=$2
BLAST_OUTPUT_FOLDER=$3
OUTPUT_FILE=$4
THREADS=$5
TAXID=$6
TRANSLATION_TABLE=$7
LOCAL_BLAST_DB=$8
PROTEIN_DB=${9:-}
TOP_N_SEQUENCES_PER_CLUSTER=${10:-0}
BLAST_SELECTION_MODE=${11:-all}
COEFFICIENTS_FILE=${12:-}
NUM_PLOT_HITS=${13:-10}
SAMPLE_SEQUENCES=${14:-}
CLUSTER_LENGTH=${15:-0}
REBLAST_OUTPUT_FILE=${16:-}

# Backward-compatible argument parsing:
#   run_blastp.sh ... LOCAL_BLAST_DB TOP_N
#   run_blastp.sh ... LOCAL_BLAST_DB PROTEIN_DB TOP_N
if [[ -z ${10:-} && $PROTEIN_DB =~ ^[0-9]+$ ]] ; then
  TOP_N_SEQUENCES_PER_CLUSTER=$PROTEIN_DB
  PROTEIN_DB=""
fi
if [[ ${10:-} =~ ^(all|top_n_per_cluster|plot_selected|plot_selected_and_top)$ ]] ; then
  TOP_N_SEQUENCES_PER_CLUSTER=$PROTEIN_DB
  PROTEIN_DB=""
  BLAST_SELECTION_MODE=${10:-all}
  COEFFICIENTS_FILE=${11:-}
  NUM_PLOT_HITS=${12:-10}
  SAMPLE_SEQUENCES=${13:-}
  CLUSTER_LENGTH=${14:-0}
  REBLAST_OUTPUT_FILE=${15:-}
fi

if [[ -z $TAXID ]] ; then
  TAXID=0
fi

# if there is no local blast db provided, it will be passed as 0
if [[ $LOCAL_BLAST_DB == "0" ]] ; then
  LOCAL_BLAST_DB=""
else
  LOCAL_BLAST_DB="--local_blast_db $LOCAL_BLAST_DB"
fi

# if PROTEIN_DB is not empty, add the "--protein_db $PROTEIN_DB" flag, otherwise set it to an empty string
if [[ -z $PROTEIN_DB ]] ; then
  PROTEIN_DB_FLAG=""
else
  PROTEIN_DB_FLAG="--protein_db $PROTEIN_DB"
fi

mkdir -p $SPLIT_TEMP_FOLDER
mkdir -p $BLAST_OUTPUT_FOLDER

python src/annotation/blast_code/run_blastp.py \
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
  --translation_table $TRANSLATION_TABLE ${LOCAL_BLAST_DB} \
  $PROTEIN_DB_FLAG # flag will be provided as --protein_db "/path/to/db" or will be empty # flag will be provided as --local_blast_db "/path/to/db" or will be empty
 

Rscript --vanilla src/annotation/blast_code/merge_blastp_with_GO.R \
  --blast_folder $BLAST_OUTPUT_FOLDER \
  --output_file $OUTPUT_FILE \
  --max_workers $THREADS

if [[ -n "$REBLAST_OUTPUT_FILE" ]] ; then
  REBLAST_QUERY_FASTA="${SPLIT_TEMP_FOLDER}_reblast_query.fasta"
  REBLAST_SPLIT_FOLDER="${SPLIT_TEMP_FOLDER}_reblast"
  REBLAST_BLAST_FOLDER="${BLAST_OUTPUT_FOLDER}_reblast"

  python src/annotation/blast_code/select_reblast_queries.py \
    --query_folder "$SPLIT_TEMP_FOLDER" \
    --annotations "$OUTPUT_FILE" \
    --output_fasta "$REBLAST_QUERY_FASTA" \
    --mode blastp

  if [[ -s "$REBLAST_QUERY_FASTA" ]] ; then
    mkdir -p "$REBLAST_SPLIT_FOLDER"
    mkdir -p "$REBLAST_BLAST_FOLDER"
    python src/annotation/blast_code/run_blastp.py \
      --input "$REBLAST_QUERY_FASTA" \
      --split_folder "$REBLAST_SPLIT_FOLDER" \
      --blast_folder "$REBLAST_BLAST_FOLDER" \
      --max_workers "$THREADS" \
      --taxid "0" \
      --translation_table "$TRANSLATION_TABLE" ${LOCAL_BLAST_DB} \
      $PROTEIN_DB_FLAG

    Rscript --vanilla src/annotation/blast_code/merge_blastp_with_GO.R \
      --blast_folder "$REBLAST_BLAST_FOLDER" \
      --output_file "$REBLAST_OUTPUT_FILE" \
      --max_workers "$THREADS"
  else
    printf "query\tidentity\tevalue\tqcovs\tqframe\tstaxids\tstitle\tNCBI_protein_accession\tUniProt_accession\tmethod\tGO\n" > "$REBLAST_OUTPUT_FILE"
  fi

  rm -f "$REBLAST_QUERY_FASTA"
  rm -rf "$REBLAST_SPLIT_FOLDER"
  rm -rf "$REBLAST_BLAST_FOLDER"
fi

rm -r $SPLIT_TEMP_FOLDER
rm -r $BLAST_OUTPUT_FOLDER
