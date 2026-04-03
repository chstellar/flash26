#!/bin/bash

# This script generates example data for the metaSPLASH pipeline.
# It downloads a sample dataset, processes it, and prepares the necessary files for testing the pipeline.
# Download the sample dataset from accessions in the resources/metadata directory.
# The dataset is a small subset of the H5N1 dataset from Nguyen et al. 2024.

# you must install fasterq-dump and prefetch from SRA Toolkit to run this script
# and have it in your PATH
# If you don't have SRA Toolkit installed, you can download it from 
# https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software

# Define the dataset name and the directory to store the example data
DATASET_NAME="H5N1-cattle"
EXAMPLE_DATA_DIR="example_data"
mkdir -p $EXAMPLE_DATA_DIR

# Download the sample dataset
accessions=$(cat resources/metadata/H5N1_accessions.txt)

$# Process each accession
echo "Downloading and processing accessions for dataset: $DATASET_NAME"
echo "Accessions to process: $accessions"
# Loop through each accession and download the data
echo "Processing accessions..."
echo "This may take a while, please be patient."
echo "Using prefetch and fasterq-dump to download and convert SRA files to FASTQ format."
echo "Make sure you have SRA Toolkit installed and in your PATH."

for accession in $accessions; do
    echo "Processing accession: $accession"
    prefetch $accession -O $EXAMPLE_DATA_DIR
    fasterq-dump --split-files --outdir $EXAMPLE_DATA_DIR $EXAMPLE_DATA_DIR/$accession/$accession.sra
    rm -r $EXAMPLE_DATA_DIR/$accession
done
echo "All accessions processed."

# create a sample sheet from the processed FASTQ files
echo "Creating sample sheet for dataset: $DATASET_NAME"
SAMPLE_SHEET="$EXAMPLE_DATA_DIR/sample_sheet.txt"
awk 'BEGIN {FS=OFS="\t"} {print $1,$1"_1.fastq"}' resources/metadata/H5N1_accessions.txt > $SAMPLE_SHEET

echo "Sample sheet created at: $SAMPLE_SHEET"

# to run SPLASH on the example, use the run_splash.sh script
echo "Example data generation complete."
