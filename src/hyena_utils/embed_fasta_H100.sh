#!/bin/bash

# This script is called by Snakemake

# Arguments: dev_embedder.py input_file output_file
PYTHON_SCRIPT=$1
SINGULARITY_IMG=$2
INPUT_FILE=$3
OUTPUT_FILE=$4

# Internal paths (now baked into the container)
export MODEL_CFG="/wdr/models/bacterial_128dim_config.yml"
export MODEL_CKPT="/wdr/models/weights.ckpt"

export PYTHON_SCRIPT=$(realpath $1)
export SINGULARITY_IMG=$(realpath $2)
export INPUT_FILE=$(realpath $3)
export OUTPUT_FILE=$(realpath $4)


# Run the Python script with the singularity image
singularity run --nv --writable-tmpfs \
    ${SINGULARITY_IMG} bash -c 'cd /wdr/hyena_wdr && python ${PYTHON_SCRIPT} \
    --model_cfg ${MODEL_CFG} --ckpt_path ${MODEL_CKPT} \
    --seq_file ${INPUT_FILE} --output_file ${OUTPUT_FILE} \
    --max_seqlen 128000 --nlayers 4 --batch_size 100'