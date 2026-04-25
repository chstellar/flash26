#!/bin/bash

# This script is called by Snakemake

PYTHON_SCRIPT=$(realpath $1)
SINGULARITY_IMG=$(realpath $2)
INPUT_FILE=$(realpath $3)
OUTPUT_FILE=$(realpath $4)

# Internal paths (now baked into the container)
export MODEL_CFG="/opt/models/bacterial_128dim_config.yml"
export MODEL_CKPT="/opt/models/weights.ckpt"

singularity exec --nv ${SINGULARITY_IMG} \
    python ${PYTHON_SCRIPT} \
        --model_cfg ${MODEL_CFG} \
        --ckpt_path ${MODEL_CKPT} \
        --seq_file ${INPUT_FILE} \
        --output_file ${OUTPUT_FILE} \
        --max_seqlen 128000 \
        --nlayers 4 \
        --batch_size 100
