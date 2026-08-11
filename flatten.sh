#!/bin/bash

set -euo pipefail

INPUT_DIR=/scratch/users/jiamuyu/proj_botryllus/flash/results
INPUT_PAT="260720-00-3ants-challenge" # can contain trailing `*` to capture multiple dirs at once
OUTPUT_DIR=/scratch/groups/horence/chester/flash2share

# derive subdirectory name from INPUT_PAT with the trailing * removed
OUTPUT_PAT="${INPUT_PAT%\*}"
OUTPUT_SUBDIR="${OUTPUT_DIR}/${OUTPUT_PAT}"

mkdir -p "$OUTPUT_SUBDIR"
> "${OUTPUT_SUBDIR}/raw_paths.txt"
cd "$INPUT_DIR"

shopt -s nullglob

for dir in $INPUT_PAT; do
    [ -d "$dir" ] || continue
    for pat in compactor.pdf heatmaps.pdf matrices.pdf summary_compactor.tsv unannotated.tsv; do
        for f in "$dir"/filter1/*/hyena/normalized/*${pat}; do
            cp "$f" "$OUTPUT_SUBDIR/"
            echo "$f" >> "${OUTPUT_SUBDIR}/raw_paths.txt"
        done
    done
done

echo "saved key flattened results from ${INPUT_DIR}/${INPUT_PAT} to ${OUTPUT_SUBDIR}"