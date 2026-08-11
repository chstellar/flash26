#!/bin/bash

set -euo pipefail

INPUT_DIR=/scratch/users/jiamuyu/proj_botryllus/flash/results
INPUT_PAT="${1:-260720-00-3ants-challenge*}"  # may contain trailing * to capture multiple dirs
OUTPUT_DIR=/scratch/groups/horence/chester/flash2share
OPTIONAL_SUFFIX="${2:-}"                        # default empty (safe under set -u)

OUTPUT_PAT="${INPUT_PAT%\*}"
OUTPUT_PAT="${OUTPUT_PAT//-/_}"
OUTPUT_PAT="${OUTPUT_PAT}${OPTIONAL_SUFFIX:+_${OPTIONAL_SUFFIX}}"
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

ml system rclone/1.73.1
rclone copy "${OUTPUT_SUBDIR}" "gdrive:${OUTPUT_PAT}" -P