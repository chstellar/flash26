#!/bin/bash

ml python/3.9.0

INPUT_TSV="results/260621-00-temnothorax-challenge/filter1/noCluster/hyena/normalized/260621-00-temnothorax-challenge_hyena_adelie_results_top20000_target1_k40_s40_trainProp0.8_nonzero_coefficients_blast_annotated.tsv"
OUTPUT_TSV="test.tsv"

if [ -z "$INPUT_TSV" ] || [ -z "$OUTPUT_TSV" ]; then
    echo "Usage: $0 <input_tsv> <output_tsv>"
    exit 1
fi

if [ ! -f "$INPUT_TSV" ]; then
    echo "Error: Input file '$INPUT_TSV' not found"
    exit 1
fi

# Process the TSV file
cut -f 1,2,15 "$INPUT_TSV" | awk 'BEGIN {FS=OFS="\t"} 
NR==1 {
    # Handle header
    print $1, $2, "gene", "note"
    next
}
{
    # Truncate column 2 (keep first two underscore-separated parts)
    split($2, a, "_")
    col2 = a[1] "_" a[2]
    
    # Extract gene field
    gene = ""
    if (match($3, /'\''gene'\'': \['\''([^'\'']+)'\''\]/, g)) {
        gene = g[1]
    }
    
    # Extract note field
    note = ""
    if (match($3, /'\''note'\'': \['\''([^'\'']+)'\''\]/, n)) {
        note = n[1]
    }
    
    print $1, col2, gene, note
}' | (head -n 1; tail -n +2 | sort -t$'\t' -k1,1 -k3,3 | uniq) > "$OUTPUT_TSV"

echo "Processing complete: $OUTPUT_TSV"