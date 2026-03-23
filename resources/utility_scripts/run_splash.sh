#!/bin/bash
# This script runs SPLASH on the files specified in sample_sheet.txt
# The sample_sheet.txt file should have two columns: sample_id and input_fastq
# The output SPLASH files will be written to the current directory
# Do not change --outname_prefix as output files will be names result.* and these 
# are the default names used by FLASH.
# Modify the threads and membory paramaters as needed. Increasing n_bins will also
# improved memory usage at the cost of speed.

SPLASH_BIN="PATH/TO/SPLASH-2.*/splash"  # Modify this path to point to your SPLASH binary 

$SPLASH_BIN --outname_prefix result \ 
    --anchor_len 27 \
    --gap_len 0 \
    --target_len 27 \
    --poly_ACGT_len 8 \
    --max_pval_opt_for_Cjs 0.05 \
    --with_effect_size_cts \
    --with_pval_asymp_opt \
    --n_threads_stage_1 16 \
    --n_threads_stage_2 16 \
    --n_bins 64 \
    --n_most_freq_targets 10 \
    --kmc_use_RAM_only_mode \
    --kmc_max_mem_GB 6 \
    --dump_sample_anchor_target_count_binary \
    --keep_top_effect_size_bin_anchors_satc \
    --keep_top_target_entropy_anchors_satc \
    --satc_merge_dump_format satc \
    sample_sheet.txt
