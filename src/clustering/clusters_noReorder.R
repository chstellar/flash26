# clusters_noReorder.R
# Daniel Cotter
# 2024-09-12

# This script takes in an anchor cluster file and a splash stats file and 
# outputs the anchor cluster file with one anchor per cluster (the first one in the file)
# and keeps the order of the anchor clusters

## import packages --------
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))

## parse arguments --------
# define command line arguments
option_list <- list(
  make_option(c("-i", "--input_anchor_clusters"), "Input anchor clusters", type="character"),
  make_option(c("-s", "--splash_stats"), "SPLASH stats file", type="character"),
  make_option(c("-o", "--output"), "Output file", type="character"),
  make_option(c("--temp_dir"), "Temporary directory to store intermediate files", 
              type="character")
)

# parse command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# check that input file exists
if (!file.exists(opt$input_anchor_clusters) | !file.exists(opt$splash_stats) | is.null(opt$output)) {
  stop("Must provide input and output files")
}

# create a temporary directory to store intermediate files
if (!is.null(opt$temp_dir)) {
  temp_dir <- ifelse(grepl("/$", opt$temp_dir),
                     opt$temp_dir, 
                     paste0(opt$temp_dir, "/"))
  system(paste("mkdir -p", temp_dir))
} else {
  temp_dir <- file.path(dirname(opt$output_prefix), "tmp/")
  system(paste("mkdir -p", temp_dir))
}

# cat to screen the input files, output files, temp dir and paramaters
cat("\n###################################################################\n")
cat("Running reorder_anchor_clusters.R with the following parameters:\n")
cat(paste("Input clusters:", opt$input_anchor_clusters))
cat("\n")
cat(paste("SPLASH stats file:", opt$splash_stats))
cat("\n")
cat(paste("Output file:", opt$output))
cat("\n")
cat(paste("Temporary directory:", temp_dir))
cat("\n")
cat("###################################################################\n\n")

# read in the anchor cluster file
cat("Reading in the anchor clusters file\n")
anchor_clusters <- fread(opt$input_anchor_clusters, 
                         header = F, col.names = c("cluster_id", "anchor"))

cat("Writing out the new anchor clusters file\n")
write_tsv(anchor_clusters, opt$output, col_names = F, quote="none")

cat("Done!\n")


