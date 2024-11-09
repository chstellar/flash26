# select_anchors_filter1.R
# Daniel Cotter
# 2024-09-12

# This script selects the most important anchors from the output of SPLASH
# it prioritizes anchors that have a high effect size and also have a large 
# number of nonzero samples. It also filters out anchors that have lookup
# table hits to artifcats. The output is a list of anchor sequences that can
# be used in downstream analyses.

# Filter 1: effect size >= 0.6
# Filter 1: number of nonzero samples > 10th percentile
# Filter 1: no lookup table hits to artifacts
# Filter 1: select top 150,000 anchors by number of nonzero samples

## import packages --------
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(Biostrings))

## parse arguments --------
# define command line arguments
option_list <- list(
  make_option(c("-i", "--input"), "Input file", type="character"),
  make_option(c("-o", "--output"), "Output file", type="character"),
  make_option(c("-n", "--num_anchors"), "Number of anchors to select",
              type="integer", default = 150000),
  make_option(c("-e", "--effect_size"), "Effect size threshold",
              type="numeric", default=0.6),
  make_option(c("-l", "--lookup_table"), "Lookup table file", type="character"),
  make_option(c("--splash_bin"), "Path to SPLASH binary folder",
              type="character", default="/oak/stanford/groups/horence/dcotter1/splash-2.6.1/"),
  make_option(c("--temp_dir"), "Temporary directory to store intermediate files", 
              type="character")
)

# parse command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# check that input file exists
if (!file.exists(opt$input) | is.null(opt$output)) {
  stop("Must provide input and output files")
}
if (!file.exists(opt$lookup_table)) {
  stop("Must provide path to lookup table")
}

# define output files 
anchors_only_out = opt$output

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
cat("Running select_important_anchors.R with the following parameters:\n")
cat(paste("Input file:", opt$input))
cat("\n")
cat(paste("Output anchors:", anchors_only_out))
cat("\n")
cat(paste("Number of anchors to select:", opt$num_anchors))
cat("\n")
cat(paste("Effect size threshold:", opt$effect_size))
cat("\n")
cat(paste("Lookup table file:", opt$lookup_table))
cat("\n")
cat(paste("SPLASH binary folder:", opt$splash_bin))
cat("\n")
cat(paste("Temporary directory:", temp_dir))
cat("\n")
cat("###################################################################\n\n")


## load the data
# read in the headers of the input file to identify the effect_size_bin column
headers <- fread(opt$input, nrows = 1, header=T)
effect_size_bin_col <- grep("effect_size_bin", names(headers))

# load the input file using awk to filter out rows with effect size < 0.6
load_cmd <- paste0("cat ", opt$input, " | awk '{OFS=\"\t\"}{if ($", effect_size_bin_col, " >= ", opt$effect_size, ") print $0}'")
if (grepl(".gz$", opt$input)) {
  load_cmd = paste0("z", load_cmd)
}

# only select the first 18 columns
cat(paste("Reading in anchors from", opt$input, "with effect size >=", opt$effect_size))
cat("\n\n")
dt <- fread(cmd=load_cmd, header=TRUE, select = 1:18)

# # filter out the bottom 10% of anchors by number nonzero samples
# sample_cutoff <- quantile(dt$number_nonzero_samples, 0.1)
# dt <- dt %>% filter(number_nonzero_samples > sample_cutoff)

## select the most important anchors -----------
anchors_to_keep <- dt %>% select(anchor, effect_size_bin, number_nonzero_samples) %>%
  arrange(desc(effect_size_bin), desc(number_nonzero_samples)) %>%
  head(opt$num_anchors)

## write the output -----------
cat(paste("Writing", nrow(anchors_to_keep), "anchors to", anchors_only_out))
cat("\n")
anchors_to_keep$anchor %>% 
  write.table(anchors_only_out, row.names=FALSE, col.names=FALSE, quote=FALSE)
